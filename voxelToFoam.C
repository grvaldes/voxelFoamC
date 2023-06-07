/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    voxelToFoam

Description
    Reads .msh file as written by Gmsh.

    Needs surface elements on mesh to be present and aligned with outside faces
    of the mesh. I.e. if the mesh is hexes, the outside faces need to be
    quads.

    Note: There is something seriously wrong with the ordering written in the
    .msh file. Normal operation is to check the ordering and invert prisms
    and hexes if found to be wrong way round.
    Use the -keepOrientation to keep the raw information.

    Note: The code now uses the element (cell,face) physical region id number
    to create cell zones and faces zones (similar to
    fluentMeshWithInternalFaces).

    A use of the cell zone information, is for field initialization with the
    "setFields" utility. see the classes:  topoSetSource, zoneToCell.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "mergePolyMesh.H"
#include "IFstream.H"
#include "cellModeller.H"
#include "repatchPolyTopoChanger.H"
#include "cellSet.H"
#include "faceSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Element type sizes
static label SIZE_NOT_DEFINED = 0;
static label NODQUAD = 4;
static label FACEHEX = 6;
static label NODHEX  = 8;

enum axisDir { dirX, dirY, dirZ };


// Comma delimited row parsing.
List<string> headerParse(string line)
{
    List<string> row;

    label pos;
    label subpos;
    label rowI(0);

    while (true)
    {
        pos = line.find(",");
        row.append(line.substr(0, pos));

        if ((subpos = row[rowI].find("=")) != -1)
        {
            row[rowI] = row[rowI].substr(subpos+1, row[rowI].length() - subpos);
        }

        line.erase(0, pos+1);
        rowI++;

        if (pos == -1)
        {
            return row;
        }
    }
}

// Comma delimited row parsing.
List<scalar> scalarParse(string line, label ListSize)
{
    List<scalar> row(ListSize);

    label pos;
    label rowI(0);

    while (true)
    {
        pos = line.find(",");

        if (ListSize == SIZE_NOT_DEFINED)
        {
            row.append(stod(line.substr(0, pos)));
        }
        else
        {
            row[rowI] = stod(line.substr(0, pos));
        }

        line.erase(0, pos+1);
        rowI++;

        if (pos == -1)
        {
            return row;
        }
    }
}

// Comma delimited row parsing.
List<label> labelParse(string line, label ListSize)
{
    List<label> row(ListSize);

    label pos;
    label rowI(0);

    while (true)
    {
        pos = line.find(",");

        if (ListSize == SIZE_NOT_DEFINED)
        {
            row.append(stoi(line.substr(0, pos)));
        }
        else
        {
            row[rowI] = stoi(line.substr(0, pos));
        }

        line.erase(0, pos+1);
        rowI++;

        if (pos == -1)
        {
            return row;
        }
    }
}


// Set which cells are associated to a point
// Copy from Foam::labelListList Foam::polyMesh::cellShapePointCells
labelListList shapePointCells
(
    const cellShapeList& cells,
    const pointField& points
)
{
    List<DynamicList<label, primitiveMesh::cellsPerPoint_>>
        pc(points.size());

    forAll(cells, i)
    {
        const labelList& labels = cells[i];

        forAll(labels, j)
        {
            label curPoint = labels[j];
            DynamicList<label, primitiveMesh::cellsPerPoint_>& curPointCells =
                pc[curPoint];

            curPointCells.append(i);
        }
    }

    labelListList pointCellAddr(pc.size());

    forAll(pc, pointi)
    {
        pointCellAddr[pointi].transfer(pc[pointi]);
    }

    return pointCellAddr;
}


// Set which cells are associated to a point
// Copy from Foam::labelListList Foam::polyMesh::cellShapePointCells
labelListList shapePointFaces
(
    const faceList& faces,
    const pointField& points
)
{
    List<DynamicList<label, primitiveMesh::facesPerPoint_>>
        pf(points.size());

    forAll(faces, i)
    {
        const labelList& labels = faces[i];

        forAll(labels, j)
        {
            label curPoint = labels[j];
            DynamicList<label, primitiveMesh::facesPerPoint_>& curPointFaces =
                pf[curPoint];

            curPointFaces.append(i);
        }
    }

    labelListList pointFaceAddr(pf.size());

    forAll(pf, pointi)
    {
        pointFaceAddr[pointi].transfer(pf[pointi]);
    }

    return pointFaceAddr;
}


// Reads mesh format
void readFileHeading(IFstream& inFile, string& line)
{
    inFile.getLine(line);
    Info<< line << "\n" << endl;
}


// Reads points and map
void readPoints
(
    IFstream& inFile, 
    string& line, 
    pointField& points, 
    Map<label>& texgenToFoam
)
{
    label pointi = 0;

    Info<< "Starting to read points at line " 
    << inFile.lineNumber()
    << endl;

    while (inFile.good())
    { 
        inFile.getLine(line);

        if (line(1) == "*")
        {
            Info<< "Finished reading nodes. Vertices read: "
                << texgenToFoam.size() << "\n"
                << endl;
            break;
        }

        List<scalar> row = scalarParse(line, 4);

        scalar mshLabel = row[0];
        scalar xVal = row[1];
        scalar yVal = row[2]; 
        scalar zVal = row[3];

        point pt;

        pt.x() = xVal;
        pt.y() = yVal;
        pt.z() = zVal;

        points.append(pt);
        texgenToFoam.insert(mshLabel, pointi);

        pointi++;
    }
}


// Reads cells and patch faces
void readCells
(
    IFstream& inFile,
    string& line,
    const pointField& points,
    const Map<label>& texgenToFoam,
    cellShapeList& cellAsShapes,
    cellList& cells
)
{
    const cellModel& hex = *(cellModeller::lookup("hex"));

    face quadPoints(NODQUAD);
    labelList hexPoints(NODHEX);

    // cellShapeList cellAsShapes;

    Info<< "Starting to read cells at line " 
        << inFile.lineNumber()
        << endl;
    
    List<string> header = headerParse(line);

    string elmType = header[1];
    label celli = 0;

    while (inFile.good())
    {
        inFile.getLine(line);

        if (line(1) == "*")
        {
            Info<< "Finished reading cells. Elements read: "
                << cells.size() << "\n"
                << endl;
            break;
        }

        List<label> row = labelParse(line, NODHEX+1);

        forAll (hexPoints, i)
        {
            hexPoints[i] = row[i+1];
        }

        forAll(hexPoints, labelI)
        {
            hexPoints[labelI] = texgenToFoam[hexPoints[labelI]];
        }

        cellAsShapes.append( cellShape(hex, hexPoints) );
        cells.append( cell( FACEHEX ) );

        celli++;
    }

    if (cells.size() == 0 || cellAsShapes.size() == 0)
    {
        FatalIOErrorInFunction(inFile)
            << "No cells read from file " << inFile.name() << nl
            << exit(FatalIOError);
    }
}


// Set faces based on cellShapes.
// Copy from void Foam::polyMesh::setTopology
// Copy from void Foam::polyMesh::initMesh(cellList& c)
void setFaces
(
    const cellShapeList& cellsAsShapes,
    const pointField& points,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    cellList& cells
)
{
    label maxFaces = 0;
    faceListList cellsFaceShapes(cellsAsShapes.size());


    // Setting all faces of mesh from cells
    forAll(cellsFaceShapes, celli)
    {
        cellsFaceShapes[celli] = cellsAsShapes[celli].faces();
        maxFaces += cellsFaceShapes[celli].size();
    }

    faces.setSize(maxFaces);

    label nFaces = 0;
    label nInternalFaces = 0;

    boolList markedFaces(maxFaces, false);
    bool found = false;

    labelListList PointCells = shapePointCells(cellsAsShapes, points);


    // Numbering faces of mesh avoiding duplicates
    forAll(cells, celli)
    {
        const faceList& currrentFaces = cellsFaceShapes[celli];
        
        labelList neighbourCells(currrentFaces.size(), -1);
        labelList faceOfNeighbourCell(currrentFaces.size(), -1);

        label nNeighbours = 0;

        forAll(currrentFaces, facei)
        {
            if (cells[celli][facei] >= 0) continue;

            found = false;
            const face& curFace = currrentFaces[facei];
            const labelList& curPoints = curFace;

            forAll(curPoints, pointi)
            {
                const labelList& curNeighbours =
                    PointCells[curPoints[pointi]];

                forAll(curNeighbours, neiI)
                {
                    label curNei = curNeighbours[neiI];

                    if (curNei > celli)
                    {
                        const faceList& searchFaces = cellsFaceShapes[curNei];

                        forAll(searchFaces, neiFacei)
                        {
                            if (searchFaces[neiFacei] == curFace)
                            {
                                found = true;

                                neighbourCells[facei] = curNei;
                                faceOfNeighbourCell[facei] = neiFacei;
                                nNeighbours++;

                                break;
                            }
                        }
                        if (found) break;
                    }
                    if (found) break;
                }
                if (found) break;
            }
        }

        for (label neiSearch = 0; neiSearch < nNeighbours; neiSearch++)
        {
            label nextNei = -1;
            label minNei = cells.size();

            forAll(neighbourCells, ncI)
            {
                if (neighbourCells[ncI] > -1 && neighbourCells[ncI] < minNei)
                {
                    nextNei = ncI;
                    minNei = neighbourCells[ncI];
                }
            }

            if (nextNei > -1)
            {
                faces[nFaces] = currrentFaces[nextNei];
                cells[celli][nextNei] = nFaces;
                cells[neighbourCells[nextNei]][faceOfNeighbourCell[nextNei]] = nFaces;
                neighbourCells[nextNei] = -1;
                nFaces++;
            }
            else
            {
                FatalErrorInFunction
                    << "Error in internal face insertion"
                    << abort(FatalError);
            }
        }
    }


    forAll(cells, celli)
    {
        labelList& curCellFaces = cells[celli];

        forAll(curCellFaces, facei)
        {
            if (curCellFaces[facei] == -1)
            {
                curCellFaces[facei] = nFaces;
                faces[nFaces] = cellsFaceShapes[celli][facei];

                nFaces++;
            }
        }
    }

    // Setting owner-neighbour pair
    faces.setSize(nFaces);
    owner.setSize(nFaces);
    neighbour.setSize(nFaces);
    
    forAll(cells, celli)
    {
        const labelList& cellfaces = cells[celli];

        forAll(cellfaces, facei)
        {
            if (cellfaces[facei] < 0)
            {
                FatalErrorInFunction
                    << "Illegal face label " << cellfaces[facei]
                    << " in cell " << celli
                    << exit(FatalError);
            }

            if (!markedFaces[cellfaces[facei]])
            {
                owner[cellfaces[facei]] = celli;
                markedFaces[cellfaces[facei]] = true;
            }
            else
            {
                neighbour[cellfaces[facei]] = celli;
                nInternalFaces++;
            }
        }
    }

    neighbour.setSize(nInternalFaces);
}


// Reads physical names
void readElSet
(
    IFstream& inFile, 
    string& line,
    Map<word>& elementSets,
    List<DynamicList<label>>& zoneCells
)
{
    label zoneI = zoneCells.size();

    List<string> header = headerParse(line);
    string regionName = header[1];

    if (regionName(3) != "All")
    {
        zoneCells.setSize(zoneI+1);
        
        Info<< "Mapping Element Set " << regionName
            << " to Foam pointZone " << zoneI 
            << " from line " << inFile.lineNumber()
            << endl;

        elementSets.insert(zoneI, string::validate<word>(regionName));
        
        while (inFile.good())
        {
            inFile.getLine(line);

            if (line(1) == "*")
            {
                Info<< "Finished reading ElSet. Elements read: "
                    << zoneCells[zoneI].size() << "\n"
                    << endl;
                break;
            }

            List<label> row = labelParse(line, SIZE_NOT_DEFINED);

            forAll (row, i)
            {
                zoneCells[zoneI].append(row[i]-1);
            }
        }
    }
    else
    {
        inFile.getLine(line);
    }
}


// Reads physical names
void readNSet
(
    IFstream& inFile, 
    string& line,
    Map<word>& pointSets,
    const Map<label>& texgenToFoam,
    List<DynamicList<label>>& zonePoints
)
{
    label zoneI = zonePoints.size();

    List<string> header = headerParse(line);
    string regionName = header[1];

    if 
    (
        regionName(11) != "Constraints" && 
        regionName(3) != "All"
    )
    {
        zonePoints.setSize(zoneI+1);
        
        Info<< "Mapping Point Set " << regionName
            << " to Foam pointZone " << zoneI 
            << " from line " << inFile.lineNumber()
            << endl;

        pointSets.insert(zoneI, string::validate<word>(regionName));
        
        while (inFile.good())
        {
            inFile.getLine(line);

            if (line(1) == "*")
            {
                Info<< "Finished reading NSet. Nodes read: "
                    << zonePoints[zoneI].size() << "\n"
                    << endl;
                break;
            }

            List<label> row = labelParse(line, SIZE_NOT_DEFINED);
            
            forAll(row, labelI)
            {
                row[labelI] = texgenToFoam[row[labelI]];
            }

            forAll (row, i)
            {
                zonePoints[zoneI].append(row[i]);
            }
        }
    }
    else
    {
        inFile.getLine(line);
    }
}


// Sets boundary properties
// Connects the boundaries from texgen to how they will be used in OpenFOAM.
// Hardcoded and ugly, I will find a solution later.
// Other functions depend on having the order FACE/EDGE/VERTEX
void setBoundaryProperties
(
    labelList& boundaryPatchIndices,
    wordList& boundaryPatchNames,
    wordList& boundaryPatchType,
    labelListList& boundaryComponents,
    List<DynamicList<label>>& zonePoints
)
{
    boundaryPatchNames[0] = "Right";
    boundaryPatchNames[1] = "Left";
    boundaryPatchNames[2] = "Back";
    boundaryPatchNames[3] = "Front";
    boundaryPatchNames[4] = "Top";
    boundaryPatchNames[5] = "Bottom";

    forAll(boundaryPatchNames, patchi)
    {
        boundaryPatchIndices[patchi] = patchi;
        boundaryPatchType[patchi] = "patch";
    }

    forAll(boundaryComponents, bc)
    {
        switch (bc)
        {
            case 0:
                boundaryComponents[bc].append(zonePoints[0]);
                boundaryComponents[bc].append(zonePoints[7]);
                boundaryComponents[bc].append(zonePoints[8]);
                boundaryComponents[bc].append(zonePoints[11]);
                boundaryComponents[bc].append(zonePoints[12]);
                boundaryComponents[bc].append(zonePoints[19]);
                boundaryComponents[bc].append(zonePoints[20]);
                boundaryComponents[bc].append(zonePoints[23]);
                boundaryComponents[bc].append(zonePoints[24]);
                break;
            case 1:
                boundaryComponents[bc].append(zonePoints[1]);
                boundaryComponents[bc].append(zonePoints[9]);
                boundaryComponents[bc].append(zonePoints[6]);
                boundaryComponents[bc].append(zonePoints[10]);
                boundaryComponents[bc].append(zonePoints[13]);
                boundaryComponents[bc].append(zonePoints[18]);
                boundaryComponents[bc].append(zonePoints[21]);
                boundaryComponents[bc].append(zonePoints[22]);
                boundaryComponents[bc].append(zonePoints[25]);
                break;
            case 2:
                boundaryComponents[bc].append(zonePoints[2]);
                boundaryComponents[bc].append(zonePoints[15]);
                boundaryComponents[bc].append(zonePoints[16]);
                boundaryComponents[bc].append(zonePoints[8]);
                boundaryComponents[bc].append(zonePoints[9]);
                boundaryComponents[bc].append(zonePoints[20]);
                boundaryComponents[bc].append(zonePoints[21]);
                boundaryComponents[bc].append(zonePoints[24]);
                boundaryComponents[bc].append(zonePoints[25]);
                break;
            case 3:
                boundaryComponents[bc].append(zonePoints[3]);
                boundaryComponents[bc].append(zonePoints[14]);
                boundaryComponents[bc].append(zonePoints[17]);
                boundaryComponents[bc].append(zonePoints[6]);
                boundaryComponents[bc].append(zonePoints[7]);
                boundaryComponents[bc].append(zonePoints[18]);
                boundaryComponents[bc].append(zonePoints[19]);
                boundaryComponents[bc].append(zonePoints[22]);
                boundaryComponents[bc].append(zonePoints[23]);
                break;
            case 4:
                boundaryComponents[bc].append(zonePoints[4]);
                boundaryComponents[bc].append(zonePoints[16]);
                boundaryComponents[bc].append(zonePoints[17]);
                boundaryComponents[bc].append(zonePoints[12]);
                boundaryComponents[bc].append(zonePoints[13]);
                boundaryComponents[bc].append(zonePoints[22]);
                boundaryComponents[bc].append(zonePoints[23]);
                boundaryComponents[bc].append(zonePoints[24]);
                boundaryComponents[bc].append(zonePoints[25]);
                break;
            case 5:
                boundaryComponents[bc].append(zonePoints[5]);
                boundaryComponents[bc].append(zonePoints[14]);
                boundaryComponents[bc].append(zonePoints[15]);
                boundaryComponents[bc].append(zonePoints[10]);
                boundaryComponents[bc].append(zonePoints[11]);
                boundaryComponents[bc].append(zonePoints[18]);
                boundaryComponents[bc].append(zonePoints[19]);
                boundaryComponents[bc].append(zonePoints[20]);
                boundaryComponents[bc].append(zonePoints[21]);
                break;
        }
    }
}


// Setting boundary elements
void setBoundaryElems
(
    const cellShapeList& cellsAsShapes,
    const pointField& points,
    const labelListList& boundaryComponents,
    cellList& cells,
    faceList& faces,
    labelList& owner,
    faceListList& boundaryFaces
)
{
    labelListList PointFaces = shapePointFaces(faces, points);
    labelListList PointCells = shapePointCells(cellsAsShapes, points);

    label nTotalFaces = faces.size();
    label nNewFaces = 0;
    Map<label> reorderedFaces;
    Map<label> reverseReorderedFaces;
    Map<label> allReorderedFaces;

    // Building boundaries
    forAll(boundaryComponents, patchi)
    {
        const labelList& patchPoints = boundaryComponents[patchi];
        Map<label> countFaces;

        forAll(patchPoints, pointi)
        {
            labelList curFaces = PointFaces[patchPoints[pointi]];

            forAll(curFaces, facei)
            {
                if (!countFaces.found(curFaces[facei]))
                {
                    countFaces.insert(curFaces[facei], 1);
                }
                else
                {
                    countFaces[curFaces[facei]]++;
                }
            }
        }

        forAllIter(Map<label>, countFaces, iter)
        {
            if (iter() == NODQUAD)
            {
                boundaryFaces[patchi].append(faces[iter.key()]);
                reorderedFaces.insert(nNewFaces,iter.key());
                reverseReorderedFaces.insert(iter.key(), nNewFaces);
                nTotalFaces--;
                nNewFaces++;
            }
        }
    }

    faceList newFaces;
    labelList newOwner;
    label ind = 0;

    forAll(faces, facei)
    {
        if(!reverseReorderedFaces.found(facei))
        {
            newFaces.append(faces[facei]);
            newOwner.append(owner[facei]);
            allReorderedFaces.insert(ind, facei);
            ind++;
        }
    }

    forAllIter(Map<label>, reorderedFaces, iter)
    {
        newFaces.append(faces[iter()]);
        newOwner.append(owner[iter()]);
        allReorderedFaces.insert(ind, iter());
        ind++;
    }

    faces = newFaces;
    owner = newOwner;

    forAll(cells, celli)
    {
        cell& curCell = cells[celli];

        forAll(curCell, facei)
        {
            curCell[facei] = allReorderedFaces[curCell[facei]];
        }
    }
}


polyMesh createMesh
(
    pointField& points,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    const labelList& boundaryPatchIndices,
    const wordList& boundaryPatchNames,
    const wordList& boundaryTypes,
    const faceListList& boundaryFaces,
    const List<DynamicList<label>>& zoneCells,
    const Map<word>& elementSets,
    const Time& runTime
)
{  
    labelList patchSizes(FACEHEX);
    labelList patchStarts(FACEHEX);
    List<cellZone*> cellZones;
    List<polyPatch*> boundaryPatches;
    label nFace = neighbour.size();

    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        move(points),
        move(faces),
        move(owner),
        move(neighbour)
    );

    forAll(boundaryFaces, patchi)
    {
        patchStarts[patchi] = nFace;
        patchSizes[patchi] = boundaryFaces[patchi].size();

        boundaryPatches.append(
            new polyPatch
            (
                boundaryPatchNames[patchi],
                patchSizes[patchi],
                patchStarts[patchi],
                boundaryPatchIndices[patchi],
                mesh.boundaryMesh(),
                boundaryTypes[patchi]
            )
        );

        nFace += patchSizes[patchi];
    }

    mesh.addPatches(boundaryPatches);

    forAll(zoneCells, zonei)
    {
        cellZones.append
        (
            new cellZone
            (
                elementSets[zonei],
                zoneCells[zonei],
                zonei,
                mesh.cellZones()
            )
        );
    }

    mesh.addZones(List<pointZone*>(0), List<faceZone*>(0), cellZones);

    return mesh;
}













int main(int argc, char *argv[])
{
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // PROGRAM INITIALIZATION
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    argList::noParallel();
    argList::validArgs.append(".inp file");
    argList::addBoolOption("debug");
    argList::addOption
    (
        "repeat",
        "vector",
        "repeats the mesh in XYZ by the specified vector "
        "- eg '(2 2 2)'"
    );

    #include "setRootCase.H"
    #include "createTime.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // VARIABLE DEFINITIONS
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    IFstream inFile(args[1]);
    bool nodesNotRead(true);
    bool elemsNotRead(true);

    const bool debug = args.optionFound("debug");
    
    vector reps;

    // Storage for points
    pointField points;
    Map<label> texgenToFoam;

    // Storage for all cells.
    cellShapeList cellAsShapes;
    cellList cells;

    // Storage for all faces.
    faceList faces;
    labelList owner;
    labelList neighbour;

    // Storage for zones.
    List<DynamicList<label>> zoneCells(0);
    List<DynamicList<label>> zonePoints(0);

    // Storage for boundaries.
    labelList boundaryPatchIndices(FACEHEX);
    wordList boundaryPatchNames(FACEHEX);
    wordList boundaryTypes(FACEHEX);
    faceListList boundaryFaces(FACEHEX);
    labelListList boundaryComponents(FACEHEX);

    // Name of zones.
    Map<word> elementSets;
    Map<word> pointSets;


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // FILE PARSING
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    string line;
    inFile.getLine(line);

    while (inFile.good())
    {   
        if (line(8) == "*Heading")
        {
            readFileHeading(inFile, line);
        }
        else if (line(5) == "*Node" && nodesNotRead)
        {   
            readPoints(inFile, line, points, texgenToFoam);
            nodesNotRead = false;
        }
        else if (line(8) == "*Element" && elemsNotRead)
        {
            readCells
            (
                inFile,
                line,
                points,
                texgenToFoam,
                cellAsShapes,
                cells
            );

            setFaces
            (
                cellAsShapes, 
                points, 
                faces, 
                owner, 
                neighbour,
                cells
            );

            elemsNotRead = false;
        }
        else if (line(6) == "*ElSet")
        {
            readElSet(inFile, line, elementSets, zoneCells);
        }
        else if (line(5) == "*NSet")
        {
            readNSet(inFile, line, pointSets, texgenToFoam, zonePoints);
        }
        else
        {
            inFile.getLine(line);
        }
    }


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // CREATING BOUNDARY ELEMENTS
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    setBoundaryProperties
    (
        boundaryPatchIndices,
        boundaryPatchNames,
        boundaryTypes,
        boundaryComponents, 
        zonePoints
    );

    if (debug)
    {
        Info<< "Points" << nl << points << "\n\n\n\n" << endl;
        Info<< "Cells" << nl << cells << "\n\n\n\n" << endl;
        Info<< "Faces" << nl << faces << "\n\n\n\n" << endl;
        Info<< "Owner" << nl << owner << "\n\n\n\n" << endl;
        Info<< "Neighbour" << nl << neighbour << "\n\n\n\n" << endl;
        Info<< "CellAsShapes" << nl << cellAsShapes << "\n\n\n\n" << endl;
        Info<< "ZoneCells" << nl << zoneCells << "\n\n\n\n" << endl;
        Info<< "ZonePoints" << nl << zonePoints << "\n\n\n\n" << endl;
        Info<< "BoundComps" << nl << boundaryComponents << "\n\n\n\n" << endl;
    }

    setBoundaryElems
    (
        cellAsShapes,
        points,
        boundaryComponents,
        cells,
        faces,
        owner,
        boundaryFaces
    );

    if (debug)
    {
        Info<< "BoundFaces" << nl << boundaryFaces << "\n\n\n\n" << endl;
        Info<< "NewFaces" << nl << faces << "\n\n\n\n" << endl;
        Info<< "NewOwner" << nl << owner << "\n\n\n\n" << endl;
        Info<< "NewCells" << nl << cells << "\n\n\n\n" << endl;
    }


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // MESH CREATION
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    polyMesh mesh = createMesh
    (
        points, 
        faces, 
        owner, 
        neighbour, 
        boundaryPatchIndices, 
        boundaryPatchNames, 
        boundaryTypes, 
        boundaryFaces, 
        zoneCells,
        elementSets, 
        runTime
    );

    mesh.write();

    Info<< "Finished creating the mesh." << endl;


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // MESH REPETITION (In case -repeat option is used)
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    if (args.optionReadIfPresent("repeat", reps))
    {
        Info<< "\nRepeating mesh of the cell." << nl <<
        "X direction: "<< reps[0] << nl <<
        "Y direction: "<< reps[1] << nl <<
        "Z direction: "<< reps[2] << nl << endl;

        const point boundSize = mesh.bounds().max() - mesh.bounds().min();

        scalar moveX = boundSize.x();
        scalar moveY = boundSize.y();
        scalar moveZ = boundSize.z();

        List<DynamicList<label>> newZoneCells; 
        labelList newBoundaryPatchIndices;  
        wordList newBoundaryPatchNames;  
        faceListList newBoundaryFaces;

        mergePolyMesh newMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                runTime.constant(),
                runTime
            )
        );

        newMesh.setAxisDir(dirX);

        if (debug)
        {
            Info<< "Points" << nl << newMesh.points() << "\n\n\n\n" << endl;
            Info<< "Cells" << nl << newMesh.cells() << "\n\n\n\n" << endl;
            Info<< "Faces" << nl << newMesh.faces() << "\n\n\n\n" << endl;
            Info<< "Owner" << nl << newMesh.faceOwner() << "\n\n\n\n" << endl;
            Info<< "Neighbour" << nl << newMesh.faceNeighbour() << "\n\n\n\n" << endl;
            Info<< "ZoneCells" << nl << newMesh.cellZones() << "\n\n\n\n" << endl;
            Info<< "Bounds" << nl << newMesh.boundaryMesh() << "\n\n\n\n" << endl;
        }

        for(label i=1; i < reps[0]; i++)
        {
            Info<< "Joining cells " << i << " and " << (i+1) <<
                " in the X direction." << endl;

            pointField addPoints = mesh.points();
            faceList addFaces = mesh.faces();
            labelList addOwner = mesh.faceOwner();
            labelList addNeighbour = mesh.faceNeighbour();

            forAll(addPoints, pointi)
            {
                addPoints[pointi].x() += moveX;
            }

            // Map<word> newBoundaryPatchNames = boundaryPatchNames;
            // newBoundaryPatchNames.erase(1);
            // newBoundaryPatchNames.insert(1, "MergingLeft");

            polyMesh addMesh = createMesh
            (
                addPoints, 
                addFaces,
                addOwner,
                addNeighbour,
                boundaryPatchIndices, 
                boundaryPatchNames, 
                boundaryTypes, 
                boundaryFaces, 
                zoneCells,
                elementSets, 
                runTime
            );

            if (debug)
            {
                Info<< "Points" << nl << addMesh.points() << "\n\n\n\n" << endl;
                Info<< "Cells" << nl << addMesh.cells() << "\n\n\n\n" << endl;
                Info<< "Faces" << nl << addMesh.faces() << "\n\n\n\n" << endl;
                Info<< "Owner" << nl << addMesh.faceOwner() << "\n\n\n\n" << endl;
                Info<< "Neighbour" << nl << addMesh.faceNeighbour() << "\n\n\n\n" << endl;
                Info<< "ZoneCells" << nl << addMesh.cellZones() << "\n\n\n\n" << endl;
                Info<< "Bounds" << nl << addMesh.boundaryMesh() << "\n\n\n\n" << endl;
            }

            newMesh.addMesh(addMesh);
            newMesh.merge();

            moveX += boundSize.x();
        }

        newMesh.setAxisDir(dirY);
        newMesh.getProps
        (
            newBoundaryPatchIndices, 
            newBoundaryPatchNames, 
            newBoundaryFaces, 
            newZoneCells
        );

        for(label j=1; j < reps[1]; j++)
        {
            Info<< "Joining cells " << j << " and " << (j+1) <<
                " in the Y direction." << endl;

            pointField addPoints = newMesh.points();
            faceList addFaces = newMesh.faces();
            labelList addOwner = newMesh.faceOwner();
            labelList addNeighbour = newMesh.faceNeighbour();
            
            // Map<word> newBoundaryPatchNames;

            // forAll(newMesh.boundaryMesh(), bd)
            // {
            //     const polyPatch& curPatch = newMesh.boundaryMesh()[bd];
            //     newBoundaryPatchNames.insert
            //         (
            //             curPatch.index(),
            //             (curPatch.name() == "Front" ? "MergingFront" : "Front")
            //         );

            // }
            // newBoundaryPatchNames.erase(3);
            // newBoundaryPatchNames.insert(3, "MergingFront");
            
            forAll(addPoints, pointi)
            {
                addPoints[pointi].y() += moveY;
            }

            polyMesh addMesh = createMesh
            (
                addPoints, 
                addFaces, 
                addOwner, 
                addNeighbour,
                newBoundaryPatchIndices, 
                newBoundaryPatchNames, 
                boundaryTypes, 
                newBoundaryFaces, 
                newZoneCells,
                elementSets, 
                runTime
            );

            newMesh.addMesh(addMesh);
            newMesh.merge();

            moveY += boundSize.y();
        }

        newMesh.setAxisDir(dirZ);
        newMesh.getProps
        (
            newBoundaryPatchIndices, 
            newBoundaryPatchNames, 
            newBoundaryFaces, 
            newZoneCells
        );

        for(label k=1; k < reps[2]; k++)
        {
            Info<< "Joining cells " << k << " and " << (k+1) <<
                " in the Z direction." << endl;

            pointField addPoints = newMesh.points();
            faceList addFaces = newMesh.faces();
            labelList addOwner = newMesh.faceOwner();
            labelList addNeighbour = newMesh.faceNeighbour();
            
            // Map<word> newBoundaryPatchNames;

            // forAll(newMesh.boundaryMesh(), bd)
            // {
            //     const polyPatch& curPatch = newMesh.boundaryMesh()[bd];
            //     newBoundaryPatchNames.insert
            //         (
            //             curPatch.index(),
            //             (curPatch.name() == "Bottom" ? "MergingBottom" : "Bottom")
            //         );

            // }
            
            forAll(addPoints, pointi)
            {
                addPoints[pointi].z() += moveZ;
            }

            polyMesh addMesh = createMesh
            (
                addPoints, 
                addFaces, 
                addOwner, 
                addNeighbour,
                newBoundaryPatchIndices, 
                newBoundaryPatchNames, 
                boundaryTypes, 
                newBoundaryFaces, 
                newZoneCells,
                elementSets, 
                runTime
            );

            newMesh.addMesh(addMesh);
            newMesh.merge();

            moveZ += boundSize.z();
        }


        newMesh.setInstance(mesh.pointsInstance());
        newMesh.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
