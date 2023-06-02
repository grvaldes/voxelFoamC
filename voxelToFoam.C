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
#include "polyMesh.H"
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

// Assign OF numbers to TG nodes.
void renumber
(
    const Map<label>& texgenToFoam,
    labelList& labels
)
{
    forAll(labels, labelI)
    {
        labels[labelI] = texgenToFoam[labels[labelI]];
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


// Copy from Foam::labelList Foam::polyMesh::facePatchFaceCells
labelList facePatchFaceCells
(
    const faceList& patchFaces,
    const labelListList& pointCells,
    const faceListList& cellsFaceShapes,
    const label patchID
)
{
    bool found;

    labelList FaceCells(patchFaces.size());

    forAll(patchFaces, fI)
    {
        found = false;

        const face& curFace = patchFaces[fI];
        const labelList& facePoints = patchFaces[fI];

        forAll(facePoints, pointi)
        {
            const labelList& facePointCells = pointCells[facePoints[pointi]];

            forAll(facePointCells, celli)
            {
                faceList cellFaces = cellsFaceShapes[facePointCells[celli]];

                forAll(cellFaces, cellFace)
                {
                    if (face::sameVertices(cellFaces[cellFace], curFace))
                    {
                        // Found the cell corresponding to this face
                        FaceCells[fI] = facePointCells[celli];

                        found = true;
                    }
                    if (found) break;
                }
                if (found) break;
            }
            if (found) break;
        }

        if (!found)
        {
            FatalErrorInFunction
                << "face " << fI << " in patch " << patchID
                << " does not have neighbour cell"
                << " face: " << patchFaces[fI]
                << abort(FatalError);
        }
    }

    return FaceCells;
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
            Info<< "Finished reading nodes. Vertices read:"
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

        renumber(texgenToFoam, hexPoints);
        cellAsShapes.append( cellShape(hex, hexPoints) );
        cells.append( cell( FACEHEX ) );

        celli++;
    }

    if (cells.size() == 0)
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
            << " to Foam cellZone " << zoneI 
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
        // regionName(6) != "Master" &&
        regionName(3) != "All"
    )
    {
        zonePoints.setSize(zoneI+1);
        
        Info<< "Mapping Point Set " << regionName
            << " to Foam cellZone " << zoneI 
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
            renumber(texgenToFoam, row);

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
    Map<word>& boundaryPatchNames,
    Map<word>& boundaryPatchType,
    Map<word>& boundaryPhysicalType,
    labelListList& boundaryComponents,
    List<DynamicList<label>>& zonePoints
)
{
    boundaryPatchNames.set(0,"Front");
    boundaryPatchNames.set(1,"Back");
    boundaryPatchNames.set(2,"Left");
    boundaryPatchNames.set(3,"Right");
    boundaryPatchNames.set(4,"Bottom");
    boundaryPatchNames.set(5,"Top");

    forAll(boundaryPatchNames, patchi)
    {
        boundaryPatchType.set(patchi, "patch");
        boundaryPhysicalType.set(patchi, "patch");
    }

    forAll(boundaryComponents, bc)
    {
        switch (bc)
        {
            case 0:
                boundaryComponents[bc].append(zonePoints[0]);
                boundaryComponents[bc].append(zonePoints[6]);
                boundaryComponents[bc].append(zonePoints[7]);
                boundaryComponents[bc].append(zonePoints[8]);
                boundaryComponents[bc].append(zonePoints[9]);
                boundaryComponents[bc].append(zonePoints[18]);
                boundaryComponents[bc].append(zonePoints[21]);
                boundaryComponents[bc].append(zonePoints[22]);
                boundaryComponents[bc].append(zonePoints[25]);
                break;
            case 1:
                boundaryComponents[bc].append(zonePoints[1]);
                boundaryComponents[bc].append(zonePoints[10]);
                boundaryComponents[bc].append(zonePoints[11]);
                boundaryComponents[bc].append(zonePoints[12]);
                boundaryComponents[bc].append(zonePoints[13]);
                boundaryComponents[bc].append(zonePoints[19]);
                boundaryComponents[bc].append(zonePoints[20]);
                boundaryComponents[bc].append(zonePoints[23]);
                boundaryComponents[bc].append(zonePoints[24]);
                break;
            case 2:
                boundaryComponents[bc].append(zonePoints[2]);
                boundaryComponents[bc].append(zonePoints[9]);
                boundaryComponents[bc].append(zonePoints[13]);
                boundaryComponents[bc].append(zonePoints[16]);
                boundaryComponents[bc].append(zonePoints[17]);
                boundaryComponents[bc].append(zonePoints[20]);
                boundaryComponents[bc].append(zonePoints[21]);
                boundaryComponents[bc].append(zonePoints[24]);
                boundaryComponents[bc].append(zonePoints[25]);
                break;
            case 3:
                boundaryComponents[bc].append(zonePoints[3]);
                boundaryComponents[bc].append(zonePoints[7]);
                boundaryComponents[bc].append(zonePoints[11]);
                boundaryComponents[bc].append(zonePoints[14]);
                boundaryComponents[bc].append(zonePoints[15]);
                boundaryComponents[bc].append(zonePoints[18]);
                boundaryComponents[bc].append(zonePoints[19]);
                boundaryComponents[bc].append(zonePoints[22]);
                boundaryComponents[bc].append(zonePoints[23]);
                break;
            case 4:
                boundaryComponents[bc].append(zonePoints[4]);
                boundaryComponents[bc].append(zonePoints[6]);
                boundaryComponents[bc].append(zonePoints[10]);
                boundaryComponents[bc].append(zonePoints[14]);
                boundaryComponents[bc].append(zonePoints[17]);
                boundaryComponents[bc].append(zonePoints[18]);
                boundaryComponents[bc].append(zonePoints[19]);
                boundaryComponents[bc].append(zonePoints[20]);
                boundaryComponents[bc].append(zonePoints[21]);
                break;
            case 5:
                boundaryComponents[bc].append(zonePoints[5]);
                boundaryComponents[bc].append(zonePoints[8]);
                boundaryComponents[bc].append(zonePoints[12]);
                boundaryComponents[bc].append(zonePoints[15]);
                boundaryComponents[bc].append(zonePoints[16]);
                boundaryComponents[bc].append(zonePoints[22]);
                boundaryComponents[bc].append(zonePoints[23]);
                boundaryComponents[bc].append(zonePoints[24]);
                boundaryComponents[bc].append(zonePoints[25]);
                break;
        }
    }
}


// Setting boundary elements
void setBoundaryElems
(
    const cellShapeList& cellsAsShapes,
    const pointField& points,
    const Map<word>& boundaryPatchNames,
    const labelListList& boundaryComponents,
    const List<DynamicList<label>>& zonePoints,
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




int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append(".inp file");

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    Foam::word regionName;

    if (args.optionReadIfPresent("region", regionName))
    {
        Foam::Info
            << "Creating polyMesh for region " << regionName << endl;
    }
    else
    {
        regionName = Foam::polyMesh::defaultRegion;
    }

    IFstream inFile(args[1]);
    bool nodesNotRead(true);
    bool elemsNotRead(true);

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
    List<DynamicList<face>> patchFaces(0);
    List<DynamicList<label>> zoneCells(0);
    List<DynamicList<label>> zonePoints(0);

    // Storage for boundaries.
    Map<word> boundaryPatchNames(FACEHEX);
    Map<word> boundaryTypes(FACEHEX);
    Map<word> boundaryPhysicalTypes(FACEHEX);
    faceListList boundaryFaces(FACEHEX);
    labelListList boundaryComponents(FACEHEX);

    // Name of zones.
    Map<word> elementSets;
    Map<word> pointSets;
    Map<word> faceSets;

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

    setBoundaryProperties
    (
        boundaryPatchNames,
        boundaryTypes,
        boundaryPhysicalTypes,
        boundaryComponents, 
        zonePoints
    );


    Info<< "Points" << nl << points << "\n\n\n\n" << endl;
    Info<< "Cells" << nl << cells << "\n\n\n\n" << endl;
    Info<< "Faces" << nl << faces << "\n\n\n\n" << endl;
    Info<< "Owner" << nl << owner << "\n\n\n\n" << endl;
    Info<< "Neighbour" << nl << neighbour << "\n\n\n\n" << endl;
    Info<< "ElSets" << nl << elementSets << "\n\n\n\n" << endl;
    Info<< "NSets" << nl << pointSets << "\n\n\n\n" << endl;
    Info<< "CellAsShapes" << nl << cellAsShapes << "\n\n\n\n" << endl;
    Info<< "TexToFoam" << nl << texgenToFoam << "\n\n\n\n" << endl;
    Info<< "ZoneCells" << nl << zoneCells << "\n\n\n\n" << endl;
    Info<< "ZonePoints" << nl << zonePoints << "\n\n\n\n" << endl;
    Info<< "BoundNames" << nl << boundaryPatchNames << "\n\n\n\n" << endl;
    Info<< "BoundTypes" << nl << boundaryTypes << "\n\n\n\n" << endl;
    Info<< "BoundPhysTypes" << nl << boundaryPhysicalTypes << "\n\n\n\n" << endl;
    Info<< "BoundComps" << nl << boundaryComponents << "\n\n\n\n" << endl;


    setBoundaryElems
    (
        cellAsShapes,
        points,
        boundaryPatchNames,
        boundaryComponents,
        zonePoints,
        cells,
        faces,
        owner
    );

    // reorderFaces
    // (
    //     points,
    //     cells,
    //     faces,
    //     owner,
    //     boundaryFaces
    // )

    Info<< "BoundFaces" << nl << boundaryFaces << "\n\n\n\n" << endl;
    Info<< "NewFaces" << nl << faces << "\n\n\n\n" << endl;
    Info<< "NewOwner" << nl << owner << "\n\n\n\n" << endl;
    Info<< "NewCells" << nl << cells << "\n\n\n\n" << endl;


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

    

    mesh.write();

    Info<< "Mesh" << nl << mesh.bounds() << "\n\n\n\n" << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
