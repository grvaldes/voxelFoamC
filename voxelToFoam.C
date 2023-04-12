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
static label NODQUAD  = 4;
static label NODHEX   = 8;


// Comma delimited row parsing.
List<string> headerParse(IStringStream& lineStr)
{
    List<string> row;

    variable tag(lineStr);
    label pos;
    label subpos;
    label rowI(0);

    while ((pos = tag.find(",")) != -1)
    {
        row.append(tag.substr(0, pos));

        if ((subpos = row[rowI].find("=")) != -1)
        {
            row[rowI] = row[rowI].substr(subpos, row[rowI].length() - subpos);
        }

        tag.erase(0, pos+1);
        rowI++;
    }

    return row;
}

// Comma delimited row parsing.
List<scalar> scalarParse(IStringStream& lineStr, label ListSize)
{
    List<scalar> row(ListSize);

    variable tag(lineStr);
    label pos;
    label rowI(0);

    while ((pos = tag.find(",")) != -1)
    {
        if (ListSize == SIZE_NOT_DEFINED)
        {
            row.append(stod(tag.substr(0, pos)));
        }
        else
        {
            row[rowI] = stod(tag.substr(0, pos));
        }

        tag.erase(0, pos+1);
        rowI++;
    }

    return row;
}

// Comma delimited row parsing.
List<label> labelParse(IStringStream& lineStr, label ListSize)
{
    List<label> row(ListSize);

    variable tag(lineStr);
    label pos;
    label rowI(0);

    while ((pos = tag.find(",")) != -1)
    {
        if (ListSize == SIZE_NOT_DEFINED)
        {
            row.append(stoi(tag.substr(0, pos)));
        }
        else
        {
            row[rowI] = stoi(tag.substr(0, pos));
        }

        tag.erase(0, pos+1);
        rowI++;
    }

    return row;
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


// Find face in pp which uses all vertices in meshF (in mesh point labels)
label findFace(const primitivePatch& pp, const labelList& meshF)
{
    const Map<label>& meshPointMap = pp.meshPointMap();

    // meshF[0] in pp labels.
    if (!meshPointMap.found(meshF[0]))
    {
        Warning<< "Not using gmsh face " << meshF
            << " since zero vertex is not on boundary of polyMesh" << endl;
        return -1;
    }

    // Find faces using first point
    const labelList& pFaces = pp.pointFaces()[meshPointMap[meshF[0]]];

    // Go through all these faces and check if there is one which uses all of
    // meshF vertices (in any order ;-)
    forAll(pFaces, i)
    {
        label facei = pFaces[i];

        const face& f = pp[facei];

        // Count uses of vertices of meshF for f
        label nMatched = 0;

        forAll(f, fp)
        {
            if (findIndex(meshF, f[fp]) != -1)
            {
                nMatched++;
            }
        }

        if (nMatched == meshF.size())
        {
            return facei;
        }
    }

    return -1;
}


// Same but find internal face. Expensive addressing.
label findInternalFace(const primitiveMesh& mesh, const labelList& meshF)
{
    const labelList& pFaces = mesh.pointFaces()[meshF[0]];

    forAll(pFaces, i)
    {
        label facei = pFaces[i];

        const face& f = mesh.faces()[facei];

        // Count uses of vertices of meshF for f
        label nMatched = 0;

        forAll(f, fp)
        {
            if (findIndex(meshF, f[fp]) != -1)
            {
                nMatched++;
            }
        }

        if (nMatched == meshF.size())
        {
            return facei;
        }
    }
    return -1;
}


// Determine whether cell is inside-out by checking for any wrong-oriented
// face.
bool correctOrientation(const pointField& points, const cellShape& shape)
{
    // Get centre of shape.
    point cc(shape.centre(points));

    // Get outwards pointing faces.
    faceList faces(shape.faces());

    forAll(faces, i)
    {
        const face& f = faces[i];

        const vector a(f.area(points));

        // Check if vector from any point on face to cc points outwards
        if (((points[f[0]] - cc) & a) < 0)
        {
            // Incorrectly oriented
            return false;
        }
    }

    return true;
}


// Reads mesh format
void readFileHeading(IFstream& inFile)
{
    string line;
    inFile.getLine(line);

    Info<< line << endl;

    inFile.getLine(line);
}


// Reads points and map
void readPoints(IFstream& inFile, pointField& points, Map<label>& texgenToFoam)
{
    Info<< "Starting to read points at line " << inFile.lineNumber() << endl;

    label pointi = 0;

    while (inFile.good())
    { 
        string line;
        inFile.getLine(line);
        IStringStream lineStr(line);
        variable tag(lineStr);

        if (tag.substr(0,1) == "*")
        {
            Info<< "Finished reading nodes. Vertices read:"
                << texgenToFoam.size()
                << endl;
            break;
        }

        List<scalar> row = scalarParse(lineStr, 4);

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
    const bool keepOrientation,
    const pointField& points,
    const Map<label>& texgenToFoam,
    IFstream& inFile,
    cellShapeList& cells
)
{
    label& lineNumber = inFile.lineNumber();

    Info<< "Starting to read cells at line " << lineNumber << endl;

    const cellModel& hex = *(cellModeller::lookup("hex"));

    face triPoints(3);
    face quadPoints(4);
    labelList tetPoints(4);
    labelList pyrPoints(5);
    labelList prismPoints(6);
    labelList hexPoints(8);

    lineNumber--;

    string line;
    inFile.getLine(line);
    IStringStream lineStr(line);

    List<string> header = headerParse(lineStr);

    string elmType = header[1];
    label celli = 0;

    while (inFile.good())
    {
        string line;
        inFile.getLine(line);
        IStringStream lineStr(line);
        variable tag(lineStr);

        if (tag.substr(0,1) == "*")
        {
            Info<< "Finished reading cells. Elements read:"
                << texgenToFoam.size()
                << endl;
            break;
        }

        List<label> row = labelParse(lineStr, NODHEX+1);

        forAll (hexPoints, i)
        {
            hexPoints[i] = row[i+1];
        }

        renumber(texgenToFoam, hexPoints);
        cells.append( cellShape(hex, hexPoints) );

        const cellShape& cell = cells[celli];

        if (!keepOrientation && !correctOrientation(points, cell))
        {
            Info<< "Inverting hex " << celli << endl;
            // Reorder hex.
            hexPoints[0] = cell[4];
            hexPoints[1] = cell[5];
            hexPoints[2] = cell[6];
            hexPoints[3] = cell[7];
            hexPoints[4] = cell[0];
            hexPoints[5] = cell[1];
            hexPoints[6] = cell[2];
            hexPoints[7] = cell[3];

            cells[celli] = cellShape(hex, hexPoints);
        }

        celli++;
    }

    Info<< "Cells:" << endl
    << "    total:" << cells.size()
    << endl;

    if (cells.size() == 0)
    {
        FatalIOErrorInFunction(inFile)
            << "No cells read from file " << inFile.name() << nl
            << exit(FatalIOError);
    }
}


// Reads physical names
void readElSet
(
    IFstream& inFile, 
    Map<word>& elementSets,
    List<DynamicList<label>>& zoneCells
)
{
    Info<< "Starting to read element set at line " << inFile.lineNumber()
        << endl;

    label zoneI = zoneCells.size();
    zoneCells.setSize(zoneI+1);

    label& lineNumber = inFile.lineNumber();
    lineNumber--;

    string line;
    inFile.getLine(line);
    IStringStream lineStr(line);

    List<string> header = headerParse(lineStr);
    string regionName = header[1];

    Info<< "Mapping Element Set " << regionName
        << " to Foam cellZone " << zoneI 
        << endl;

    elementSets.insert(zoneI, string::validate<word>(regionName));
    
    while (inFile.good())
    {
        string line;
        inFile.getLine(line);
        IStringStream lineStr(line);
        variable tag(lineStr);

        if (tag.substr(0,1) == "*")
        {
            Info<< "Finished reading ElSet. Elements read:"
                << zoneCells[zoneI].size()
                << endl;
            break;
        }

        List<label> row = labelParse(lineStr, SIZE_NOT_DEFINED);

        for (int i = 1; i < row.size(); i++)
        {
            zoneCells[zoneI].append(row[i]-1);
        }

        inFile.getLine(line);
        IStringStream tagStr(line);
    }
}


// Reads physical names
void readNSet
(
    IFstream& inFile, 
    Map<word>& pointSets,
    const Map<label>& texgenToFoam,
    List<DynamicList<label>>& zonePoints
)
{
    Info<< "Starting to read point set at line " << inFile.lineNumber()
        << endl;

    label zoneI = zonePoints.size();
    zonePoints.setSize(zoneI+1);

    label& lineNumber = inFile.lineNumber();
    lineNumber--;

    string line;
    inFile.getLine(line);
    IStringStream lineStr(line);

    List<string> header = headerParse(lineStr);
    string regionName = header[1];

    Info<< "Mapping Point Set " << regionName
        << " to Foam cellZone " << zoneI 
        << endl;

    pointSets.insert(zoneI, string::validate<word>(regionName));
    
    while (inFile.good())
    {
        string line;
        inFile.getLine(line);
        IStringStream lineStr(line);
        variable tag(lineStr);

        if (tag.substr(0,1) == "*")
        {
            Info<< "Finished reading NSet. Nodes read:"
                << zonePoints[zoneI]
                << endl;
            break;
        }

        List<label> row = labelParse(lineStr, SIZE_NOT_DEFINED);
        renumber(texgenToFoam, row);

        for (int i = 1; i < row.size(); i++)
        {
            zonePoints[zoneI].append(row[i]);
        }
    }
}



int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append(".msh file");
    argList::addBoolOption
    (
        "keepOrientation",
        "retain raw orientation for prisms/hexs"
    );

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

    const bool keepOrientation = args.optionFound("keepOrientation");
    IFstream inFile(args[1]);
    bool nodesNotRead(true);
    bool elemsNotRead(true);

    // Storage for points
    pointField points;
    Map<label> texgenToFoam;

    // Storage for all cells.
    cellShapeList cells;

    // Storage for zones.
    List<DynamicList<face>> patchFaces(0);
    List<DynamicList<label>> zoneCells(0);
    List<DynamicList<label>> zonePoints(0);

    // Name of zones.
    Map<word> elementSets;
    Map<word> pointSets;

    string line;
    inFile.getLine(line);
    IStringStream lineStr(line);

    while (inFile.good())
    {
        variable tag(lineStr);

        if (tag.substr(0, 8) == "*Heading")
        {
            readFileHeading(inFile);
        }
        else if (tag.substr(0, 5) == "*Node")
        {
            if (nodesNotRead)
            {
                readPoints(inFile, points, texgenToFoam);
                nodesNotRead = false;
            }
        }
        else if (tag.substr(0, 8) == "*Element")
        {
            if (elemsNotRead)
            {
                readCells
                (
                    keepOrientation,
                    points,
                    texgenToFoam,
                    inFile,
                    cells
                );
                elemsNotRead = false;
            }
        }
        else if (tag.substr(0, 6) == "*Elset")
        {
            readElSet(inFile, elementSets, zoneCells);
        }
        else if (tag.substr(0, 5) == "*NSet")
        {
            readNSet(inFile, pointSets, texgenToFoam, zonePoints);
        }
        else
        {
            string line;
            inFile.getLine(line);
            IStringStream lineStr(line);
        }
    }


    label nValidCellZones = 0;

    forAll(zoneCells, zoneI)
    {
        if (zoneCells[zoneI].size())
        {
            nValidCellZones++;
        }
    }


    // // Problem is that the orientation of the patchFaces does not have to
    // // be consistent with the outwards orientation of the mesh faces. So
    // // we have to construct the mesh in two stages:
    // // 1. define mesh with all boundary faces in one patch
    // // 2. use the read patchFaces to find the corresponding boundary face
    // //    and repatch it.


    // // Create correct number of patches
    // // (but without any faces in it)
    // faceListList boundaryFaces(patchFaces.size());

    // wordList boundaryPatchNames(boundaryFaces.size());

    // forAll(boundaryPatchNames, patchi)
    // {
    //     label physReg = patchToNset[patchi];

    //     Map<word>::const_iterator iter = elementSets.find(physReg);

    //     if (iter != elementSets.end())
    //     {
    //         boundaryPatchNames[patchi] = iter();
    //     }
    //     else
    //     {
    //         boundaryPatchNames[patchi] = word("patch") + name(patchi);
    //     }
    //     Info<< "Patch " << patchi << " gets name "
    //         << boundaryPatchNames[patchi] << endl;
    // }
    // Info<< endl;

    // wordList boundaryPatchTypes(boundaryFaces.size(), polyPatch::typeName);
    // word defaultFacesName = "defaultFaces";
    // word defaultFacesType = polyPatch::typeName;
    // wordList boundaryPatchPhysicalTypes
    // (
    //     boundaryFaces.size(),
    //     polyPatch::typeName
    // );

    // polyMesh mesh
    // (
    //     IOobject
    //     (
    //         regionName,
    //         runTime.constant(),
    //         runTime
    //     ),
    //     move(points),
    //     cells,
    //     boundaryFaces,
    //     boundaryPatchNames,
    //     boundaryPatchTypes,
    //     defaultFacesName,
    //     defaultFacesType,
    //     boundaryPatchPhysicalTypes
    // );

    // repatchPolyTopoChanger repatcher(mesh);

    // // Now use the patchFaces to patch up the outside faces of the mesh.

    // // Get the patch for all the outside faces (= default patch added as last)
    // const polyPatch& pp = mesh.boundaryMesh().last();

    // // Storage for faceZones.
    // List<DynamicList<label>> zoneFaces(patchFaces.size());


    // // Go through all the patchFaces and find corresponding face in pp.
    // forAll(patchFaces, patchi)
    // {
    //     const DynamicList<face>& pFaces = patchFaces[patchi];

    //     Info<< "Finding faces of patch " << patchi << endl;

    //     forAll(pFaces, i)
    //     {
    //         const face& f = pFaces[i];

    //         // Find face in pp using all vertices of f.
    //         label patchFacei = findFace(pp, f);

    //         if (patchFacei != -1)
    //         {
    //             label meshFacei = pp.start() + patchFacei;

    //             repatcher.changePatchID(meshFacei, patchi);
    //         }
    //         else
    //         {
    //             // Maybe internal face? If so add to faceZone with same index
    //             // - might be useful.
    //             label meshFacei = findInternalFace(mesh, f);

    //             if (meshFacei != -1)
    //             {
    //                 zoneFaces[patchi].append(meshFacei);
    //             }
    //             else
    //             {
    //                 WarningInFunction
    //                     << "Could not match gmsh face " << f
    //                     << " to any of the interior or exterior faces"
    //                     << " that share the same 0th point" << endl;
    //             }
    //         }
    //     }
    // }
    // Info<< nl;

    // // Face zones
    // label nValidFaceZones = 0;

    // Info<< "FaceZones:" << nl
    //     << "Zone\tSize" << endl;

    // forAll(zoneFaces, zoneI)
    // {
    //     zoneFaces[zoneI].shrink();

    //     const labelList& zFaces = zoneFaces[zoneI];

    //     if (zFaces.size())
    //     {
    //         nValidFaceZones++;

    //         Info<< "    " << zoneI << '\t' << zFaces.size() << endl;
    //     }
    // }
    // Info<< endl;


    // // Get polyMesh to write to constant

    // runTime.setTime(instant(runTime.constant()), 0);

    // repatcher.repatch();

    // List<cellZone*> cz;
    // List<faceZone*> fz;

    // // Construct and add the zones. Note that cell ordering does not change
    // // because of repatch() and neither does internal faces so we can
    // // use the zoneCells/zoneFaces as is.

    // if (nValidCellZones > 0)
    // {
    //     cz.setSize(nValidCellZones);

    //     nValidCellZones = 0;

    //     forAll(zoneCells, zoneI)
    //     {
    //         if (zoneCells[zoneI].size())
    //         {
    //             label physReg = zoneToElset[zoneI];

    //             Map<word>::const_iterator iter = elementSets.find(physReg);

    //             word zoneName = "cellZone_" + name(zoneI);
    //             if (iter != elementSets.end())
    //             {
    //                 zoneName = iter();
    //             }

    //             Info<< "Writing zone " << zoneI << " to cellZone "
    //                 << zoneName << " and cellSet"
    //                 << endl;

    //             cellSet cset(mesh, zoneName, zoneCells[zoneI]);
    //             cset.write();

    //             cz[nValidCellZones] = new cellZone
    //             (
    //                 zoneName,
    //                 zoneCells[zoneI],
    //                 nValidCellZones,
    //                 mesh.cellZones()
    //             );
    //             nValidCellZones++;
    //         }
    //     }
    // }

    // if (nValidFaceZones > 0)
    // {
    //     fz.setSize(nValidFaceZones);

    //     nValidFaceZones = 0;

    //     forAll(zoneFaces, zoneI)
    //     {
    //         if (zoneFaces[zoneI].size())
    //         {
    //             label physReg = patchToNset[zoneI];

    //             Map<word>::const_iterator iter = elementSets.find(physReg);

    //             word zoneName = "faceZone_" + name(zoneI);
    //             if (iter != elementSets.end())
    //             {
    //                 zoneName = iter();
    //             }

    //             Info<< "Writing zone " << zoneI << " to faceZone "
    //                 << zoneName << " and faceSet"
    //                 << endl;

    //             faceSet fset(mesh, zoneName, zoneFaces[zoneI]);
    //             fset.write();

    //             fz[nValidFaceZones] = new faceZone
    //             (
    //                 zoneName,
    //                 zoneFaces[zoneI],
    //                 boolList(zoneFaces[zoneI].size(), true),
    //                 nValidFaceZones,
    //                 mesh.faceZones()
    //             );
    //             nValidFaceZones++;
    //         }
    //     }
    // }

    // if (cz.size() || fz.size())
    // {
    //     mesh.addZones(List<pointZone*>(0), fz, cz);
    // }

    // // Remove empty defaultFaces
    // label defaultPatchID = mesh.boundaryMesh().findPatchID(defaultFacesName);
    // if (mesh.boundaryMesh()[defaultPatchID].size() == 0)
    // {
    //     List<polyPatch*> newPatchPtrList((mesh.boundaryMesh().size() - 1));
    //     label newPatchi = 0;
    //     forAll(mesh.boundaryMesh(), patchi)
    //     {
    //         if (patchi != defaultPatchID)
    //         {
    //             const polyPatch& patch = mesh.boundaryMesh()[patchi];

    //             newPatchPtrList[newPatchi] = patch.clone
    //             (
    //                 mesh.boundaryMesh(),
    //                 newPatchi,
    //                 patch.size(),
    //                 patch.start()
    //             ).ptr();

    //             newPatchi++;
    //         }
    //     }
    //     repatcher.changePatches(newPatchPtrList);
    // }

    // mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
