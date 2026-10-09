// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#include "VtkMeshExporter.hpp"

#include "IMeshDAO.hpp"
#include "IOUtilities.hpp"
#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
namespace External
{
namespace MeshUtilities
{
// ***************************************************************************
void ExportMeshToVTU(const Gedim::IMeshDAO &mesh, const string &exportFolder, const string &fileName, const bool &separateFile)
{
    string cell0DsFolder = exportFolder;
    string cell1DsFolder = exportFolder;
    string cell2DsFolder = exportFolder;
    string cell3DsFolder = exportFolder;

    if (separateFile)
    {
        Gedim::Output::CreateFolder(exportFolder + "/Cell0Ds");
        Gedim::Output::CreateFolder(exportFolder + "/Cell1Ds");
        Gedim::Output::CreateFolder(exportFolder + "/Cell2Ds");
        Gedim::Output::CreateFolder(exportFolder + "/Cell3Ds");

        cell0DsFolder = exportFolder + "/Cell0Ds";
        cell1DsFolder = exportFolder + "/Cell1Ds";
        cell2DsFolder = exportFolder + "/Cell2Ds";
        cell3DsFolder = exportFolder + "/Cell3Ds";
    }

    // Export Cell0Ds
    if (mesh.Cell0DTotalNumber() > 0)
    {
        vector<double> id(mesh.Cell0DTotalNumber());
        vector<double> marker(mesh.Cell0DTotalNumber());
        vector<double> active(mesh.Cell0DTotalNumber());

        for (unsigned int g = 0; g < mesh.Cell0DTotalNumber(); g++)
        {
            id[g] = g;
            marker[g] = mesh.Cell0DMarker(g);
            active[g] = mesh.Cell0DIsActive(g);
        }

        vector<VTPProperty> properties(3 + mesh.Cell0DNumberDoubleProperties());
        vector<vector<double>> propertyValues(mesh.Cell0DNumberDoubleProperties());

        properties[0] = {"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()};
        properties[1] = {"Marker", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(marker.size()), marker.data()};
        properties[2] = {"Active", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(active.size()), active.data()};

        for (unsigned int p = 0; p < mesh.Cell0DNumberDoubleProperties(); p++)
        {
            propertyValues[p].resize(mesh.Cell0DTotalNumber());
            for (unsigned int g = 0; g < mesh.Cell0DTotalNumber(); g++)
            {
                propertyValues[p][g] = mesh.Cell0DDoublePropertySize(g, p) == 1 ? mesh.Cell0DDoublePropertyValue(g, p, 0) : 0.0;
            }

            properties[3 + p] = {mesh.Cell0DDoublePropertyId(p),
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(propertyValues[p].size()),
                                 propertyValues[p].data()};
        }

        Gedim::VTKUtilities vtpUtilities;
        vtpUtilities.AddPoints(mesh.Cell0DsCoordinates(), properties);
        vtpUtilities.Export(cell0DsFolder + "/" + "Cell0Ds_" + fileName + ".vtu");
    }

    // Export Cell1Ds
    if (mesh.Cell1DTotalNumber() > 0)
    {
        vector<double> id(mesh.Cell1DTotalNumber());
        vector<double> marker(mesh.Cell1DTotalNumber());
        vector<double> active(mesh.Cell1DTotalNumber());

        for (unsigned int g = 0; g < mesh.Cell1DTotalNumber(); g++)
        {
            id[g] = g;
            marker[g] = mesh.Cell1DMarker(g);
            active[g] = mesh.Cell1DIsActive(g);
        }

        vector<VTPProperty> properties(3 + mesh.Cell1DNumberDoubleProperties());
        vector<vector<double>> propertyValues(mesh.Cell1DNumberDoubleProperties());

        properties[0] = {"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()};
        properties[1] = {"Marker", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(marker.size()), marker.data()};
        properties[2] = {"Active", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(active.size()), active.data()};

        for (unsigned int p = 0; p < mesh.Cell1DNumberDoubleProperties(); p++)
        {
            propertyValues[p].resize(mesh.Cell1DTotalNumber());
            for (unsigned int g = 0; g < mesh.Cell1DTotalNumber(); g++)
            {
                propertyValues[p][g] = mesh.Cell1DDoublePropertySize(g, p) == 1 ? mesh.Cell1DDoublePropertyValue(g, p, 0) : 0.0;
            }

            properties[3 + p] = {mesh.Cell1DDoublePropertyId(p),
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(propertyValues[p].size()),
                                 propertyValues[p].data()};
        }

        Gedim::VTKUtilities vtpUtilities;
        vtpUtilities.AddSegments(mesh.Cell0DsCoordinates(), mesh.Cell1DsExtremes(), properties);
        vtpUtilities.Export(cell1DsFolder + "/" + "Cell1Ds_" + fileName + ".vtu");
    }

    // Export Cell2Ds
    if (mesh.Cell2DTotalNumber() > 0)
    {
        vector<double> id(mesh.Cell2DTotalNumber());
        vector<double> marker(mesh.Cell2DTotalNumber());
        vector<double> active(mesh.Cell2DTotalNumber());

        for (unsigned int g = 0; g < mesh.Cell2DTotalNumber(); g++)
        {
            id[g] = g;
            marker[g] = mesh.Cell2DMarker(g);
            active[g] = mesh.Cell2DIsActive(g);
        }

        vector<VTPProperty> properties(3 + mesh.Cell2DNumberDoubleProperties());
        vector<vector<double>> propertyValues(mesh.Cell2DNumberDoubleProperties());

        properties[0] = {"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()};
        properties[1] = {"Marker", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(marker.size()), marker.data()};
        properties[2] = {"Active", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(active.size()), active.data()};

        for (unsigned int p = 0; p < mesh.Cell2DNumberDoubleProperties(); p++)
        {
            propertyValues[p].resize(mesh.Cell2DTotalNumber());
            for (unsigned int g = 0; g < mesh.Cell2DTotalNumber(); g++)
            {
                propertyValues[p][g] = mesh.Cell2DDoublePropertySize(g, p) == 1 ? mesh.Cell2DDoublePropertyValue(g, p, 0) : 0.0;
            }

            properties[3 + p] = {mesh.Cell2DDoublePropertyId(p),
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(propertyValues[p].size()),
                                 propertyValues[p].data()};
        }

        Gedim::VTKUtilities vtpUtilities;
        vtpUtilities.AddPolygons(mesh.Cell0DsCoordinates(), mesh.Cell2DsVertices(), properties);
        vtpUtilities.Export(cell2DsFolder + "/" + "Cell2Ds_" + fileName + ".vtu");
    }

    // Export Cell3Ds
    if (mesh.Cell3DTotalNumber() > 0)
    {
        vector<double> id(mesh.Cell3DTotalNumber());
        vector<double> marker(mesh.Cell3DTotalNumber());
        vector<double> active(mesh.Cell3DTotalNumber());

        for (unsigned int g = 0; g < mesh.Cell3DTotalNumber(); g++)
        {
            id[g] = g;
            marker[g] = mesh.Cell3DMarker(g);
            active[g] = mesh.Cell3DIsActive(g);
        }

        vector<VTPProperty> properties(3 + mesh.Cell3DNumberDoubleProperties());
        vector<vector<double>> propertyValues(mesh.Cell3DNumberDoubleProperties());

        properties[0] = {"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()};
        properties[1] = {"Marker", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(marker.size()), marker.data()};
        properties[2] = {"Active", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(active.size()), active.data()};

        for (unsigned int p = 0; p < mesh.Cell3DNumberDoubleProperties(); p++)
        {
            propertyValues[p].resize(mesh.Cell3DTotalNumber());
            for (unsigned int g = 0; g < mesh.Cell3DTotalNumber(); g++)
            {
                propertyValues[p][g] = mesh.Cell3DDoublePropertySize(g, p) == 1 ? mesh.Cell3DDoublePropertyValue(g, p, 0) : 0.0;
            }

            properties[3 + p] = {mesh.Cell3DDoublePropertyId(p),
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(propertyValues[p].size()),
                                 propertyValues[p].data()};
        }

        Gedim::VTKUtilities vtpUtilities;
        vtpUtilities.AddPolyhedrons(mesh.Cell0DsCoordinates(), mesh.Cell3DsFacesVertices(), properties);
        vtpUtilities.Export(cell3DsFolder + "/" + "Cell3Ds_" + fileName + ".vtu");
    }
}
// ***************************************************************************
void ExportCell2DToVTU(const Gedim::IMeshDAO &,
                       const unsigned int &cell2DIndex,
                       const Eigen::MatrixXd &cell2DVertices,
                       const vector<Eigen::Matrix3d> &cell2DTriangulations,
                       const double &cell2DArea,
                       const Eigen::Vector3d &cell2DCentroid,
                       const string &exportFolder)
{
    {
        Gedim::VTKUtilities vtpUtilities;

        vector<double> id(1, cell2DIndex);
        vector<double> area(1, cell2DArea);

        // Export cell2D
        vtpUtilities.AddPolygon(
            cell2DVertices,
            {{"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()},
             {"Area", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(area.size()), area.data()}});

        vtpUtilities.Export(exportFolder + "/" + "Cell2D_" + to_string(cell2DIndex) + ".vtu");
    }

    {
        Gedim::VTKUtilities vtpUtilities;

        // Export cell2D triangulation
        for (unsigned int t = 0; t < cell2DTriangulations.size(); t++)
        {
            vector<double> id(1, t);
            vtpUtilities.AddPolygon(cell2DTriangulations[t],
                                    {{"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()}});
        }

        vtpUtilities.Export(exportFolder + "/" + "Cell2D_" + to_string(cell2DIndex) + "_Triangles" + ".vtu");
    }

    {
        Gedim::VTKUtilities vtpUtilities;

        // Export cell2D centroid
        vtpUtilities.AddPoint(cell2DCentroid);

        vtpUtilities.Export(exportFolder + "/" + "Cell2D_" + to_string(cell2DIndex) + "_Centroid" + ".vtu");
    }
}
// ***************************************************************************
void ExportCell3DToVTU(const Gedim::GeometryUtilities &geometryUtilities,
                       const Gedim::IMeshDAO &mesh,
                       const unsigned int &cell3DIndex,
                       const Eigen::MatrixXd &cell3DVertices,
                       const std::vector<Eigen::MatrixXd> &cell3DTetrahedrons,
                       const std::vector<std::vector<Eigen::Matrix3d>> &cell3DFaces3DTriangulations,
                       const double &cell3DVolume,
                       const Eigen::Vector3d &cell3DCentroid,
                       const std::vector<Eigen::Vector3d> &cell3DFacesTranslation,
                       const std::vector<Eigen::Matrix3d> &cell3DFacesRotationMatrix,
                       const std::vector<double> &cell3DFacesArea,
                       const std::vector<Eigen::MatrixXd> &cell3DFaces2DVertices,
                       const std::vector<Eigen::MatrixXd> &cell3DFaces3DVertices,
                       const std::vector<Eigen::VectorXd> &cell3DFacesEdgeLengths,
                       const std::vector<std::vector<bool>> &cell3DFacesEdgeDirections,
                       const std::vector<Eigen::MatrixXd> &cell3DFacesEdges2DTangent,
                       const std::vector<Eigen::MatrixXd> &cell3DFacesEdges2DNormal,
                       const std::vector<Eigen::Vector3d> &cell3DFacesNormals,
                       const std::vector<bool> &cell3DFacesNormalDirections,
                       const std::vector<Eigen::Vector3d> &cell3DFaces2DCentroids,
                       const std::string &exportFolder)
{
    {
        const Gedim::GeometryUtilities::Polyhedron cell3D = Gedim::MeshUtilities::MeshCell3DToPolyhedron(mesh, cell3DIndex);

        {
            Gedim::VTKUtilities vtpUtilities;
            const vector<double> id(1, cell3DIndex);
            const vector<double> volume(1, cell3DVolume);

            // Export cell3D vertices
            vtpUtilities.AddPoints(cell3DVertices);

            vtpUtilities.Export(exportFolder + "/" + "Cell3D_" + to_string(cell3DIndex) + "_Vertices" + ".vtu");
        }

        {
            Gedim::VTKUtilities vtpUtilities;
            const vector<double> id(1, cell3DIndex);
            const vector<double> volume(1, cell3DVolume);

            // Export cell3D edges
            vtpUtilities.AddSegments(cell3D.Vertices, cell3D.Edges);

            vtpUtilities.Export(exportFolder + "/" + "Cell3D_" + to_string(cell3DIndex) + "_Edges" + ".vtu");
        }

        {
            Gedim::VTKUtilities vtpUtilities;
            const vector<double> id(1, cell3DIndex);
            const vector<double> volume(1, cell3DVolume);

            // Export cell3D
            vtpUtilities.AddPolyhedron(
                cell3D.Vertices,
                cell3D.Edges,
                cell3D.Faces,
                {{"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()},
                 {"Volume", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(volume.size()), volume.data()}});

            vtpUtilities.Export(exportFolder + "/" + "Cell3D_" + to_string(cell3DIndex) + ".vtu");
        }
    }

    {
        Gedim::VTKUtilities vtpUtilities;

        // Export cell3D tetrahedrons
        for (unsigned int t = 0; t < cell3DTetrahedrons.size(); t++)
        {
            const Gedim::GeometryUtilities::Polyhedron tetra =
                geometryUtilities.CreateTetrahedronWithVertices(cell3DTetrahedrons[t].col(0),
                                                                cell3DTetrahedrons[t].col(1),
                                                                cell3DTetrahedrons[t].col(2),
                                                                cell3DTetrahedrons[t].col(3));

            const vector<double> id(1, t);
            vtpUtilities.AddPolyhedron(
                tetra.Vertices,
                tetra.Edges,
                tetra.Faces,
                {{"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()}});
        }

        vtpUtilities.Export(exportFolder + "/" + "Cell3D_" + to_string(cell3DIndex) + "_Tetra" + ".vtu");
    }

    {
        Gedim::VTKUtilities vtpUtilities;

        unsigned int totT = 0;
        for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
        {
            const std::vector<Eigen::Matrix3d> &cell2DTriangulations = cell3DFaces3DTriangulations[f];

            // Export cell2D triangulation
            for (unsigned int t = 0; t < cell2DTriangulations.size(); t++)
            {
                const vector<double> face_local_id(1, f);
                const vector<double> face_global_id(1, mesh.Cell3DFace(cell3DIndex, f));
                const vector<double> triangle_local_id(1, t);
                const vector<double> triangle_global_id(1, totT++);
                vtpUtilities.AddPolygon(cell2DTriangulations[t],
                                        {{"face_global_id",
                                          Gedim::VTPProperty::Formats::Cells,
                                          static_cast<unsigned int>(face_global_id.size()),
                                          face_global_id.data()},
                                         {"face_local_id",
                                          Gedim::VTPProperty::Formats::Cells,
                                          static_cast<unsigned int>(face_local_id.size()),
                                          face_local_id.data()},
                                         {"triangle_global_id",
                                          Gedim::VTPProperty::Formats::Cells,
                                          static_cast<unsigned int>(triangle_global_id.size()),
                                          triangle_global_id.data()},
                                         {"triangle_local_id",
                                          Gedim::VTPProperty::Formats::Cells,
                                          static_cast<unsigned int>(triangle_local_id.size()),
                                          triangle_local_id.data()}});
            }
        }

        vtpUtilities.Export(exportFolder + "/" + "Cell3D_" + to_string(cell3DIndex) + "_Faces_Triangles" + ".vtu");
    }

    {
        for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
        {
            const Eigen::MatrixXd &face2DVertices = cell3DFaces2DVertices[f];
            const Eigen::MatrixXd &face3DVertices = cell3DFaces3DVertices[f];
            const Eigen::VectorXd &faceEdgeLengths = cell3DFacesEdgeLengths[f];
            const std::vector<bool> &faceEdgeDirections = cell3DFacesEdgeDirections[f];
            const Eigen::Matrix3d &faceRotation = cell3DFacesRotationMatrix[f];
            const Eigen::Vector3d &faceTranslation = cell3DFacesTranslation[f];
            const Eigen::MatrixXd &faceEdgesNormal = cell3DFacesEdges2DNormal[f];
            const Eigen::MatrixXd &faceEdgesTangent = cell3DFacesEdges2DTangent[f];

            {
                Gedim::VTKUtilities vtpUtilities;
                vtpUtilities.AddPoints(face3DVertices);
                vtpUtilities.Export(exportFolder + "/" + "Cell3D_" + to_string(cell3DIndex) + "_Face_" + to_string(f) +
                                    "_Vertices" + ".vtu");
            }

            {
                Gedim::VTKUtilities vtpUtilities;
                for (unsigned int e = 0; e < faceEdgeDirections.size(); e++)
                {
                    // Export cell2D centroids
                    const vector<double> edge_length(1, faceEdgeLengths[e]);
                    const vector<double> edge_direction(1, faceEdgeDirections[e] ? 1.0 : -1.0);
                    const vector<double> face_local_id(1, f);
                    const vector<double> face_global_id(1, mesh.Cell3DFace(cell3DIndex, f));

                    const Eigen::Vector3d edgeO = face3DVertices.col(e);
                    const Eigen::Vector3d edgeE = face3DVertices.col((e + 1) % faceEdgeDirections.size());

                    vtpUtilities.AddSegment(edgeO,
                                            edgeE,
                                            {{"edge_length",
                                              Gedim::VTPProperty::Formats::Cells,
                                              static_cast<unsigned int>(edge_length.size()),
                                              edge_length.data()},
                                             {"edge_direction",
                                              Gedim::VTPProperty::Formats::Cells,
                                              static_cast<unsigned int>(edge_direction.size()),
                                              edge_direction.data()},
                                             {"face_global_id",
                                              Gedim::VTPProperty::Formats::Cells,
                                              static_cast<unsigned int>(face_global_id.size()),
                                              face_global_id.data()},
                                             {"face_local_id",
                                              Gedim::VTPProperty::Formats::Cells,
                                              static_cast<unsigned int>(face_local_id.size()),
                                              face_local_id.data()}});
                }

                vtpUtilities.Export(exportFolder + "/" + "Cell3D_" + to_string(cell3DIndex) + "_Face_" + to_string(f) +
                                    "_Edges" + ".vtu");
            }

            {
                Gedim::VTKUtilities vtpUtilitiesNormal, vtpUtilitiesTangent;
                for (unsigned int e = 0; e < faceEdgeDirections.size(); e++)
                {
                    // Export cell2D centroids
                    const vector<double> edge_direction(1, faceEdgeDirections[e] ? 1.0 : -1.0);
                    const vector<double> face_local_id(1, f);
                    const vector<double> face_global_id(1, mesh.Cell3DFace(cell3DIndex, f));

                    const double sizeNormal =
                        std::max(std::max(std::sqrt(cell3DFacesArea[f]), std::cbrt(cell3DVolume)), faceEdgeLengths[e]);

                    const Eigen::Vector3d edgeO2D = face2DVertices.col(e);
                    const Eigen::Vector3d edgeE2D = face2DVertices.col((e + 1) % faceEdgeDirections.size());

                    const Eigen::Vector3d edgeNormal2DO = 0.5 * (edgeO2D + edgeE2D);
                    const Eigen::Vector3d edgeNormal2DE = edgeNormal2DO + sizeNormal * faceEdgesNormal.col(e);

                    const Eigen::Vector3d edgeTangent2DO = edgeO2D;
                    const Eigen::Vector3d edgeTangent2DE = edgeO2D + faceEdgesTangent.col(e);

                    const Eigen::Vector3d edgeNormal3DO =
                        geometryUtilities.RotatePointsFrom2DTo3D(edgeNormal2DO, faceRotation, faceTranslation);
                    const Eigen::Vector3d edgeNormal3DE =
                        geometryUtilities.RotatePointsFrom2DTo3D(edgeNormal2DE, faceRotation, faceTranslation);

                    const Eigen::Vector3d edgeTangent3DO =
                        geometryUtilities.RotatePointsFrom2DTo3D(edgeTangent2DO, faceRotation, faceTranslation);
                    const Eigen::Vector3d edgeTangent3DE =
                        geometryUtilities.RotatePointsFrom2DTo3D(edgeTangent2DE, faceRotation, faceTranslation);

                    vtpUtilitiesNormal.AddSegment(edgeTangent3DO,
                                                  edgeTangent3DE,
                                                  {{"edge_direction",
                                                    Gedim::VTPProperty::Formats::Cells,
                                                    static_cast<unsigned int>(edge_direction.size()),
                                                    edge_direction.data()}});
                    vtpUtilitiesTangent.AddSegment(edgeNormal3DO, edgeNormal3DE);
                }

                vtpUtilitiesNormal.Export(exportFolder + "/" + "Cell3D_" + to_string(cell3DIndex) + "_Face_" +
                                          to_string(f) + "_EdgeTangents" + ".vtu");
                vtpUtilitiesTangent.Export(exportFolder + "/" + "Cell3D_" + to_string(cell3DIndex) + "_Face_" +
                                           to_string(f) + "_EdgeNormals" + ".vtu");
            }
        }
    }

    {
        Gedim::VTKUtilities vtpUtilities;

        for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
        {
            const Eigen::Vector3d &faceCentroid = geometryUtilities.RotatePointsFrom2DTo3D(cell3DFaces2DCentroids[f],
                                                                                           cell3DFacesRotationMatrix[f],
                                                                                           cell3DFacesTranslation[f]);

            // Export cell2D centroids
            const vector<double> face_local_id(1, f);
            const vector<double> face_global_id(1, mesh.Cell3DFace(cell3DIndex, f));
            vtpUtilities.AddPoint(faceCentroid,
                                  {{"face_global_id",
                                    Gedim::VTPProperty::Formats::Cells,
                                    static_cast<unsigned int>(face_global_id.size()),
                                    face_global_id.data()},
                                   {"face_local_id",
                                    Gedim::VTPProperty::Formats::Cells,
                                    static_cast<unsigned int>(face_local_id.size()),
                                    face_local_id.data()}});
        }

        vtpUtilities.Export(exportFolder + "/" + "Cell3D_" + to_string(cell3DIndex) + "_Faces_Centroid" + ".vtu");
    }

    {
        Gedim::VTKUtilities vtpUtilities;

        for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
        {
            const Eigen::Vector3d &faceCentroid = geometryUtilities.RotatePointsFrom2DTo3D(cell3DFaces2DCentroids[f],
                                                                                           cell3DFacesRotationMatrix[f],
                                                                                           cell3DFacesTranslation[f]);
            const double sizeNormal = std::max(std::sqrt(cell3DFacesArea[f]), std::cbrt(cell3DVolume));
            const Eigen::Vector3d &faceNormal = cell3DFacesNormals[f];
            const double faceDirection = cell3DFacesNormalDirections[f] ? 1.0 : -1.0;

            // Export cell2D centroids
            const vector<double> face_local_id(1, f);
            const vector<double> face_global_id(1, mesh.Cell3DFace(cell3DIndex, f));
            vtpUtilities.AddSegment(faceCentroid,
                                    faceCentroid + faceDirection * sizeNormal * faceNormal,
                                    {{"face_global_id",
                                      Gedim::VTPProperty::Formats::Cells,
                                      static_cast<unsigned int>(face_global_id.size()),
                                      face_global_id.data()},
                                     {"face_local_id",
                                      Gedim::VTPProperty::Formats::Cells,
                                      static_cast<unsigned int>(face_local_id.size()),
                                      face_local_id.data()}});
        }

        vtpUtilities.Export(exportFolder + "/" + "Cell3D_" + to_string(cell3DIndex) + "_Faces_Normal" + ".vtu");
    }

    {
        Gedim::VTKUtilities vtpUtilities;

        // Export cell3D centroid
        vtpUtilities.AddPoint(cell3DCentroid);

        vtpUtilities.Export(exportFolder + "/" + "Cell3D_" + to_string(cell3DIndex) + "_Centroid" + ".vtu");
    }
}
} // namespace MeshUtilities
} // namespace External
} // namespace Gedim

// ***************************************************************************
namespace Gedim
{
namespace External
{
namespace GeometryUtilities
{
void ExportPolyhedronToVTU(const Eigen::MatrixXd &polyhedronVertices,
                           const Eigen::MatrixXi &polyhedronEdges,
                           const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                           const std::string &exportFolder)
{
    vector<vector<unsigned int>> facesVertices(polyhedronFaces.size());
    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
        facesVertices[f].resize(polyhedronFaces[f].cols());
        for (unsigned int v = 0; v < polyhedronFaces[f].cols(); v++)
            facesVertices[f][v] = polyhedronFaces[f](0, v);
    }

    {
        // export vertices
        VTKUtilities exporter;
        vector<double> verticesIndex(polyhedronVertices.cols());
        for (unsigned int v = 0; v < polyhedronVertices.cols(); v++)
            verticesIndex[v] = v;

        exporter.AddPoints(polyhedronVertices,
                           {{"Index",
                             Gedim::VTPProperty::Formats::Cells,
                             static_cast<unsigned int>(verticesIndex.size()),
                             verticesIndex.data()}});
        exporter.Export(exportFolder + "/Polyhedron_Vertices.vtu");
    }

    {
        // export edges
        VTKUtilities exporter;
        vector<double> edgesIndex(polyhedronEdges.cols());
        for (unsigned int e = 0; e < polyhedronEdges.cols(); e++)
            edgesIndex[e] = e;

        exporter.AddSegments(
            polyhedronVertices,
            polyhedronEdges,
            {{"Index", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(edgesIndex.size()), edgesIndex.data()}});
        exporter.Export(exportFolder + "/Polyhedron_Edges.vtu");
    }

    {
        // export faces
        VTKUtilities exporter;
        vector<double> facesIndex(polyhedronFaces.size());
        for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
            facesIndex[f] = f;

        exporter.AddPolygons(
            polyhedronVertices,
            facesVertices,
            {{"Index", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(facesIndex.size()), facesIndex.data()}});
        exporter.Export(exportFolder + "/Polyhedron_Faces.vtu");
    }

    {
        // export polyhedron
        VTKUtilities exporter;
        exporter.AddPolyhedron(polyhedronVertices, facesVertices);
        exporter.Export(exportFolder + "/Polyhedron.vtu");
    }
}
// ***************************************************************************
void ExportPolyhedronToVTU(const Gedim::GeometryUtilities &geometryUtilities,
                           const unsigned int &index,
                           const Eigen::MatrixXd &polyhedronVertices,
                           const Eigen::MatrixXi &polyhedronEdges,
                           const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                           const std::vector<Eigen::MatrixXd> &polyhedronTetra,
                           const double &polyhedronVolume,
                           const Eigen::Vector3d &polyhedronCentroid,
                           const std::vector<Eigen::MatrixXd> &polyhedronFaces3DVertices,
                           const std::vector<double> &polyhedronFacesArea,
                           const std::vector<Eigen::Vector3d> &polyhedronFaces2DCentroid,
                           const std::vector<Eigen::Vector3d> &polyhedronFacesTranslation,
                           const std::vector<Eigen::Matrix3d> &polyhedronFacesRotationMatrix,
                           const std::vector<std::vector<Eigen::Matrix3d>> &polyhedronFaces3DTriangles,
                           const std::vector<Eigen::Vector3d> &polyhedronFaces3DInternalPoint,
                           const std::vector<Eigen::Vector3d> &polyhedronFaces3DNormal,
                           const std::vector<bool> &polyhedronFaces3DNormalDirection,
                           const std::string &exportFolder)
{
    {
        Gedim::VTKUtilities exporter;

        vector<double> id(1, index);
        vector<double> volume(1, polyhedronVolume);

        // Export Polyhedron
        exporter.AddPolyhedron(
            polyhedronVertices,
            polyhedronEdges,
            polyhedronFaces,
            {{"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()},
             {"Volume", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(volume.size()), volume.data()}});

        exporter.Export(exportFolder + "/" + "Polyhedron.vtu");
    }

    {
        Gedim::VTKUtilities exporter;

        // Export Polyhedron tetra
        for (unsigned int t = 0; t < polyhedronTetra.size(); t++)
        {
            vector<double> id(1, t);

            const Gedim::GeometryUtilities::Polyhedron polyhedron =
                geometryUtilities.CreateTetrahedronWithVertices(polyhedronTetra[t].col(0),
                                                                polyhedronTetra[t].col(1),
                                                                polyhedronTetra[t].col(2),
                                                                polyhedronTetra[t].col(3));

            exporter.AddPolyhedron(polyhedron.Vertices,
                                   polyhedron.Edges,
                                   polyhedron.Faces,
                                   {{"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()}});
        }

        exporter.Export(exportFolder + "/" + "Polyhedron_Tetra.vtu");
    }

    {
        Gedim::VTKUtilities exporter;

        for (unsigned int f = 0; f < polyhedronFaces3DVertices.size(); f++)
        {
            vector<double> id(1, f);
            vector<double> area(1, polyhedronFacesArea[f]);

            // Export faces
            exporter.AddPolygon(
                polyhedronFaces3DVertices[f],
                {{"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()},
                 {"Area", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(area.size()), area.data()}});
        }

        exporter.Export(exportFolder + "/" + "Polyhedron_Faces.vtu");
    }

    {
        Gedim::VTKUtilities exporter;

        // Export faces triangles
        unsigned int numFaceTriangles = 0;
        for (unsigned int f = 0; f < polyhedronFaces3DTriangles.size(); f++)
        {
            vector<double> face(1, f);

            for (unsigned int t = 0; t < polyhedronFaces3DTriangles[f].size(); t++)
            {
                vector<double> id(1, numFaceTriangles++);
                exporter.AddPolygon(
                    polyhedronFaces3DTriangles[f][t],
                    {{"Face", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(face.size()), face.data()},
                     {"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()}});
            }
        }

        exporter.Export(exportFolder + "/" + "Polyhedron_FacesTriangles.vtu");
    }

    {
        Gedim::VTKUtilities exporter;

        // Export faces normal
        for (unsigned int f = 0; f < polyhedronFaces3DNormal.size(); f++)
        {
            vector<double> face(1, f);
            vector<double> normalDirection(1, polyhedronFaces3DNormalDirection[f] ? 1.0 : -1.0);

            exporter.AddSegment(polyhedronFaces3DInternalPoint[f],
                                polyhedronFaces3DInternalPoint[f] + normalDirection[0] * polyhedronFaces3DNormal[f],
                                {{"Face", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(face.size()), face.data()},
                                 {"NormalDirection",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(normalDirection.size()),
                                  normalDirection.data()}});
        }

        exporter.Export(exportFolder + "/" + "Polyhedron_FacesNormal.vtu");
    }

    {
        Gedim::VTKUtilities exporter;

        // Export Polyhedron centroid
        exporter.AddPoint(polyhedronCentroid);

        exporter.Export(exportFolder + "/" + "Polyhedron_Centroid.vtu");
    }

    {
        Gedim::VTKUtilities exporter;

        // Export faces centroids
        for (unsigned int f = 0; f < polyhedronFaces2DCentroid.size(); f++)
        {
            vector<double> face(1, f);

            const Eigen::Vector3d rotatedCentroid = geometryUtilities.RotatePointsFrom2DTo3D(polyhedronFaces2DCentroid[f],
                                                                                             polyhedronFacesRotationMatrix[f],
                                                                                             polyhedronFacesTranslation[f]);
            exporter.AddPoint(rotatedCentroid,
                              {{"Face", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(face.size()), face.data()}});
        }

        exporter.Export(exportFolder + "/" + "Polyhedron_FacesCentroid.vtu");
    }
}

// ***************************************************************************
void ExportPolygonToVTU(const unsigned int &index,
                        const Eigen::MatrixXd &polygon,
                        const std::vector<Eigen::Matrix3d> &polygon_triangles,
                        const double &polygon_volume,
                        const Eigen::Vector3d &polygon_centroid,
                        const Eigen::MatrixXd &polygon_edges_centroid,
                        const Eigen::MatrixXd &polygon_edges_normal,
                        const std::vector<bool> &polygon_edges_normal_direction,
                        const std::string &exportFolder)
{
    {
        Gedim::VTKUtilities exporter;

        vector<double> id(1, index);
        vector<double> volume(1, polygon_volume);

        // Export Polyhedron
        exporter.AddPolygon(
            polygon,
            {{"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()},
             {"Area", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(volume.size()), volume.data()}});

        exporter.Export(exportFolder + "/" + "Polygon.vtu");
    }

    {
        Gedim::VTKUtilities exporter;

        // Export polygon triangles
        for (unsigned int t = 0; t < polygon_triangles.size(); t++)
        {
            vector<double> id(1, t);

            const Eigen::MatrixXd triangle = polygon_triangles.at(t);

            exporter.AddPolygon(triangle,
                                {{"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()}});
        }

        exporter.Export(exportFolder + "/" + "Polygon_Triangles.vtu");
    }

    {
        Gedim::VTKUtilities exporter;

        // Export edges normal
        for (unsigned int e = 0; e < polygon.cols(); e++)
        {
            vector<double> edge(1, e);
            vector<double> normalDirection(1, polygon_edges_normal_direction[e] ? 1.0 : -1.0);

            const Eigen::Vector3d normal_origin = polygon_edges_centroid.col(e);
            const Eigen::Vector3d normal_end = polygon_edges_centroid.col(e) + normalDirection[0] * polygon_edges_normal.col(e);

            exporter.AddSegment(normal_origin,
                                normal_end,
                                {{"Face", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(edge.size()), edge.data()},
                                 {"NormalDirection",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(normalDirection.size()),
                                  normalDirection.data()}});
        }

        exporter.Export(exportFolder + "/" + "Polygon_Edges_Normal.vtu");
    }

    {
        Gedim::VTKUtilities exporter;

        // Export centroid
        exporter.AddPoint(polygon_centroid);

        exporter.Export(exportFolder + "/" + "Polygon_Centroid.vtu");
    }
}
// ***************************************************************************
void ExportPlanteToVTU(const Gedim::GeometryUtilities &geometryUtilities,
                       const Eigen::Vector3d &plane_origin,
                       const Eigen::Vector3d &plane_normal,
                       const Eigen::Vector3d &plane_translation,
                       const Eigen::Matrix3d &plane_rotation_matrix,
                       const string &exportFolder)
{
    Eigen::Matrix3d plane_triangle_2D;
    plane_triangle_2D.col(0) << 1.0, 0.0, 0.0;
    plane_triangle_2D.col(1) << -0.5, +sqrt(3.0) / 2.0, 0.0;
    plane_triangle_2D.col(2) << -0.5, -sqrt(3.0) / 2.0, 0.0;

    const auto plane_triangle_3D =
        geometryUtilities.RotatePointsFrom2DTo3D(plane_triangle_2D, plane_rotation_matrix, plane_translation);

    {
        Gedim::VTKUtilities exporter;
        exporter.AddPolygon(plane_triangle_3D);
        exporter.Export(exportFolder + "/" + "plane.vtu");
    }

    {
        Gedim::VTKUtilities exporter;
        exporter.AddPoint(plane_origin);
        exporter.Export(exportFolder + "/" + "plane_origin.vtu");
    }

    {
        Gedim::VTKUtilities exporter;

        // Export Polyhedron
        exporter.AddSegment(plane_origin, plane_origin + plane_normal);
        exporter.Export(exportFolder + "/" + "plane_normal.vtu");
    }
}
// ***************************************************************************
} // namespace GeometryUtilities

} // namespace External

} // namespace Gedim
