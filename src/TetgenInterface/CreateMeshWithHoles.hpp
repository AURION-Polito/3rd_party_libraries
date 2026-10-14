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

#ifndef __CreateMeshWithHoles_hpp_
#define __CreateMeshWithHoles_hpp_

#include "GeometryUtilities.hpp"
#include "GraphUtilities.hpp"
#include "IMeshDAO.hpp"

// #define EXPORT_TEST_VTU

using namespace std;
using namespace Eigen;

namespace Gedim
{
namespace External
{
namespace MeshUtilities
{
struct AgglomerateMeshByFunctionResult final
{
    std::vector<unsigned int> Cell3Ds_mark;
    std::vector<std::vector<Eigen::MatrixXd>> Cell3Ds_tetra_vertices;
    std::vector<std::vector<Eigen::Matrix3d>> Cell2Ds_triangles_vertices;
};

struct SplitPolyhedronByFunctionResult final
{
    enum struct Types
    {
        NoMarks = 0,
        AllMarked = 1,
        Splitted = 2
    };

    struct Vertex final
    {
        Eigen::Vector3d Coordinates;
        unsigned int Original_vertex_index;
        unsigned int Marker;
        bool IsNew;
    };

    struct Edge final
    {
        std::array<unsigned int, 2> Extremes;
        unsigned int Original_edge_index;
        unsigned int Marker;
        bool IsNew;
    };

    struct Face final
    {
        Eigen::MatrixXi Extremes;
        unsigned int Original_face_index;
        unsigned int Marker;
        bool IsNew;
        std::vector<Eigen::Matrix3d> TriaglesVertices;
    };

    struct Polyhedron final
    {
        std::vector<unsigned int> Vertices;
        std::vector<unsigned int> Edges;
        std::vector<unsigned int> Faces;
        std::vector<Eigen::MatrixXd> TetrahedronVertices;
        unsigned int Mark;
    };

    Types Type;
    std::vector<unsigned int> Original_vertex_to_split_vertex;
    std::vector<std::vector<unsigned int>> Original_edge_to_split_edges;
    std::vector<std::vector<unsigned int>> Original_face_to_split_faces;
    std::vector<Vertex> Split_vertices;
    std::vector<Edge> Split_edges;
    std::vector<Face> Split_faces;
    std::vector<Polyhedron> Split_polyhedrons;
};

Eigen::MatrixXd FilterConstrainedPoints(const Gedim::GeometryUtilities &geometryUtilities,
                                        const Eigen::MatrixXd &polyhedron_vertices,
                                        const Eigen::MatrixXd &polyhedron_additional_constrained_points);

bool AgglomerateMeshCellByFunction(const Gedim::GeometryUtilities &geometryUtilities,
                                   const Gedim::GraphUtilities &graphUtilities,
                                   Gedim::IMeshDAO &mesh,
                                   const Eigen::MatrixXd &cell3D_Vertices,
                                   const Eigen::MatrixXi &cell3D_Edges,
                                   const std::vector<Eigen::MatrixXi> &cell3D_Faces,
                                   std::vector<Eigen::MatrixXd> &cell3Ds_additional_constrained_points,
                                   std::vector<Eigen::Vector3d> &cell3Ds_centroid,
                                   std::vector<double> &cell3Ds_volume,
                                   std::vector<unsigned int> &cell3Ds_mark,
                                   std::vector<std::vector<Eigen::MatrixXd>> &cell3Ds_tetra_vertices,
                                   std::vector<std::vector<Eigen::Matrix3d>> &cell2Ds_triangles_vertices,
                                   const unsigned int cell3D_index,
                                   const std::function<Eigen::VectorXi(const Eigen::MatrixXd &)> &marking_function,
                                   const double &resolution);

SplitPolyhedronByFunctionResult SplitPolyhedronByFunction(const Gedim::GeometryUtilities &geometryUtilities,
                                                          const Gedim::GraphUtilities &graphUtilities,
                                                          const Eigen::MatrixXd &polyhedron_vertices,
                                                          const Eigen::MatrixXi &polyhedron_edges,
                                                          const std::vector<Eigen::MatrixXi> &polyhedron_faces,
                                                          const Eigen::MatrixXd &polyhedron_additional_constrained_points,
                                                          const double &polyhedron_volume,
                                                          const std::function<Eigen::VectorXi(const Eigen::MatrixXd &)> &marking_function,
                                                          const double &resolution,
                                                          const unsigned int polyhedron_index);

AgglomerateMeshByFunctionResult AgglomerateMeshByFunction(const Gedim::GeometryUtilities &geometryUtilities,
                                                          const Gedim::GraphUtilities &graphUtilities,
                                                          Gedim::MeshMatrices &mesh_data,
                                                          const std::function<Eigen::VectorXi(const Eigen::MatrixXd &)> &marking_function,
                                                          const double &resolution,
                                                          const std::vector<Eigen::MatrixXd> &cell3Ds_constrained_points = {});

void CreateMeshWithHoles(const Gedim::GeometryUtilities &geometry_utilities,
                         const std::vector<Gedim::GeometryUtilities::Polyhedron> &domains,
                         Gedim::MeshMatrices &full_mesh_data,
                         Gedim::IMeshDAO &full_mesh,
                         IMeshDAO &mesh);

} // namespace MeshUtilities
} // namespace External
} // namespace Gedim

#endif
