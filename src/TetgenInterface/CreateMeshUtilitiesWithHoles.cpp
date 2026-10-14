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

#include "CreateMeshUtilitiesWithHoles.hpp"

#include "CommonUtilities.hpp"
#include "Eigen_Utilities.hpp"
#include "GraphUtilities.hpp"
#include "MeshUtilities.hpp"
#include "TetgenInterface.hpp"

#ifdef EXPORT_TEST_VTU
#include "VTKUtilities.hpp"
#include "VtkMeshExporter.hpp"
#endif

using namespace std;
using namespace Eigen;

namespace Gedim
{
namespace External
{
namespace MeshUtilities
{
// ***************************************************************************
Eigen::MatrixXd FilterConstrainedPoints(const Gedim::GeometryUtilities &geometryUtilities,
                                        const Eigen::MatrixXd &polyhedron_vertices,
                                        const Eigen::MatrixXd &polyhedron_additional_constrained_points)
{
    std::list<Eigen::Vector3d> filtered_points;

    for (unsigned int p = 0; p < polyhedron_additional_constrained_points.cols(); p++)
    {
        const Eigen::Vector3d point = polyhedron_additional_constrained_points.col(p);

        const auto find_point_index = geometryUtilities.FindPointInPoints(polyhedron_vertices, point);

        if (find_point_index.empty())
            filtered_points.push_back(point);
    }

    return Gedim::Eigen_Utilities::CoordinatesToMatrix(filtered_points);
}
// ***************************************************************************
SplitPolyhedronByFunctionResult SplitPolyhedronByFunction(const Gedim::GeometryUtilities &geometryUtilities,
                                                          const Gedim::GraphUtilities &graphUtilities,
                                                          const Eigen::MatrixXd &polyhedron_vertices,
                                                          const Eigen::MatrixXi &polyhedron_edges,
                                                          const std::vector<Eigen::MatrixXi> &polyhedron_faces,
                                                          const Eigen::MatrixXd &polyhedron_additional_constrained_points,
                                                          const double &polyhedron_volume,
                                                          const std::function<Eigen::VectorXi(const Eigen::MatrixXd &)> &marking_function,
                                                          const double &resolution,
                                                          const unsigned int polyhedron_index)
{
    Gedim::Utilities::Unused(polyhedron_index);

#ifdef EXPORT_TEST_VTU
    std::string mesh_3D_export_folder = "./TEST/Polyhedron_" + std::to_string(polyhedron_index);
    Gedim::Output::CreateFolder(mesh_3D_export_folder);
#endif

    const unsigned int num_polyhedron_vertices = polyhedron_vertices.cols();
    const unsigned int num_polyhedron_edges = polyhedron_edges.cols();
    const unsigned int num_polyhedron_faces = polyhedron_faces.size();

    const double mesh_max_volume = polyhedron_volume * resolution;

#ifdef EXPORT_TEST_VTU
    {
        Gedim::External::GeometryUtilities::ExportPolyhedronToVTU(polyhedron_vertices, polyhedron_edges, polyhedron_faces, mesh_3D_export_folder);
    }
#endif

    Gedim::MeshMatrices original_mesh_data;
    Gedim::MeshMatricesDAO original_mesh(original_mesh_data);

    if (geometryUtilities.IsValueZero(mesh_max_volume, geometryUtilities.Tolerance3D()))
    {
        const Eigen::MatrixXd filtered_constrained_points =
            FilterConstrainedPoints(geometryUtilities, polyhedron_vertices, polyhedron_additional_constrained_points);

        const unsigned int num_additional_constrained_points = filtered_constrained_points.cols();
        Eigen::MatrixXd constraneid_points(3, num_polyhedron_vertices + num_additional_constrained_points);
        constraneid_points.block(0, 0, 3, num_polyhedron_vertices) = polyhedron_vertices;
        constraneid_points.block(0, num_polyhedron_vertices, 3, num_additional_constrained_points) = filtered_constrained_points;

        std::vector<unsigned int> constrained_points_marker(num_polyhedron_vertices + num_additional_constrained_points);
        for (unsigned int v = 0; v < num_polyhedron_vertices; v++)
            constrained_points_marker[v] = v + 1;
        for (unsigned int c_p = 0; c_p < num_additional_constrained_points; ++c_p)
            constrained_points_marker[num_polyhedron_vertices + c_p] = 0;

        Gedim::External::MeshUtilities::CreateDelaunayMesh3D(constraneid_points, constrained_points_marker, original_mesh);

        const auto polyhedron_edges_length = geometryUtilities.PolyhedronEdgesLength(polyhedron_vertices, polyhedron_edges);
        const auto polyhedron_edges_tangent = geometryUtilities.PolyhedronEdgeTangents(polyhedron_vertices, polyhedron_edges);
        const auto polyhedron_faces_vertices = geometryUtilities.PolyhedronFaceVertices(polyhedron_vertices, polyhedron_faces);
        const auto polyhedron_faces_normal = geometryUtilities.PolyhedronFaceNormals(polyhedron_faces_vertices);
        const auto polyhedron_faces_translation = geometryUtilities.PolyhedronFaceTranslations(polyhedron_faces_vertices);
        const auto polyhedron_faces_rotation_matrix =
            geometryUtilities.PolyhedronFaceRotationMatrices(polyhedron_faces_vertices, polyhedron_faces_normal, polyhedron_faces_translation);
        const auto polyhedron_faces_vertices_2D =
            geometryUtilities.PolyhedronFaceRotatedVertices(polyhedron_faces_vertices, polyhedron_faces_translation, polyhedron_faces_rotation_matrix);

        std::vector<unsigned int> vertices_marker(num_polyhedron_vertices);
        for (unsigned int v = 0; v < num_polyhedron_vertices; v++)
            vertices_marker[v] = v + 1;

        std::vector<unsigned int> edges_marker(num_polyhedron_edges);
        for (unsigned int e = 0; e < num_polyhedron_edges; e++)
            edges_marker[e] = num_polyhedron_vertices + e + 1;

        std::vector<unsigned int> faces_marker(num_polyhedron_faces);
        for (unsigned int f = 0; f < num_polyhedron_faces; f++)
            faces_marker[f] = num_polyhedron_vertices + num_polyhedron_edges + f + 1;

        std::vector<Eigen::Vector3d> cell1Ds_centroid(original_mesh.Cell1DTotalNumber());
        for (unsigned int c = 0; c < original_mesh.Cell1DTotalNumber(); ++c)
        {
            cell1Ds_centroid[c] = geometryUtilities.SegmentBarycenter(original_mesh.Cell1DOriginCoordinates(c),
                                                                      original_mesh.Cell1DEndCoordinates(c));
        }

        std::vector<Eigen::Vector3d> cell2Ds_centroid(original_mesh.Cell2DTotalNumber());
        for (unsigned int c = 0; c < original_mesh.Cell2DTotalNumber(); ++c)
            cell2Ds_centroid[c] = geometryUtilities.SimplexBarycenter(original_mesh.Cell2DVerticesCoordinates(c));

        Gedim::MeshUtilities::ChangePolyhedronMeshMarkers(geometryUtilities,
                                                          polyhedron_vertices,
                                                          polyhedron_edges,
                                                          polyhedron_faces,
                                                          polyhedron_edges_tangent,
                                                          polyhedron_edges_length,
                                                          polyhedron_faces_normal,
                                                          polyhedron_faces_vertices,
                                                          polyhedron_faces_vertices_2D,
                                                          polyhedron_faces_translation,
                                                          polyhedron_faces_rotation_matrix,
                                                          vertices_marker,
                                                          edges_marker,
                                                          faces_marker,
                                                          cell1Ds_centroid,
                                                          cell2Ds_centroid,
                                                          original_mesh);
    }
    else
    {
        Gedim::External::MeshUtilities::CreateTetrahedralMesh(polyhedron_vertices, polyhedron_edges, polyhedron_faces, mesh_max_volume, original_mesh);
    }

    Gedim::MeshMatrices local_mesh_data = original_mesh_data;
    Gedim::MeshMatricesDAO local_mesh(local_mesh_data);

#ifdef EXPORT_TEST_VTU
    {
        Gedim::External::MeshUtilities::ExportMeshToVTU(local_mesh, mesh_3D_export_folder, "LocalMesh3D");
    }
#endif

    Gedim::MeshUtilities::ComputeCell0DCell3DNeighbours(local_mesh);
    Gedim::MeshUtilities::ComputeCell0DCell1DNeighbours(local_mesh);
    Gedim::MeshUtilities::ComputeCell0DCell2DNeighbours(local_mesh);
    Gedim::MeshUtilities::ComputeCell1DCell2DNeighbours(local_mesh);
    Gedim::MeshUtilities::ComputeCell2DCell3DNeighbours(local_mesh);

    std::vector<Eigen::MatrixXd> cell3Ds_mark_points(local_mesh.Cell3DTotalNumber());

    if (polyhedron_additional_constrained_points.size() > 0)
    {
        for (unsigned int c = 0; c < local_mesh.Cell3DTotalNumber(); ++c)
            cell3Ds_mark_points[c] = geometryUtilities.SimplexBarycenter(local_mesh.Cell3DVerticesCoordinates(c));
    }
    else
    {
        for (unsigned int c = 0; c < local_mesh.Cell3DTotalNumber(); ++c)
        {
            const unsigned int num_cell_vertices = local_mesh.Cell3DNumberVertices(c);
            const auto cell3D_vertices = local_mesh.Cell3DVerticesCoordinates(c);
            cell3Ds_mark_points[c].resize(3, num_cell_vertices + 1);
            cell3Ds_mark_points[c].block(0, 0, 3, num_cell_vertices) = cell3D_vertices;
            cell3Ds_mark_points[c].col(num_cell_vertices) << geometryUtilities.SimplexBarycenter(cell3D_vertices);
        }
    }

    auto cell3Ds_mark = Gedim::MeshUtilities::MarkCells(marking_function, cell3Ds_mark_points, 0);

    {
        unsigned int num_marked_cells = 0;
        for (unsigned int c = 0; c < local_mesh.Cell3DTotalNumber(); ++c)
        {
            if (cell3Ds_mark.at(c) > 0)
                num_marked_cells++;
        }

        if (num_marked_cells == 0)
            return {SplitPolyhedronByFunctionResult::Types::NoMarks, {}, {}, {}, {}, {}, {}, {}};

        if (num_marked_cells == local_mesh.Cell3DTotalNumber())
            return {SplitPolyhedronByFunctionResult::Types::AllMarked, {}, {}, {}, {}, {}, {}, {}};
    }

#ifdef EXPORT_TEST_VTU
    {
        {
            Gedim::VTKUtilities exporter;
            std::vector<double> marks(cell3Ds_mark.begin(), cell3Ds_mark.end());

            exporter.AddPolyhedrons(
                local_mesh.Cell0DsCoordinates(),
                local_mesh.Cell3DsFacesVertices(),
                {{"Mark", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(marks.size()), marks.data()}});
            exporter.Export(mesh_3D_export_folder + "/cell3Ds_mark_original.vtu");
        }
    }
#endif

    {
        std::array<std::list<unsigned int>, 2> cell3Ds_mark_lists;
        for (unsigned int c = 0; c < local_mesh.Cell3DTotalNumber(); ++c)
        {
            if (cell3Ds_mark.at(c) == 0)
                cell3Ds_mark_lists[0].push_back(c);
            else
                cell3Ds_mark_lists[1].push_back(c);
        }

        std::vector<unsigned int> cell3Ds_index(cell3Ds_mark_lists.at(0).begin(), cell3Ds_mark_lists.at(0).end());
        std::unordered_map<unsigned int, unsigned int> cell3Ds;
        for (const unsigned int c_3D : cell3Ds_mark_lists.at(0))
            cell3Ds.insert(std::make_pair(c_3D, cell3Ds.size()));

        std::vector<std::list<unsigned int>> cell3Ds_adj_lst(cell3Ds.size());

        for (const auto &c_3D : cell3Ds)
        {
            for (unsigned int f = 0; f < local_mesh.Cell3DNumberFaces(c_3D.first); f++)
            {
                const unsigned int c_2D = local_mesh.Cell3DFace(c_3D.first, f);

                for (unsigned int n = 0; n < local_mesh.Cell2DNumberNeighbourCell3D(c_2D); n++)
                {
                    if (!local_mesh.Cell2DHasNeighbourCell3D(c_2D, n))
                        continue;

                    const unsigned int c_3D_n = local_mesh.Cell2DNeighbourCell3D(c_2D, n);

                    if (c_3D_n == c_3D.first)
                        continue;

                    const auto c_3D_n_pos = cell3Ds.find(c_3D_n);
                    if (c_3D_n_pos == cell3Ds.end())
                        continue;

                    cell3Ds_adj_lst[c_3D.second].push_back(c_3D_n_pos->second);
                }
            }
        }

        std::vector<std::vector<unsigned int>> cell3Ds_adj(cell3Ds.size());
        for (unsigned int c = 0; c < cell3Ds_adj_lst.size(); c++)
        {
            if (cell3Ds_adj_lst.at(c).empty())
                continue;

            cell3Ds_adj[c] = std::vector<unsigned int>(cell3Ds_adj_lst.at(c).begin(), cell3Ds_adj_lst.at(c).end());
        }

        const auto connected_components = graphUtilities.ComputeStronglyConnectedComponents(cell3Ds_adj);

        if (connected_components.size() > 2)
        {
            std::cout << "Connected components - case not managed yet" << std::endl;
            return {SplitPolyhedronByFunctionResult::Types::NoMarks, {}, {}, {}, {}, {}, {}, {}};
        }

        if (connected_components.size() == 2)
        {
            const auto &min_connected = connected_components.at(0).size() < connected_components.at(1).size()
                                            ? connected_components.at(0)
                                            : connected_components.at(1);
            const auto &max_connected = connected_components.at(0).size() < connected_components.at(1).size()
                                            ? connected_components.at(1)
                                            : connected_components.at(0);

            cell3Ds_mark_lists[0].clear();
            for (const unsigned int cell3D_pos : max_connected)
            {
                const unsigned int c_3D = cell3Ds_index.at(cell3D_pos);
                cell3Ds_mark_lists[0].push_back(c_3D);
                cell3Ds_mark[c_3D] = 0;
            }

            for (const unsigned int cell3D_pos : min_connected)
            {
                const unsigned int c_3D = cell3Ds_index.at(cell3D_pos);
                cell3Ds_mark_lists[1].push_back(c_3D);
                cell3Ds_mark[c_3D] = 1;
            }
        }
    }

#ifdef EXPORT_TEST_VTU
    {
        {
            Gedim::VTKUtilities exporter;
            std::vector<double> marks(cell3Ds_mark.begin(), cell3Ds_mark.end());

            exporter.AddPolyhedrons(
                local_mesh.Cell0DsCoordinates(),
                local_mesh.Cell3DsFacesVertices(),
                {{"Mark", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(marks.size()), marks.data()}});
            exporter.Export(mesh_3D_export_folder + "/cell3Ds_mark_connected.vtu");
        }
    }
#endif

    {
        std::set<unsigned int> border_points;
        for (unsigned int c_0D = 0; c_0D < local_mesh.Cell0DTotalNumber(); c_0D++)
        {
            if (!local_mesh.Cell0DIsActive(c_0D))
                continue;

            const unsigned int cell0D_marker = local_mesh.Cell0DMarker(c_0D);

            if (cell0D_marker < num_polyhedron_vertices + num_polyhedron_edges + 1)
                continue;

            if (cell0D_marker > num_polyhedron_vertices + num_polyhedron_edges + num_polyhedron_faces)
                continue;

            bool point_to_check = (local_mesh.Cell0DNumberNeighbourCell2D(c_0D) > 0);
            int cell3D_mark = 2;
            for (unsigned int n_f = 0; n_f < local_mesh.Cell0DNumberNeighbourCell2D(c_0D); n_f++)
            {
                if (!local_mesh.Cell0DHasNeighbourCell2D(c_0D, n_f))
                    continue;

                const unsigned int c_2D = local_mesh.Cell0DNeighbourCell2D(c_0D, n_f);

                if (local_mesh.Cell2DMarker(c_2D) != cell0D_marker)
                    continue;

                for (unsigned int n_c = 0; n_c < local_mesh.Cell2DNumberNeighbourCell3D(c_2D); n_c++)
                {
                    if (!local_mesh.Cell2DHasNeighbourCell3D(c_2D, n_c))
                        continue;

                    const unsigned int c_3D = local_mesh.Cell2DNeighbourCell3D(c_2D, n_c);

                    if (cell3D_mark == 2)
                        cell3D_mark = cell3Ds_mark[c_3D];

                    if (cell3Ds_mark[c_3D] != cell3D_mark)
                    {
                        point_to_check = false;
                        break;
                    }
                }

                if (!point_to_check)
                    break;
            }

            if (!point_to_check)
                continue;

            for (unsigned int n_c = 0; n_c < local_mesh.Cell0DNumberNeighbourCell3D(c_0D); n_c++)
            {
                if (!local_mesh.Cell0DHasNeighbourCell3D(c_0D, n_c))
                    continue;

                const unsigned int c_3D = local_mesh.Cell0DNeighbourCell3D(c_0D, n_c);

                if (cell3Ds_mark[c_3D] != cell3D_mark)
                {
                    std::cout << "point isolated - case not managed yet, return empty" << std::endl;
                    return {SplitPolyhedronByFunctionResult::Types::NoMarks, {}, {}, {}, {}, {}, {}, {}};
                }
            }
        }
    }

#ifdef EXPORT_TEST_VTU
    {
        {
            Gedim::VTKUtilities exporter;
            std::vector<double> marks(cell3Ds_mark.begin(), cell3Ds_mark.end());

            exporter.AddPolyhedrons(
                local_mesh.Cell0DsCoordinates(),
                local_mesh.Cell3DsFacesVertices(),
                {{"Mark", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(marks.size()), marks.data()}});
            exporter.Export(mesh_3D_export_folder + "/cell3Ds_mark_final.vtu");
        }
    }
#endif

    std::array<std::list<unsigned int>, 2> cell3Ds_mark_lists;
    for (unsigned int c = 0; c < local_mesh.Cell3DTotalNumber(); ++c)
    {
        if (cell3Ds_mark.at(c) == 0)
            cell3Ds_mark_lists[0].push_back(c);
        else
            cell3Ds_mark_lists[1].push_back(c);
    }

    std::vector<std::vector<unsigned int>> mesh_cell3Ds_original_cell3Ds(local_mesh.Cell3DTotalNumber());
    for (unsigned int c = 0; c < local_mesh.Cell3DTotalNumber(); ++c)
    {
        mesh_cell3Ds_original_cell3Ds.at(c) = std::vector<unsigned int>({c});
    }

    std::vector<std::vector<unsigned int>> mesh_cell2Ds_original_cell2Ds(local_mesh.Cell2DTotalNumber());
    for (unsigned int c = 0; c < local_mesh.Cell2DTotalNumber(); ++c)
    {
        mesh_cell2Ds_original_cell2Ds.at(c) = std::vector<unsigned int>({c});
    }

    std::vector<std::vector<unsigned int>> mesh_cell1Ds_original_cell1Ds(local_mesh.Cell1DTotalNumber());
    for (unsigned int c = 0; c < local_mesh.Cell1DTotalNumber(); ++c)
    {
        mesh_cell1Ds_original_cell1Ds.at(c) = std::vector<unsigned int>({c});
    }

    for (const auto &cell3Ds_mark_list : cell3Ds_mark_lists)
    {
        const std::unordered_set<unsigned int> cell3Ds(cell3Ds_mark_list.rbegin(), cell3Ds_mark_list.rend());

        const auto agglomerationInfo = Gedim::MeshUtilities::AgglomerateCell3Ds(cell3Ds, local_mesh);

        Gedim::MeshUtilities::AgglomerateCell3Ds(cell3Ds,
                                                 agglomerationInfo.AgglomerateCell3DVertices,
                                                 agglomerationInfo.AgglomerateCell3DEdges,
                                                 agglomerationInfo.AgglomerateCell3DFaces,
                                                 agglomerationInfo.SubCell3DsRemovedVertices,
                                                 agglomerationInfo.SubCell3DsRemovedEdges,
                                                 agglomerationInfo.SubCell3DsRemovedFaces,
                                                 local_mesh,
                                                 mesh_cell3Ds_original_cell3Ds);
    }

    {
        Gedim::MeshMatrices temp_mesh_data = local_mesh_data;
        Gedim::MeshMatricesDAO temp_mesh(temp_mesh_data);
        local_mesh_data = Gedim::MeshMatrices();

        const auto filter_active_mesh = Gedim::MeshUtilities::FilterActiveMesh(temp_mesh);
        const auto extract_mesh_data = Gedim::MeshUtilities::ExtractMesh3D(filter_active_mesh.Cell0Ds,
                                                                           filter_active_mesh.Cell1Ds,
                                                                           filter_active_mesh.Cell2Ds,
                                                                           filter_active_mesh.Cell3Ds,
                                                                           temp_mesh,
                                                                           local_mesh);

        std::vector<std::vector<unsigned int>> active_mesh_cell1Ds_original_cell1Ds(local_mesh.Cell1DTotalNumber());
        for (unsigned int c1D_index = 0; c1D_index < local_mesh.Cell1DTotalNumber(); c1D_index++)
        {
            const unsigned int active_c1D_temp_mesh = extract_mesh_data.NewCell1DToOldCell1D.at(c1D_index);
            active_mesh_cell1Ds_original_cell1Ds[c1D_index] = mesh_cell1Ds_original_cell1Ds.at(active_c1D_temp_mesh);
        }
        mesh_cell1Ds_original_cell1Ds = active_mesh_cell1Ds_original_cell1Ds;

        std::vector<std::vector<unsigned int>> active_mesh_cell2Ds_original_cell2Ds(local_mesh.Cell2DTotalNumber());
        for (unsigned int c2D_index = 0; c2D_index < local_mesh.Cell2DTotalNumber(); c2D_index++)
        {
            const unsigned int active_c2D_temp_mesh = extract_mesh_data.NewCell2DToOldCell2D.at(c2D_index);
            active_mesh_cell2Ds_original_cell2Ds[c2D_index] = mesh_cell2Ds_original_cell2Ds.at(active_c2D_temp_mesh);
        }
        mesh_cell2Ds_original_cell2Ds = active_mesh_cell2Ds_original_cell2Ds;

        std::vector<std::vector<unsigned int>> active_mesh_cell3Ds_original_cell3Ds(local_mesh.Cell3DTotalNumber());
        for (unsigned int c3D_index = 0; c3D_index < local_mesh.Cell3DTotalNumber(); c3D_index++)
        {
            const unsigned int active_c3D_temp_mesh = extract_mesh_data.NewCell3DToOldCell3D.at(c3D_index);
            active_mesh_cell3Ds_original_cell3Ds[c3D_index] = mesh_cell3Ds_original_cell3Ds.at(active_c3D_temp_mesh);
        }
        mesh_cell3Ds_original_cell3Ds = active_mesh_cell3Ds_original_cell3Ds;
    }

#ifdef EXPORT_TEST_VTU
    {
        Gedim::External::MeshUtilities::ExportMeshToVTU(local_mesh, mesh_3D_export_folder, "AgglomeratedMesh_Cell3Ds");
    }
#endif

    {
        Gedim::MeshUtilities::ComputeCell0DCell3DNeighbours(local_mesh);
        Gedim::MeshUtilities::ComputeCell0DCell1DNeighbours(local_mesh);
        Gedim::MeshUtilities::ComputeCell1DCell2DNeighbours(local_mesh);
        Gedim::MeshUtilities::ComputeCell2DCell3DNeighbours(local_mesh);

        const unsigned int num_agglomerated_cell3Ds = local_mesh.Cell3DTotalNumber();
        for (unsigned int c_3D = 0; c_3D < num_agglomerated_cell3Ds; ++c_3D)
        {
            std::vector<std::list<unsigned int>> cell2Ds_mark_lists(num_polyhedron_faces);

            std::unordered_set<unsigned int> new_cell3D_vertices;
            std::unordered_set<unsigned int> new_cell3D_edges;
            std::list<unsigned int> new_cell3D_faces;

            for (unsigned int f = 0; f < local_mesh.Cell3DNumberFaces(c_3D); ++f)
            {
                const unsigned c_2D = local_mesh.Cell3DFace(c_3D, f);
                const unsigned int cell2D_marker = local_mesh.Cell2DMarker(c_2D);

                if (cell2D_marker == 0 || cell2D_marker > num_polyhedron_vertices + num_polyhedron_edges + num_polyhedron_faces)
                {
                    for (unsigned int v = 0; v < local_mesh.Cell2DNumberVertices(c_2D); v++)
                    {
                        const unsigned int c_0D = local_mesh.Cell2DVertex(c_2D, v);
                        if (new_cell3D_vertices.find(c_0D) == new_cell3D_vertices.end())
                            new_cell3D_vertices.insert(c_0D);
                    }

                    for (unsigned int e = 0; e < local_mesh.Cell2DNumberEdges(c_2D); e++)
                    {
                        const unsigned int c_1D = local_mesh.Cell2DEdge(c_2D, e);
                        if (new_cell3D_edges.find(c_1D) == new_cell3D_edges.end())
                            new_cell3D_edges.insert(c_1D);
                    }

                    new_cell3D_faces.push_back(c_2D);
                    continue;
                }

                const unsigned int cell2Ds_mark_position = cell2D_marker - num_polyhedron_vertices - num_polyhedron_edges - 1;
                cell2Ds_mark_lists[cell2Ds_mark_position].push_back(c_2D);
            }

            std::list<std::unordered_set<unsigned int>> cell2Ds_to_agglomerate;
            {
                for (const auto &cell2Ds : cell2Ds_mark_lists)
                {
                    std::vector<unsigned int> cell2Ds_index(cell2Ds.begin(), cell2Ds.end());
                    std::unordered_map<unsigned int, unsigned int> cell2Ds_map;
                    for (const unsigned int c_2D : cell2Ds)
                        cell2Ds_map.insert(std::make_pair(c_2D, cell2Ds_map.size()));
                    std::vector<std::list<unsigned int>> cell2Ds_adj_lst(cell2Ds.size());

                    for (const auto &c_2D : cell2Ds_map)
                    {
                        for (unsigned int e = 0; e < local_mesh.Cell2DNumberEdges(c_2D.first); e++)
                        {
                            const unsigned int c_1D = local_mesh.Cell2DEdge(c_2D.first, e);

                            for (unsigned int n = 0; n < local_mesh.Cell1DNumberNeighbourCell2D(c_1D); n++)
                            {
                                if (!local_mesh.Cell1DHasNeighbourCell2D(c_1D, n))
                                    continue;

                                const unsigned int c_2D_n = local_mesh.Cell1DNeighbourCell2D(c_1D, n);

                                if (c_2D_n == c_2D.first)
                                    continue;

                                const auto c_2D_n_pos = cell2Ds_map.find(c_2D_n);
                                if (c_2D_n_pos == cell2Ds_map.end())
                                    continue;

                                cell2Ds_adj_lst[c_2D.second].push_back(c_2D_n_pos->second);
                            }
                        }
                    }

                    std::vector<std::vector<unsigned int>> cell2Ds_adj(cell2Ds_adj_lst.size());
                    for (unsigned int c = 0; c < cell2Ds_adj_lst.size(); c++)
                    {
                        cell2Ds_adj[c] = std::vector<unsigned int>(cell2Ds_adj_lst.at(c).begin(), cell2Ds_adj_lst.at(c).end());
                    }

                    const auto connected_components = graphUtilities.ComputeStronglyConnectedComponents(cell2Ds_adj);

                    for (const auto &connected_component : connected_components)
                    {
                        cell2Ds_to_agglomerate.push_back(std::unordered_set<unsigned int>());
                        auto &cell2Ds_to_fill = cell2Ds_to_agglomerate.back();

                        for (const auto &c_2D_pos : connected_component)
                            cell2Ds_to_fill.insert(cell2Ds_index.at(c_2D_pos));
                    }
                }
            }

            for (const auto &cell2Ds : cell2Ds_to_agglomerate)
            {
                if (cell2Ds.empty())
                    continue;

                if (cell2Ds.size() == 1)
                {
                    const unsigned int c_2D = *cell2Ds.begin();

                    for (unsigned int v = 0; v < local_mesh.Cell2DNumberVertices(c_2D); v++)
                    {
                        const unsigned int c_0D = local_mesh.Cell2DVertex(c_2D, v);
                        if (new_cell3D_vertices.find(c_0D) == new_cell3D_vertices.end())
                            new_cell3D_vertices.insert(c_0D);
                    }

                    for (unsigned int e = 0; e < local_mesh.Cell2DNumberEdges(c_2D); e++)
                    {
                        const unsigned int c_1D = local_mesh.Cell2DEdge(c_2D, e);
                        if (new_cell3D_edges.find(c_1D) == new_cell3D_edges.end())
                            new_cell3D_edges.insert(c_1D);
                    }

                    new_cell3D_faces.push_back(c_2D);
                    continue;
                }

                {
                    std::unordered_set<unsigned int> aggl_cell0Ds;
                    std::unordered_set<unsigned int> aggl_cell1Ds;

                    for (const auto &c_2D : cell2Ds)
                    {
                        for (unsigned int v = 0; v < local_mesh.Cell2DNumberVertices(c_2D); v++)
                        {
                            const unsigned int c_0D = local_mesh.Cell2DVertex(c_2D, v);

                            if (aggl_cell0Ds.find(c_0D) == aggl_cell0Ds.end())
                                aggl_cell0Ds.insert(c_0D);
                        }

                        for (unsigned int e = 0; e < local_mesh.Cell2DNumberEdges(c_2D); e++)
                        {
                            const unsigned int c_1D = local_mesh.Cell2DEdge(c_2D, e);

                            if (aggl_cell1Ds.find(c_1D) == aggl_cell1Ds.end())
                                aggl_cell1Ds.insert(c_1D);
                        }
                    }

                    const unsigned int num_cell0Ds_aggl = aggl_cell0Ds.size();
                    const unsigned int num_cell1Ds_aggl = aggl_cell1Ds.size();
                    const unsigned int num_cell2Ds_aggl = cell2Ds.size();

                    const int euler_char = num_cell0Ds_aggl - num_cell1Ds_aggl + num_cell2Ds_aggl;

                    if (euler_char != 1)
                    {
                        std::cout << "Euler characteristic face: " << euler_char << " - case not managed yet" << std::endl;
                        return {SplitPolyhedronByFunctionResult::Types::NoMarks, {}, {}, {}, {}, {}, {}, {}};
                    }
                }

                const auto agglomerationInfo = Gedim::MeshUtilities::AgglomerateCell2Ds(geometryUtilities, cell2Ds, local_mesh);

                const unsigned int agglomerated_cell2D_index =
                    Gedim::MeshUtilities::AgglomerateCell2Ds(cell2Ds,
                                                             agglomerationInfo.AgglomerateCell2DVertices,
                                                             agglomerationInfo.AgglomerateCell2DEdges,
                                                             agglomerationInfo.SubCell2DsRemovedVertices,
                                                             agglomerationInfo.SubCell2DsRemovedEdges,
                                                             local_mesh,
                                                             mesh_cell2Ds_original_cell2Ds);

                for (unsigned int v = 0; v < local_mesh.Cell2DNumberVertices(agglomerated_cell2D_index); v++)
                {
                    const unsigned int c_0D = local_mesh.Cell2DVertex(agglomerated_cell2D_index, v);
                    if (new_cell3D_vertices.find(c_0D) == new_cell3D_vertices.end())
                        new_cell3D_vertices.insert(c_0D);
                }

                for (unsigned int e = 0; e < local_mesh.Cell2DNumberEdges(agglomerated_cell2D_index); e++)
                {
                    const unsigned int c_1D = local_mesh.Cell2DEdge(agglomerated_cell2D_index, e);
                    if (new_cell3D_edges.find(c_1D) == new_cell3D_edges.end())
                        new_cell3D_edges.insert(c_1D);
                }

                new_cell3D_faces.push_back(agglomerated_cell2D_index);
            }

            Gedim::MeshUtilities::SplitCell3D(
                c_3D,
                {std::vector<unsigned int>(new_cell3D_vertices.begin(), new_cell3D_vertices.end())},
                {std::vector<unsigned int>(new_cell3D_edges.begin(), new_cell3D_edges.end())},
                {std::vector<unsigned int>(new_cell3D_faces.begin(), new_cell3D_faces.end())},
                local_mesh);

            mesh_cell3Ds_original_cell3Ds.push_back(mesh_cell3Ds_original_cell3Ds.at(c_3D));
        }
        mesh_cell3Ds_original_cell3Ds.shrink_to_fit();
    }

    {
        Gedim::MeshMatrices temp_mesh_data = local_mesh_data;
        Gedim::MeshMatricesDAO temp_mesh(temp_mesh_data);
        local_mesh_data = Gedim::MeshMatrices();

        const auto filter_active_mesh = Gedim::MeshUtilities::FilterActiveMesh(temp_mesh);
        const auto extract_mesh_data = Gedim::MeshUtilities::ExtractMesh3D(filter_active_mesh.Cell0Ds,
                                                                           filter_active_mesh.Cell1Ds,
                                                                           filter_active_mesh.Cell2Ds,
                                                                           filter_active_mesh.Cell3Ds,
                                                                           temp_mesh,
                                                                           local_mesh);

        std::vector<std::vector<unsigned int>> active_mesh_cell1Ds_original_cell1Ds(local_mesh.Cell1DTotalNumber());
        for (unsigned int c1D_index = 0; c1D_index < local_mesh.Cell1DTotalNumber(); c1D_index++)
        {
            const unsigned int active_c1D_temp_mesh = extract_mesh_data.NewCell1DToOldCell1D.at(c1D_index);
            active_mesh_cell1Ds_original_cell1Ds[c1D_index] = mesh_cell1Ds_original_cell1Ds.at(active_c1D_temp_mesh);
        }
        mesh_cell1Ds_original_cell1Ds = active_mesh_cell1Ds_original_cell1Ds;

        std::vector<std::vector<unsigned int>> active_mesh_cell2Ds_original_cell2Ds(local_mesh.Cell2DTotalNumber());
        for (unsigned int c2D_index = 0; c2D_index < local_mesh.Cell2DTotalNumber(); c2D_index++)
        {
            const unsigned int active_c2D_temp_mesh = extract_mesh_data.NewCell2DToOldCell2D.at(c2D_index);
            active_mesh_cell2Ds_original_cell2Ds[c2D_index] = mesh_cell2Ds_original_cell2Ds.at(active_c2D_temp_mesh);
        }
        mesh_cell2Ds_original_cell2Ds = active_mesh_cell2Ds_original_cell2Ds;

        std::vector<std::vector<unsigned int>> active_mesh_cell3Ds_original_cell3Ds(local_mesh.Cell3DTotalNumber());
        for (unsigned int c3D_index = 0; c3D_index < local_mesh.Cell3DTotalNumber(); c3D_index++)
        {
            const unsigned int active_c3D_temp_mesh = extract_mesh_data.NewCell3DToOldCell3D.at(c3D_index);
            active_mesh_cell3Ds_original_cell3Ds[c3D_index] = mesh_cell3Ds_original_cell3Ds.at(active_c3D_temp_mesh);
        }
        mesh_cell3Ds_original_cell3Ds = active_mesh_cell3Ds_original_cell3Ds;
    }

#ifdef EXPORT_TEST_VTU
    {
        Gedim::External::MeshUtilities::ExportMeshToVTU(local_mesh, mesh_3D_export_folder, "AgglomeratedMesh_Cell2Ds");
    }
#endif

    {
        Gedim::MeshUtilities::ComputeCell0DCell3DNeighbours(local_mesh);
        Gedim::MeshUtilities::ComputeCell0DCell1DNeighbours(local_mesh);
        Gedim::MeshUtilities::ComputeCell1DCell2DNeighbours(local_mesh);
        Gedim::MeshUtilities::ComputeCell2DCell3DNeighbours(local_mesh);

        const unsigned int num_agglomerated_cell2Ds = local_mesh.Cell2DTotalNumber();
        std::vector<unsigned int> cell1Ds_agglomerated_cell1D(local_mesh.Cell1DTotalNumber(), 0);
        std::vector<unsigned int> cell2Ds_agglomerated_cell2D(local_mesh.Cell2DTotalNumber(), 0);

        for (unsigned int c_2D = 0; c_2D < num_agglomerated_cell2Ds; ++c_2D)
        {
            std::list<unsigned int> new_cell2D_vertices_list;
            std::list<unsigned int> new_cell2D_edges_list;

            bool to_update = false;
            unsigned int last_marker = 0;
            unsigned int last_agglomerated_cell1D = 0;
            std::list<std::pair<std::unordered_set<unsigned int>, unsigned int>> cell1Ds_to_agglomerate_lists;

            const unsigned int c_2D_num_edges = local_mesh.Cell2DNumberEdges(c_2D);
            unsigned int starting_e = 0;
            while (starting_e < c_2D_num_edges && local_mesh.Cell0DMarker(local_mesh.Cell2DVertex(c_2D, starting_e)) > num_polyhedron_vertices)
            {
                starting_e++;
            }

            if (starting_e == c_2D_num_edges)
                continue;

            for (unsigned int f_e = 0; f_e < c_2D_num_edges; ++f_e)
            {
                const unsigned int edge_position = (starting_e + f_e) % c_2D_num_edges;

                const unsigned c_1D = local_mesh.Cell2DEdge(c_2D, edge_position);
                const unsigned int cell1D_marker = local_mesh.Cell1DMarker(c_1D);

                if (cell1D_marker == 0 || cell1D_marker > num_polyhedron_vertices + num_polyhedron_edges)
                {
                    new_cell2D_vertices_list.push_back(local_mesh.Cell2DVertex(c_2D, edge_position));
                    new_cell2D_edges_list.push_back(c_1D);
                    last_agglomerated_cell1D = 0;
                    last_marker = 0;
                    continue;
                }

                const unsigned int agglomerated_cell1D = cell1Ds_agglomerated_cell1D.at(c_1D);

                if (agglomerated_cell1D > 0 && agglomerated_cell1D == last_agglomerated_cell1D)
                {
                    continue;
                }

                if (agglomerated_cell1D > 0 && agglomerated_cell1D != last_agglomerated_cell1D)
                {
                    to_update = true;
                    last_agglomerated_cell1D = agglomerated_cell1D;
                    last_marker = local_mesh.Cell1DMarker(agglomerated_cell1D);
                    new_cell2D_vertices_list.push_back(local_mesh.Cell2DVertex(c_2D, edge_position));
                    new_cell2D_edges_list.push_back(agglomerated_cell1D);
                    continue;
                }

                const unsigned int origin_vertex_index = local_mesh.Cell2DVertex(c_2D, edge_position);

                if (cell1Ds_to_agglomerate_lists.empty() || (local_mesh.Cell0DNumberNeighbourCell1D(origin_vertex_index) > 2) ||
                    (last_marker > 0 && cell1D_marker != last_marker))
                {
                    to_update = true;
                    last_marker = cell1D_marker;
                    cell1Ds_to_agglomerate_lists.push_back(
                        std::make_pair(std::unordered_set<unsigned int>(), new_cell2D_edges_list.size()));
                    new_cell2D_vertices_list.push_back(local_mesh.Cell2DVertex(c_2D, edge_position));
                    new_cell2D_edges_list.push_back(c_1D);
                }

                auto &cell1Ds_to_agglomerate_list = cell1Ds_to_agglomerate_lists.back();
                cell1Ds_to_agglomerate_list.first.insert(c_1D);
            }

            if (!to_update)
                continue;

            std::vector<unsigned int> new_cell2D_vertices(new_cell2D_vertices_list.begin(), new_cell2D_vertices_list.end());
            std::vector<unsigned int> new_cell2D_edges(new_cell2D_edges_list.begin(), new_cell2D_edges_list.end());

            for (const auto &cell1Ds_to_agglomerate_list : cell1Ds_to_agglomerate_lists)
            {
                if (cell1Ds_to_agglomerate_list.first.empty())
                    continue;

                if (cell1Ds_to_agglomerate_list.first.size() == 1)
                    continue;

                const auto agglomerationInfo =
                    Gedim::MeshUtilities::AgglomerateCell1Ds(cell1Ds_to_agglomerate_list.first, local_mesh);

                const unsigned int agglomerated_cell1D_index =
                    Gedim::MeshUtilities::AgglomerateCell1Ds(cell1Ds_to_agglomerate_list.first,
                                                             agglomerationInfo.AgglomerateCell1DVertices,
                                                             agglomerationInfo.SubCell1DsRemovedVertices,
                                                             local_mesh,
                                                             mesh_cell1Ds_original_cell1Ds);

                for (const auto &c_1D : cell1Ds_to_agglomerate_list.first)
                    cell1Ds_agglomerated_cell1D[c_1D] = agglomerated_cell1D_index;

                new_cell2D_edges[cell1Ds_to_agglomerate_list.second] = agglomerated_cell1D_index;
            }

            Eigen::MatrixXi new_cell2d(2, new_cell2D_vertices.size());
            for (unsigned int v = 0; v < new_cell2D_vertices.size(); v++)
                new_cell2d.col(v) << new_cell2D_vertices.at(v), new_cell2D_edges.at(v);

            const auto new_cell2Ds = Gedim::MeshUtilities::SplitCell2D(c_2D, {new_cell2d}, local_mesh);

            mesh_cell2Ds_original_cell2Ds.push_back(mesh_cell2Ds_original_cell2Ds.at(c_2D));
            cell2Ds_agglomerated_cell2D[c_2D] = new_cell2Ds.at(0);
        }
        mesh_cell2Ds_original_cell2Ds.shrink_to_fit();

        const unsigned int num_agglomerated_cell3Ds = local_mesh.Cell3DTotalNumber();
        for (unsigned int c_3D = 0; c_3D < num_agglomerated_cell3Ds; ++c_3D)
        {
            if (!local_mesh.Cell3DIsActive(c_3D))
                continue;

            std::unordered_set<unsigned int> cell3D_new_cell0Ds;
            std::unordered_set<unsigned int> cell3D_new_cell1Ds;
            std::unordered_set<unsigned int> cell3D_new_cell2Ds;

            for (unsigned int f = 0; f < local_mesh.Cell3DNumberFaces(c_3D); f++)
            {
                const unsigned int original_c2D = local_mesh.Cell3DFace(c_3D, f);
                unsigned int aggl_c2D = cell2Ds_agglomerated_cell2D.at(original_c2D);

                if (aggl_c2D == 0)
                    aggl_c2D = original_c2D;

                if (cell3D_new_cell2Ds.find(aggl_c2D) != cell3D_new_cell2Ds.end())
                    continue;

                cell3D_new_cell2Ds.insert(aggl_c2D);

                for (unsigned int e = 0; e < local_mesh.Cell2DNumberEdges(aggl_c2D); e++)
                {
                    const unsigned int aggl_c1D = local_mesh.Cell2DEdge(aggl_c2D, e);

                    if (cell3D_new_cell1Ds.find(aggl_c1D) != cell3D_new_cell1Ds.end())
                        continue;

                    cell3D_new_cell1Ds.insert(aggl_c1D);
                }

                for (unsigned int v = 0; v < local_mesh.Cell2DNumberVertices(aggl_c2D); v++)
                {
                    const unsigned int aggl_c0D = local_mesh.Cell2DVertex(aggl_c2D, v);

                    if (cell3D_new_cell0Ds.find(aggl_c0D) != cell3D_new_cell0Ds.end())
                        continue;

                    cell3D_new_cell0Ds.insert(aggl_c0D);
                }
            }

            Gedim::MeshUtilities::SplitCell3D(
                c_3D,
                {std::vector<unsigned int>(cell3D_new_cell0Ds.begin(), cell3D_new_cell0Ds.end())},
                {std::vector<unsigned int>(cell3D_new_cell1Ds.begin(), cell3D_new_cell1Ds.end())},
                {std::vector<unsigned int>(cell3D_new_cell2Ds.begin(), cell3D_new_cell2Ds.end())},
                local_mesh);

            mesh_cell3Ds_original_cell3Ds.push_back(mesh_cell3Ds_original_cell3Ds.at(c_3D));
        }

        mesh_cell3Ds_original_cell3Ds.shrink_to_fit();
    }

    {
        Gedim::MeshMatrices temp_mesh_data = local_mesh_data;
        Gedim::MeshMatricesDAO temp_mesh(temp_mesh_data);
        local_mesh_data = Gedim::MeshMatrices();

        const auto filter_active_mesh = Gedim::MeshUtilities::FilterActiveMesh(temp_mesh);
        const auto extract_mesh_data = Gedim::MeshUtilities::ExtractMesh3D(filter_active_mesh.Cell0Ds,
                                                                           filter_active_mesh.Cell1Ds,
                                                                           filter_active_mesh.Cell2Ds,
                                                                           filter_active_mesh.Cell3Ds,
                                                                           temp_mesh,
                                                                           local_mesh);

        std::vector<std::vector<unsigned int>> active_mesh_cell1Ds_original_cell1Ds(local_mesh.Cell1DTotalNumber());
        for (unsigned int c1D_index = 0; c1D_index < local_mesh.Cell1DTotalNumber(); c1D_index++)
        {
            const unsigned int active_c1D_temp_mesh = extract_mesh_data.NewCell1DToOldCell1D.at(c1D_index);
            active_mesh_cell1Ds_original_cell1Ds[c1D_index] = mesh_cell1Ds_original_cell1Ds.at(active_c1D_temp_mesh);
        }
        mesh_cell1Ds_original_cell1Ds = active_mesh_cell1Ds_original_cell1Ds;

        std::vector<std::vector<unsigned int>> active_mesh_cell2Ds_original_cell2Ds(local_mesh.Cell2DTotalNumber());
        for (unsigned int c2D_index = 0; c2D_index < local_mesh.Cell2DTotalNumber(); c2D_index++)
        {
            const unsigned int active_c2D_temp_mesh = extract_mesh_data.NewCell2DToOldCell2D.at(c2D_index);
            active_mesh_cell2Ds_original_cell2Ds[c2D_index] = mesh_cell2Ds_original_cell2Ds.at(active_c2D_temp_mesh);
        }
        mesh_cell2Ds_original_cell2Ds = active_mesh_cell2Ds_original_cell2Ds;

        std::vector<std::vector<unsigned int>> active_mesh_cell3Ds_original_cell3Ds(local_mesh.Cell3DTotalNumber());
        for (unsigned int c3D_index = 0; c3D_index < local_mesh.Cell3DTotalNumber(); c3D_index++)
        {
            const unsigned int active_c3D_temp_mesh = extract_mesh_data.NewCell3DToOldCell3D.at(c3D_index);
            active_mesh_cell3Ds_original_cell3Ds[c3D_index] = mesh_cell3Ds_original_cell3Ds.at(active_c3D_temp_mesh);
        }
        mesh_cell3Ds_original_cell3Ds = active_mesh_cell3Ds_original_cell3Ds;
    }

#ifdef EXPORT_TEST_VTU
    {
        {
            Gedim::External::MeshUtilities::ExportMeshToVTU(local_mesh, mesh_3D_export_folder, "AgglomeratedMesh_Cell1Ds");
        }

        {
            std::unordered_set<unsigned int> original_cell3Ds_index;
            Gedim::VTKUtilities exporter;

            double total_volume = 0.0;
            for (unsigned int c_3D = 0; c_3D < local_mesh.Cell3DTotalNumber(); c_3D++)
            {
                if (!local_mesh.Cell3DIsActive(c_3D))
                    continue;

                for (const auto &original_cell3D_index : mesh_cell3Ds_original_cell3Ds.at(c_3D))
                {
                    if (original_cell3Ds_index.find(original_cell3D_index) != original_cell3Ds_index.end())
                        throw std::runtime_error("Double tetra found!");

                    original_cell3Ds_index.insert(original_cell3D_index);

                    std::vector<double> final_cell3D(1, c_3D);

                    const auto original_polyhedron =
                        Gedim::MeshUtilities::MeshCell3DToPolyhedron(original_mesh, original_cell3D_index);

                    total_volume += geometryUtilities.PolyhedronVolumeByInternalIntegral({original_polyhedron.Vertices});

                    exporter.AddPolyhedron(original_polyhedron.Vertices,
                                           original_polyhedron.Edges,
                                           original_polyhedron.Faces,
                                           {{"Final_Cell3D",
                                             Gedim::VTPProperty::Formats::Cells,
                                             static_cast<unsigned int>(final_cell3D.size()),
                                             final_cell3D.data()}});
                }
            }

            exporter.Export(mesh_3D_export_folder + "/final_original_cell3Ds.vtu");

            Gedim::Output::Assert(geometryUtilities.AreValuesEqual(polyhedron_volume, total_volume, geometryUtilities.Tolerance3D()));
        }

        {
            std::unordered_set<unsigned int> original_cell2Ds_index;
            Gedim::VTKUtilities exporter;

            for (unsigned int c_2D = 0; c_2D < local_mesh.Cell2DTotalNumber(); c_2D++)
            {
                if (!local_mesh.Cell2DIsActive(c_2D))
                    continue;

                for (const auto &original_cell2D_index : mesh_cell2Ds_original_cell2Ds.at(c_2D))
                {
                    if (original_cell2Ds_index.find(original_cell2D_index) != original_cell2Ds_index.end())
                        throw std::runtime_error("Double triangle found!");

                    original_cell2Ds_index.insert(original_cell2D_index);

                    std::vector<double> final_cell2D(1, c_2D);

                    exporter.AddPolygon(original_mesh.Cell2DVerticesCoordinates(original_cell2D_index),
                                        {{"Final_Cell2D",
                                          Gedim::VTPProperty::Formats::Cells,
                                          static_cast<unsigned int>(final_cell2D.size()),
                                          final_cell2D.data()}});
                }
            }

            exporter.Export(mesh_3D_export_folder + "/final_original_cell2Ds.vtu");
        }

        {
            std::unordered_set<unsigned int> original_cell1Ds_index;
            Gedim::VTKUtilities exporter;

            for (unsigned int c_1D = 0; c_1D < local_mesh.Cell1DTotalNumber(); c_1D++)
            {
                if (!local_mesh.Cell1DIsActive(c_1D))
                    continue;

                for (const auto &original_cell1D_index : mesh_cell1Ds_original_cell1Ds.at(c_1D))
                {
                    if (original_cell1Ds_index.find(original_cell1D_index) != original_cell1Ds_index.end())
                        throw std::runtime_error("Double segment found!");

                    original_cell1Ds_index.insert(original_cell1D_index);

                    std::vector<double> final_cell1D(1, c_1D);

                    exporter.AddSegment(original_mesh.Cell1DOriginCoordinates(original_cell1D_index),
                                        original_mesh.Cell1DEndCoordinates(original_cell1D_index),
                                        {{"Final_Cell1D",
                                          Gedim::VTPProperty::Formats::Cells,
                                          static_cast<unsigned int>(final_cell1D.size()),
                                          final_cell1D.data()}});
                }
            }

            exporter.Export(mesh_3D_export_folder + "/final_original_cell1Ds.vtu");
        }
    }
#endif

    SplitPolyhedronByFunctionResult result;

    result.Type = SplitPolyhedronByFunctionResult::Types::Splitted;
    result.Original_vertex_to_split_vertex.resize(num_polyhedron_vertices);
    result.Split_vertices.resize(local_mesh.Cell0DTotalNumber());

    for (unsigned int v = 0; v < local_mesh.Cell0DTotalNumber(); v++)
    {
        const unsigned int vertex_marker = local_mesh.Cell0DMarker(v);

        if (vertex_marker == 0 || vertex_marker > num_polyhedron_vertices)
        {
            result.Split_vertices[v] = {local_mesh.Cell0DCoordinates(v), 0, vertex_marker, true};

            continue;
        }

        const unsigned int vertex_position = vertex_marker - 1;
        result.Original_vertex_to_split_vertex[vertex_position] = v;
        result.Split_vertices[v] = {local_mesh.Cell0DCoordinates(v), vertex_position, vertex_marker, false};
    }

    std::vector<std::list<unsigned int>> original_edge_to_split_edges(num_polyhedron_edges);
    result.Split_edges.resize(local_mesh.Cell1DTotalNumber());

    for (unsigned int e = 0; e < local_mesh.Cell1DTotalNumber(); e++)
    {
        const unsigned int edge_marker = local_mesh.Cell1DMarker(e);
        const std::array<unsigned int, 2> edge_extreme = {local_mesh.Cell1DOrigin(e), local_mesh.Cell1DEnd(e)};

        if (edge_marker == 0 || edge_marker > num_polyhedron_vertices + num_polyhedron_edges)
        {
            result.Split_edges[e] = {edge_extreme, 0, edge_marker, true};
            continue;
        }

        const unsigned int edge_position = edge_marker - num_polyhedron_vertices - 1;
        original_edge_to_split_edges[edge_position].push_back(e);
        result.Split_edges[e] = {edge_extreme, edge_position, edge_marker, false};
    }

    result.Original_edge_to_split_edges.resize(num_polyhedron_edges);
    for (unsigned int e = 0; e < num_polyhedron_edges; e++)
    {
        result.Original_edge_to_split_edges[e] =
            std::vector<unsigned int>(original_edge_to_split_edges.at(e).begin(), original_edge_to_split_edges.at(e).end());
    }

    std::vector<std::list<unsigned int>> original_face_to_split_faces(num_polyhedron_faces);
    result.Split_faces.resize(local_mesh.Cell2DTotalNumber());

    for (unsigned int f = 0; f < local_mesh.Cell2DTotalNumber(); f++)
    {
        const unsigned int face_marker = local_mesh.Cell2DMarker(f);

        Eigen::MatrixXi face_extreme(2, local_mesh.Cell2DNumberVertices(f));
        for (unsigned int v = 0; v < local_mesh.Cell2DNumberVertices(f); v++)
            face_extreme.col(v) << local_mesh.Cell2DVertex(f, v), local_mesh.Cell2DEdge(f, v);

        std::list<Eigen::MatrixXd> triangles;
        for (const auto &original_c2D_index : mesh_cell2Ds_original_cell2Ds.at(f))
            triangles.push_back(original_mesh.Cell2DVerticesCoordinates(original_c2D_index));

        if (face_marker == 0 || face_marker > num_polyhedron_vertices + num_polyhedron_edges + num_polyhedron_faces)
        {
            result.Split_faces[f] = {face_extreme,
                                     0,
                                     face_marker,
                                     true,
                                     std::vector<Eigen::Matrix3d>(triangles.begin(), triangles.end())};
            continue;
        }

        const unsigned int face_position = face_marker - num_polyhedron_vertices - num_polyhedron_edges - 1;
        original_face_to_split_faces[face_position].push_back(f);
        result.Split_faces[f] = {face_extreme,
                                 face_position,
                                 face_marker,
                                 false,
                                 std::vector<Eigen::Matrix3d>(triangles.begin(), triangles.end())};
    }

    result.Original_face_to_split_faces.resize(num_polyhedron_faces);
    for (unsigned int e = 0; e < num_polyhedron_faces; e++)
    {
        result.Original_face_to_split_faces[e] =
            std::vector<unsigned int>(original_face_to_split_faces.at(e).begin(), original_face_to_split_faces.at(e).end());
    }

    result.Split_polyhedrons.resize(local_mesh.Cell3DTotalNumber());

    for (unsigned int p = 0; p < local_mesh.Cell3DTotalNumber(); p++)
    {
        std::list<Eigen::MatrixXd> tetrahedrons;
        unsigned int c3D_marker;

        for (const auto &original_c3D_index : mesh_cell3Ds_original_cell3Ds.at(p))
        {
            c3D_marker = cell3Ds_mark.at(original_c3D_index);
            const auto tetra_vertices = original_mesh.Cell3DVerticesCoordinates(original_c3D_index);

            if (geometryUtilities.IsValueZero(geometryUtilities.SimplexMeasure(tetra_vertices), geometryUtilities.Tolerance3D()))
                continue;

            tetrahedrons.push_back(tetra_vertices);
        }

        result.Split_polyhedrons[p] = {local_mesh.Cell3DVertices(p),
                                       local_mesh.Cell3DEdges(p),
                                       local_mesh.Cell3DFaces(p),
                                       std::vector<Eigen::MatrixXd>(tetrahedrons.begin(), tetrahedrons.end()),
                                       c3D_marker};
    }

    return result;
}
// ***************************************************************************
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
                                   const double &resolution)
{
    const auto &polyhedron_vertices = cell3D_Vertices;
    const auto &polyhedron_edges = cell3D_Edges;
    const auto &polyhedron_faces = cell3D_Faces;
    const double &cell3D_volume = cell3Ds_volume.at(cell3D_index);
    const auto &cell3D_additional_constrained_points = cell3Ds_additional_constrained_points.at(cell3D_index);

    const unsigned int num_polyhedron_vertices = polyhedron_vertices.cols();
    const unsigned int num_polyhedron_edges = polyhedron_edges.cols();
    const unsigned int num_polyhedron_faces = polyhedron_faces.size();

    const auto split_result = SplitPolyhedronByFunction(geometryUtilities,
                                                        graphUtilities,
                                                        polyhedron_vertices,
                                                        polyhedron_edges,
                                                        polyhedron_faces,
                                                        cell3D_additional_constrained_points,
                                                        cell3D_volume,
                                                        marking_function,
                                                        resolution,
                                                        cell3D_index);

    if (split_result.Split_polyhedrons.size() == 0)
    {
        switch (split_result.Type)
        {
        case SplitPolyhedronByFunctionResult::Types::NoMarks:
            cell3Ds_mark[cell3D_index] = 0;
            break;
        case SplitPolyhedronByFunctionResult::Types::AllMarked:
            cell3Ds_mark[cell3D_index] = 1;
            break;
        default:
            throw std::runtime_error("Not managed type");
        }

        return false;
    }

    if (split_result.Split_polyhedrons.size() == 1)
    {
        switch (split_result.Type)
        {
        case SplitPolyhedronByFunctionResult::Types::NoMarks:
            cell3Ds_mark[cell3D_index] = 0;
            break;
        case SplitPolyhedronByFunctionResult::Types::AllMarked:
            cell3Ds_mark[cell3D_index] = 1;
            break;
        default:
            throw std::runtime_error("Not managed type");
        }

        return false;
    }

    std::vector<unsigned int> mesh_cell0Ds_index(split_result.Split_vertices.size());
    for (unsigned int v = 0; v < split_result.Split_vertices.size(); v++)
    {
        const auto &vertex = split_result.Split_vertices.at(v);

        if (!vertex.IsNew)
        {
            mesh_cell0Ds_index[v] = mesh.Cell3DVertex(cell3D_index, vertex.Original_vertex_index);
            continue;
        }

        mesh_cell0Ds_index[v] = mesh.Cell0DAppend(1);
        const unsigned int &newCell0DIndex = mesh_cell0Ds_index[v];
        mesh.Cell0DInsertCoordinates(newCell0DIndex, vertex.Coordinates);
        mesh.Cell0DSetState(newCell0DIndex, true);

        if (vertex.Marker == 0 || vertex.Marker > num_polyhedron_vertices + num_polyhedron_edges + num_polyhedron_faces)
        {
            mesh.Cell0DSetMarker(newCell0DIndex, mesh.Cell3DMarker(cell3D_index));
        }
        else if (vertex.Marker < num_polyhedron_vertices + 1)
        {
            const unsigned int vertex_pos = vertex.Marker - 1;
            mesh.Cell0DSetMarker(newCell0DIndex, mesh.Cell0DMarker(mesh.Cell3DVertex(cell3D_index, vertex_pos)));
        }
        else if (vertex.Marker < num_polyhedron_vertices + num_polyhedron_edges + 1)
        {
            const unsigned int edge_pos = vertex.Marker - num_polyhedron_vertices - 1;
            mesh.Cell0DSetMarker(newCell0DIndex, mesh.Cell1DMarker(mesh.Cell3DEdge(cell3D_index, edge_pos)));
        }
        else if (vertex.Marker < num_polyhedron_vertices + num_polyhedron_edges + num_polyhedron_faces + 1)
        {
            const unsigned int face_pos = vertex.Marker - num_polyhedron_vertices - num_polyhedron_edges - 1;
            mesh.Cell0DSetMarker(newCell0DIndex, mesh.Cell2DMarker(mesh.Cell3DFace(cell3D_index, face_pos)));
        }
    }

    std::vector<unsigned int> mesh_cell1Ds_index(split_result.Split_edges.size());
    for (unsigned int e = 0; e < split_result.Split_edges.size(); e++)
    {
        const auto &edge = split_result.Split_edges.at(e);

        if (!edge.IsNew)
            continue;

        mesh_cell1Ds_index[e] = mesh.Cell1DAppend(1);
        const unsigned int &newCell1DIndex = mesh_cell1Ds_index.at(e);
        mesh.Cell1DInsertExtremes(newCell1DIndex,
                                  mesh_cell0Ds_index.at(edge.Extremes.at(0)),
                                  mesh_cell0Ds_index.at(edge.Extremes.at(1)));
        mesh.Cell1DSetState(newCell1DIndex, true);

        if (edge.Marker < num_polyhedron_vertices + 1 || edge.Marker > num_polyhedron_vertices + num_polyhedron_edges + num_polyhedron_faces)
        {
            mesh.Cell1DSetMarker(newCell1DIndex, mesh.Cell3DMarker(cell3D_index));
        }
        else if (edge.Marker < num_polyhedron_vertices + num_polyhedron_edges + 1)
        {
            const unsigned int edge_pos = edge.Marker - num_polyhedron_vertices - 1;
            mesh.Cell1DSetMarker(newCell1DIndex, mesh.Cell1DMarker(mesh.Cell3DEdge(cell3D_index, edge_pos)));
        }
        else if (edge.Marker < num_polyhedron_vertices + num_polyhedron_edges + num_polyhedron_faces + 1)
        {
            const unsigned int face_pos = edge.Marker - num_polyhedron_vertices - num_polyhedron_edges - 1;
            mesh.Cell1DSetMarker(newCell1DIndex, mesh.Cell2DMarker(mesh.Cell3DFace(cell3D_index, face_pos)));
        }
    }

    for (unsigned int e = 0; e < split_result.Original_edge_to_split_edges.size(); e++)
    {
        const auto &split_edges = split_result.Original_edge_to_split_edges.at(e);

        const unsigned int original_cell1D = mesh.Cell3DEdge(cell3D_index, e);

        if (split_edges.size() == 1)
        {
            mesh_cell1Ds_index[split_edges.at(0)] = original_cell1D;
            continue;
        }

        const unsigned int new_edges_size = split_edges.size();
        Eigen::MatrixXi sub_cell1Ds(2, new_edges_size);
        for (unsigned int s_e = 0; s_e < new_edges_size; ++s_e)
        {
            const auto &sub_edge = split_result.Split_edges.at(split_edges.at(s_e));

            sub_cell1Ds.col(s_e) << mesh_cell0Ds_index.at(sub_edge.Extremes.at(0)),
                mesh_cell0Ds_index.at(sub_edge.Extremes.at(1));
        }

        const auto new_cell1Ds_index = Gedim::MeshUtilities::SplitCell1D(original_cell1D, sub_cell1Ds, mesh);

        for (unsigned int s_e = 0; s_e < new_edges_size; ++s_e)
            mesh_cell1Ds_index[split_edges.at(s_e)] = new_cell1Ds_index.at(s_e);
    }

    std::set<unsigned int> updated_Cell2Ds;
    std::vector<unsigned int> mesh_cell2Ds_index(split_result.Split_faces.size());
    for (unsigned int f = 0; f < split_result.Split_faces.size(); f++)
    {
        const auto &face = split_result.Split_faces.at(f);

        if (!face.IsNew)
            continue;

        mesh_cell2Ds_index[f] = mesh.Cell2DAppend(1);
        const unsigned int &newCell2DIndex = mesh_cell2Ds_index.at(f);

        Eigen::MatrixXi face_extremes(2, face.Extremes.cols());
        for (unsigned int e = 0; e < face.Extremes.cols(); e++)
        {
            face_extremes(0, e) = mesh_cell0Ds_index.at(face.Extremes(0, e));
            face_extremes(1, e) = mesh_cell1Ds_index.at(face.Extremes(1, e));
        }

        mesh.Cell2DAddVerticesAndEdges(newCell2DIndex, face_extremes);
        mesh.Cell2DSetState(newCell2DIndex, true);

        if (face.Marker < num_polyhedron_vertices + num_polyhedron_edges + 1 ||
            face.Marker > num_polyhedron_vertices + num_polyhedron_edges + num_polyhedron_faces)
        {
            mesh.Cell2DSetMarker(newCell2DIndex, mesh.Cell3DMarker(cell3D_index));
        }
        else if (face.Marker < num_polyhedron_vertices + num_polyhedron_edges + num_polyhedron_faces + 1)
        {
            const unsigned int face_pos = face.Marker - num_polyhedron_vertices - num_polyhedron_edges - 1;
            mesh.Cell2DSetMarker(newCell2DIndex, mesh.Cell2DMarker(mesh.Cell3DFace(cell3D_index, face_pos)));
        }

        cell2Ds_triangles_vertices.push_back(face.TriaglesVertices);
    }

    for (unsigned int f = 0; f < split_result.Original_face_to_split_faces.size(); f++)
    {
        const auto &split_faces = split_result.Original_face_to_split_faces.at(f);

        const unsigned int original_cell2D = mesh.Cell3DFace(cell3D_index, f);

        const unsigned int new_faces_size = split_faces.size();
        std::vector<Eigen::MatrixXi> sub_cell2Ds(new_faces_size);
        for (unsigned int s_f = 0; s_f < new_faces_size; ++s_f)
        {
            const auto &sub_face = split_result.Split_faces.at(split_faces.at(s_f));
            auto &sub_cell2D_extreme = sub_cell2Ds[s_f];

            sub_cell2D_extreme.resize(2, sub_face.Extremes.cols());
            for (unsigned int e = 0; e < sub_face.Extremes.cols(); e++)
            {
                sub_cell2D_extreme(0, e) = mesh_cell0Ds_index.at(sub_face.Extremes(0, e));
                sub_cell2D_extreme(1, e) = mesh_cell1Ds_index.at(sub_face.Extremes(1, e));
            }

            cell2Ds_triangles_vertices.push_back(sub_face.TriaglesVertices);
        }

        const auto new_cell2Ds_index = Gedim::MeshUtilities::SplitCell2D(original_cell2D, sub_cell2Ds, mesh);

        updated_Cell2Ds.insert(original_cell2D);

        for (unsigned int s_f = 0; s_f < new_faces_size; ++s_f)
            mesh_cell2Ds_index[split_faces.at(s_f)] = new_cell2Ds_index.at(s_f);
    }

    const unsigned int new_polyhedra_size = split_result.Split_polyhedrons.size();
    std::vector<std::vector<unsigned int>> sub_vertices(new_polyhedra_size);
    std::vector<std::vector<unsigned int>> sub_edges(new_polyhedra_size);
    std::vector<std::vector<unsigned int>> sub_faces(new_polyhedra_size);

    for (unsigned int s_p = 0; s_p < new_polyhedra_size; ++s_p)
    {
        const auto &split_polyhedron = split_result.Split_polyhedrons.at(s_p);

        std::vector<unsigned int> &cell_vertices = sub_vertices.at(s_p);
        std::vector<unsigned int> &cell_edges = sub_edges.at(s_p);
        std::vector<unsigned int> &cell_faces = sub_faces.at(s_p);

        cell_vertices.resize(split_polyhedron.Vertices.size());
        for (unsigned int v = 0; v < split_polyhedron.Vertices.size(); v++)
            cell_vertices[v] = mesh_cell0Ds_index.at(split_polyhedron.Vertices.at(v));

        cell_edges.resize(split_polyhedron.Edges.size());
        for (unsigned int e = 0; e < split_polyhedron.Edges.size(); e++)
            cell_edges[e] = mesh_cell1Ds_index.at(split_polyhedron.Edges.at(e));

        cell_faces.resize(split_polyhedron.Faces.size());
        for (unsigned int f = 0; f < split_polyhedron.Faces.size(); f++)
            cell_faces[f] = mesh_cell2Ds_index.at(split_polyhedron.Faces.at(f));

        unsigned int num_tetra = split_polyhedron.TetrahedronVertices.size();
        double new_cell3D_volume = 0.0;
        Eigen::MatrixXd tetra_cetroid(3, num_tetra);
        Eigen::VectorXd tetra_volume(num_tetra);
        for (unsigned int t = 0; t < num_tetra; ++t)
        {
            const auto &tetra = split_polyhedron.TetrahedronVertices.at(t);
            tetra_cetroid.col(t) = geometryUtilities.PolyhedronBarycenter(tetra);
            tetra_volume[t] = geometryUtilities.PolyhedronVolumeByInternalIntegral({tetra});
            new_cell3D_volume += tetra_volume[t];
        }

        cell3Ds_volume.push_back(new_cell3D_volume);
        cell3Ds_centroid.push_back(tetra_cetroid * tetra_volume / new_cell3D_volume);

        cell3Ds_mark.push_back(split_polyhedron.Mark);
        cell3Ds_tetra_vertices.push_back(split_polyhedron.TetrahedronVertices);
        cell3Ds_additional_constrained_points.push_back(cell3Ds_centroid.back());
    }

    const auto mesh_cell3Ds_index = Gedim::MeshUtilities::SplitCell3D(cell3D_index, sub_vertices, sub_edges, sub_faces, mesh);

    for (unsigned int e = 0; e < split_result.Original_edge_to_split_edges.size(); e++)
    {
        const auto &split_edges = split_result.Original_edge_to_split_edges.at(e);

        const unsigned int original_cell1D = mesh.Cell3DEdge(cell3D_index, e);

        std::list<unsigned int> updatedCell1Ds;
        mesh.Cell1DUpdatedCell1Ds(original_cell1D, updatedCell1Ds);

        if (split_edges.size() == 1)
            continue;

        const unsigned int original_cell1D_origin = mesh.Cell1DOrigin(original_cell1D);
        const unsigned int original_cell1D_end = mesh.Cell1DEnd(original_cell1D);

        const unsigned int new_edges_size = split_edges.size();
        std::vector<unsigned int> updated_cell1Ds(new_edges_size);
        std::vector<unsigned int> new_cell0Ds(new_edges_size + 1);

        std::vector<bool> edge_selected(new_edges_size, false);
        unsigned int new_edge_origin = original_cell1D_origin;
        unsigned int new_edge_selected = 0;

        new_cell0Ds[0] = original_cell1D_origin;
        while (new_edge_selected < new_edges_size)
        {
            for (unsigned int s_e = 0; s_e < new_edges_size; ++s_e)
            {
                if (edge_selected.at(s_e))
                    continue;

                const auto &sub_edge = split_result.Split_edges.at(split_edges.at(s_e));
                const unsigned int sub_edge_origin = mesh_cell0Ds_index.at(sub_edge.Extremes[0]);
                const unsigned int sub_edge_end = mesh_cell0Ds_index.at(sub_edge.Extremes[1]);

                if (sub_edge_origin == new_edge_origin)
                {
                    updated_cell1Ds[new_edge_selected] = mesh_cell1Ds_index.at(split_edges.at(s_e));
                    new_cell0Ds[new_edge_selected + 1] = sub_edge_end;

                    edge_selected[s_e] = true;
                    new_edge_origin = sub_edge_end;
                    new_edge_selected++;
                    break;
                }
                else if (sub_edge_end == new_edge_origin)
                {
                    updated_cell1Ds[new_edge_selected] = mesh_cell1Ds_index.at(split_edges.at(s_e));
                    new_cell0Ds[new_edge_selected + 1] = sub_edge_origin;

                    edge_selected[s_e] = true;
                    new_edge_origin = sub_edge_origin;
                    new_edge_selected++;
                    break;
                }
            }
        }
        Gedim::Output::Assert(new_cell0Ds[new_edges_size] == original_cell1D_end);

        for (unsigned int n_p = 0; n_p < mesh.Cell1DNumberNeighbourCell2D(original_cell1D); n_p++)
        {
            if (!mesh.Cell1DHasNeighbourCell2D(original_cell1D, n_p))
                continue;

            const unsigned int c_2D_n = mesh.Cell1DNeighbourCell2D(original_cell1D, n_p);

            std::list<unsigned int> updatedCell2Ds;
            mesh.Cell2DUpdatedCell2Ds(c_2D_n, updatedCell2Ds);

            if (updatedCell2Ds.size() > 1 || updatedCell2Ds.front() != c_2D_n)
                continue;

            const unsigned int localFaceEdgeIndex = mesh.Cell2DFindEdge(c_2D_n, original_cell1D);
            if (localFaceEdgeIndex == mesh.Cell2DNumberVertices(c_2D_n))
                throw std::runtime_error("Not correct edge and face");

            const unsigned int numVertices = mesh.Cell2DNumberVertices(c_2D_n);
            const std::vector<unsigned int> originalVertices = mesh.Cell2DVertices(c_2D_n);
            const std::vector<unsigned int> originalEdges = mesh.Cell2DEdges(c_2D_n);

            // create new face 2D
            Eigen::MatrixXi newFace = Eigen::MatrixXi::Zero(2, numVertices + new_edges_size - 1);

            for (unsigned int e = 0; e < localFaceEdgeIndex; e++)
                newFace.col(e) << originalVertices.at(e), originalEdges.at(e);

            const unsigned int origin_new_edge = mesh.Cell2DVertex(c_2D_n, localFaceEdgeIndex);

            if (mesh.Cell1DOrigin(updated_cell1Ds.at(0)) == origin_new_edge || mesh.Cell1DEnd(updated_cell1Ds.at(0)) == origin_new_edge)
            {
                for (unsigned int s_e = 0; s_e < new_edges_size; s_e++)
                    newFace.col(localFaceEdgeIndex + s_e) << new_cell0Ds.at(s_e), updated_cell1Ds.at(s_e);
            }
            else if (mesh.Cell1DOrigin(updated_cell1Ds.at(new_edges_size - 1)) == origin_new_edge ||
                     mesh.Cell1DEnd(updated_cell1Ds.at(new_edges_size - 1)) == origin_new_edge)
            {
                for (unsigned int s_e = 0; s_e < new_edges_size; s_e++)
                    newFace.col(localFaceEdgeIndex + s_e) << new_cell0Ds.at(new_edges_size - s_e),
                        updated_cell1Ds.at(new_edges_size - s_e - 1);
            }
            else
                throw std::runtime_error("Not managed edge split two");

            for (unsigned int e = localFaceEdgeIndex + 1; e < numVertices; e++)
                newFace.col(e + new_edges_size - 1) << originalVertices.at(e), originalEdges.at(e);

            const auto newCell2DIndices = Gedim::MeshUtilities::SplitCell2D(c_2D_n, {newFace}, mesh);

            cell2Ds_triangles_vertices.push_back(cell2Ds_triangles_vertices.at(c_2D_n));

            Gedim::Output::Assert(newCell2DIndices.size() == 1);
            updated_Cell2Ds.insert(c_2D_n);
        }
    }

    for (const unsigned int original_cell2D : updated_Cell2Ds)
    {
        for (unsigned int n_p = 0; n_p < mesh.Cell2DNumberNeighbourCell3D(original_cell2D); n_p++)
        {
            if (!mesh.Cell2DHasNeighbourCell3D(original_cell2D, n_p))
                continue;

            const unsigned int c_3D_n = mesh.Cell2DNeighbourCell3D(original_cell2D, n_p);

            if (c_3D_n == cell3D_index)
                continue;

            std::list<unsigned int> updatedCell3Ds;
            mesh.Cell3DUpdatedCell3Ds(c_3D_n, updatedCell3Ds);

            if (updatedCell3Ds.size() > 1 || updatedCell3Ds.front() != c_3D_n)
                continue;

            std::unordered_set<unsigned int> new_cell3D_n_vertices;
            std::unordered_set<unsigned int> new_cell3D_n_edges;
            std::unordered_set<unsigned int> new_cell3D_n_faces;

            const auto &original_cell3D_n_faces = mesh.Cell3DFaces(c_3D_n);

            for (const unsigned int c_2D : original_cell3D_n_faces)
            {
                std::list<unsigned int> updatedCell2Ds;
                mesh.Cell2DUpdatedCell2Ds(c_2D, updatedCell2Ds);

                for (const unsigned int new_c_2D : updatedCell2Ds)
                {
                    if (new_cell3D_n_faces.find(new_c_2D) != new_cell3D_n_faces.end())
                        continue;

                    new_cell3D_n_faces.insert(new_c_2D);

                    for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(new_c_2D); v++)
                    {
                        const unsigned int new_c0D = mesh.Cell2DVertex(new_c_2D, v);

                        if (new_cell3D_n_vertices.find(new_c0D) != new_cell3D_n_vertices.end())
                            continue;

                        new_cell3D_n_vertices.insert(new_c0D);
                    }

                    for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(new_c_2D); e++)
                    {
                        const unsigned int new_c1D = mesh.Cell2DEdge(new_c_2D, e);

                        if (new_cell3D_n_edges.find(new_c1D) != new_cell3D_n_edges.end())
                            continue;

                        new_cell3D_n_edges.insert(new_c1D);
                    }
                }
            }

            Gedim::MeshUtilities::SplitCell3D(
                c_3D_n,
                {std::vector<unsigned int>(new_cell3D_n_vertices.begin(), new_cell3D_n_vertices.end())},
                {std::vector<unsigned int>(new_cell3D_n_edges.begin(), new_cell3D_n_edges.end())},
                {std::vector<unsigned int>(new_cell3D_n_faces.begin(), new_cell3D_n_faces.end())},
                mesh);

            cell3Ds_volume.push_back(cell3Ds_volume.at(c_3D_n));
            cell3Ds_centroid.push_back(cell3Ds_centroid.at(c_3D_n));
            cell3Ds_mark.push_back(cell3Ds_mark.at(c_3D_n));

            cell3Ds_tetra_vertices.push_back(cell3Ds_tetra_vertices.at(c_3D_n));
            cell3Ds_additional_constrained_points.push_back(cell3Ds_additional_constrained_points.at(c_3D_n));
        }
    }

    cell2Ds_triangles_vertices.shrink_to_fit();

    cell3Ds_volume.shrink_to_fit();
    cell3Ds_centroid.shrink_to_fit();
    cell3Ds_mark.shrink_to_fit();
    cell3Ds_tetra_vertices.shrink_to_fit();
    cell3Ds_additional_constrained_points.shrink_to_fit();

#ifdef EXPORT_TEST_VTU
    {
        std::string mesh_3D_export_folder = "./TEST/Cell3D_" + std::to_string(cell3D_index);
        Gedim::Output::CreateFolder(mesh_3D_export_folder);
        Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, mesh_3D_export_folder, "Almost_Mesh3D");
    }
#endif

    return true;
}
// ***************************************************************************
AgglomerateMeshByFunctionResult AgglomerateMeshByFunction(const Gedim::GeometryUtilities &geometryUtilities,
                                                          const Gedim::GraphUtilities &graphUtilities,
                                                          Gedim::MeshMatrices &mesh_data,
                                                          const std::function<Eigen::VectorXi(const Eigen::MatrixXd &)> &marking_function,
                                                          const double &resolution,
                                                          const std::vector<Eigen::MatrixXd> &cell3Ds_constrained_points)
{
    Gedim::MeshMatricesDAO mesh(mesh_data);
    std::list<unsigned int> to_agglomerate_cell3Ds_lst;

    Gedim::MeshUtilities::ComputeCell2DCell3DNeighbours(mesh);

    const auto mesh_geometric_data = Gedim::MeshUtilities::FillMesh3DGeometricData(geometryUtilities, mesh);
    std::vector<Eigen::Vector3d> cell3Ds_centroid(mesh.Cell3DTotalNumber(), Eigen::Vector3d::Zero());
    std::vector<double> cell3Ds_volume(mesh.Cell3DTotalNumber(), 0.0);
    std::vector<unsigned int> cell3Ds_mark(mesh.Cell3DTotalNumber(), 0);
    std::vector<std::vector<Eigen::MatrixXd>> cell3Ds_tetra_vertices(mesh.Cell3DTotalNumber());
    std::vector<std::vector<Eigen::Matrix3d>> cell2Ds_triangles_vertices(mesh.Cell2DTotalNumber());
    std::vector<Eigen::MatrixXd> cell3Ds_additional_constrained_points(mesh.Cell3DTotalNumber());

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
        if (!mesh.Cell3DIsActive(c))
            continue;

        const auto &cell3D_vertices = mesh_geometric_data.Cell3DsVertices.at(c);
        const auto &cell3D_centroid = mesh_geometric_data.Cell3DsCentroids.at(c);
        const auto &cell3D_faces_centroid_2D = mesh_geometric_data.Cell3DsFaces2DCentroids.at(c);
        const auto &cell3D_constrained_points = cell3Ds_constrained_points.empty() ? Eigen::MatrixXd()
                                                                                   : cell3Ds_constrained_points.at(c);
        Eigen::MatrixXd cell3D_faces_centroid_3D(3, cell3D_faces_centroid_2D.size());

        cell3Ds_volume[c] = mesh_geometric_data.Cell3DsVolumes.at(c);
        cell3Ds_centroid[c] = cell3D_centroid;
        cell3Ds_tetra_vertices[c] = mesh_geometric_data.Cell3DsTetrahedronPoints.at(c);

        for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(c); f++)
        {
            const unsigned int c2D_index = mesh.Cell3DFace(c, f);
            cell2Ds_triangles_vertices[c2D_index] = mesh_geometric_data.Cell3DsFaces3DTriangulations.at(c).at(f);
        }

        if (cell3D_constrained_points.size() > 0)
        {
            cell3Ds_additional_constrained_points[c] = cell3D_constrained_points;
            to_agglomerate_cell3Ds_lst.push_back(c);
            continue;
        }
        else
            cell3Ds_additional_constrained_points[c] = cell3D_centroid;

        for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(c); f++)
        {
            cell3D_faces_centroid_3D.col(f) =
                geometryUtilities.RotatePointsFrom2DTo3D(cell3D_faces_centroid_2D.at(f),
                                                         mesh_geometric_data.Cell3DsFacesRotationMatrices.at(c).at(f),
                                                         mesh_geometric_data.Cell3DsFacesTranslations.at(c).at(f));
        }

        Eigen::MatrixXd polyhedron_test_points(3, cell3D_vertices.cols() + 1 + cell3D_faces_centroid_3D.cols());

        polyhedron_test_points.block(0, 0, 3, cell3D_vertices.cols()) = cell3D_vertices;
        polyhedron_test_points.block(0, cell3D_vertices.cols(), 3, 1) = cell3D_centroid;
        polyhedron_test_points.block(0, cell3D_vertices.cols() + 1, 3, cell3D_faces_centroid_3D.cols()) = cell3D_faces_centroid_3D;

        const auto vertices_markers = marking_function(polyhedron_test_points);

        if (vertices_markers.sum() == 0)
            continue;

        if (vertices_markers.sum() == polyhedron_test_points.cols())
        {
            cell3Ds_mark[c] = 1;
            continue;
        }

        to_agglomerate_cell3Ds_lst.push_back(c);
    }

    std::vector<unsigned int> to_agglomerate_cell3Ds(to_agglomerate_cell3Ds_lst.begin(), to_agglomerate_cell3Ds_lst.end());

    for (unsigned int c = 0; c < to_agglomerate_cell3Ds.size(); c++)
    {
        const unsigned int cell3D_index = to_agglomerate_cell3Ds.at(c);

        Gedim::MeshUtilities::ComputeCell1DCell2DNeighbours(mesh);
        Gedim::MeshUtilities::ComputeCell1DCell3DNeighbours(mesh);
        Gedim::MeshUtilities::ComputeCell2DCell3DNeighbours(mesh);

        // std::cout << "Split " << cell3D_index << " - " << c << "/" << to_agglomerate_cell3Ds.size() << std::endl;

        const auto cell3D_polyhedron = Gedim::MeshUtilities::MeshCell3DToPolyhedron(mesh, cell3D_index);

        const auto agglomeration_result = AgglomerateMeshCellByFunction(geometryUtilities,
                                                                        graphUtilities,
                                                                        mesh,
                                                                        cell3D_polyhedron.Vertices,
                                                                        cell3D_polyhedron.Edges,
                                                                        cell3D_polyhedron.Faces,
                                                                        cell3Ds_additional_constrained_points,
                                                                        cell3Ds_centroid,
                                                                        cell3Ds_volume,
                                                                        cell3Ds_mark,
                                                                        cell3Ds_tetra_vertices,
                                                                        cell2Ds_triangles_vertices,
                                                                        cell3D_index,
                                                                        marking_function,
                                                                        resolution);
        if (!agglomeration_result)
            continue;

#ifdef CHECK_MESH
        Gedim::MeshUtilities::CheckMesh3DConfiguration check;
        check.Cell2D_CheckConvexity = false;
        check.Cell3D_CheckConvexity = false;
        check.Cell2D_CheckMeasure = false;
        check.Cell3D_CheckMeasure = false;
        meshUtilities.CheckMesh3D(check, geometryUtilities, mesh);
#endif

        {
            Gedim::MeshMatrices temp_mesh_data = mesh_data;
            Gedim::MeshMatricesDAO temp_mesh(temp_mesh_data);
            mesh_data = Gedim::MeshMatrices();

            const auto temp_centroid = cell3Ds_centroid;
            const auto temp_volume = cell3Ds_volume;
            const auto temp_mark = cell3Ds_mark;
            const auto temp_tetra_vertices = cell3Ds_tetra_vertices;
            const auto temp_triangles_vertices = cell2Ds_triangles_vertices;
            const auto temp_additional_constrained_points = cell3Ds_additional_constrained_points;

            const auto filter_active_mesh = Gedim::MeshUtilities::FilterActiveMesh(temp_mesh);
            const auto extract_mesh_data = Gedim::MeshUtilities::ExtractMesh3D(filter_active_mesh.Cell0Ds,
                                                                               filter_active_mesh.Cell1Ds,
                                                                               filter_active_mesh.Cell2Ds,
                                                                               filter_active_mesh.Cell3Ds,
                                                                               temp_mesh,
                                                                               mesh);

            cell3Ds_volume.resize(mesh.Cell3DTotalNumber());
            cell3Ds_centroid.resize(mesh.Cell3DTotalNumber());
            cell3Ds_mark.resize(mesh.Cell3DTotalNumber());
            cell3Ds_tetra_vertices.resize(mesh.Cell3DTotalNumber());
            cell3Ds_additional_constrained_points.resize(mesh.Cell3DTotalNumber());
            for (unsigned int c_3D = 0; c_3D < mesh.Cell3DTotalNumber(); ++c_3D)
            {
                const auto old_c3D = extract_mesh_data.NewCell3DToOldCell3D.at(c_3D);

                cell3Ds_volume.at(c_3D) = temp_volume.at(old_c3D);
                cell3Ds_centroid.at(c_3D) = temp_centroid.at(old_c3D);
                cell3Ds_mark.at(c_3D) = temp_mark.at(old_c3D);
                cell3Ds_tetra_vertices.at(c_3D) = temp_tetra_vertices.at(old_c3D);
                cell3Ds_additional_constrained_points.at(c_3D) = temp_additional_constrained_points.at(old_c3D);
            }

            cell2Ds_triangles_vertices.resize(mesh.Cell2DTotalNumber());
            for (unsigned int c_2D = 0; c_2D < mesh.Cell2DTotalNumber(); ++c_2D)
            {
                const auto old_c2D = extract_mesh_data.NewCell2DToOldCell2D.at(c_2D);

                cell2Ds_triangles_vertices.at(c_2D) = temp_triangles_vertices.at(old_c2D);
            }

            for (unsigned int n_c = c + 1; n_c < to_agglomerate_cell3Ds.size(); n_c++)
            {
                const unsigned int original_cell3D_index = to_agglomerate_cell3Ds.at(n_c);

                std::list<unsigned int> updatedCell3Ds;
                temp_mesh.Cell3DUpdatedCell3Ds(original_cell3D_index, updatedCell3Ds);

                Gedim::Output::Assert(updatedCell3Ds.size() == 1);
                to_agglomerate_cell3Ds[n_c] = extract_mesh_data.OldCell3DToNewCell3D.at(updatedCell3Ds.front());
            }
        }
    }

    return {cell3Ds_mark, cell3Ds_tetra_vertices, cell2Ds_triangles_vertices};
}
// ***************************************************************************
void CreateMeshWithHoles(const Gedim::GeometryUtilities &geometry_utilities,
                         const std::vector<Gedim::GeometryUtilities::Polyhedron> &domains,
                         Gedim::MeshMatrices &full_mesh_data,
                         Gedim::IMeshDAO &full_mesh,
                         Gedim::IMeshDAO &mesh)
{
    Gedim::GraphUtilities graphUtilities;

    std::vector<Eigen::MatrixXd> cell2Ds_vertices(full_mesh.Cell2DTotalNumber());
    std::vector<Eigen::MatrixXd> cell2Ds_2D_vertices(full_mesh.Cell2DTotalNumber());
    std::vector<Eigen::Vector3d> cell2Ds_normal(full_mesh.Cell2DTotalNumber());
    std::vector<Eigen::Vector3d> cell2Ds_translation(full_mesh.Cell2DTotalNumber());
    std::vector<Eigen::Matrix3d> cell2Ds_rotation_matrix(full_mesh.Cell2DTotalNumber());
    std::vector<Eigen::MatrixXd> cell2Ds_bounding_box(full_mesh.Cell2DTotalNumber());
    for (unsigned int p = 0; p < full_mesh.Cell2DTotalNumber(); p++)
    {
        cell2Ds_vertices[p] = full_mesh.Cell2DVerticesCoordinates(p);
        cell2Ds_normal[p] = geometry_utilities.PolygonNormal(cell2Ds_vertices[p]);
        cell2Ds_bounding_box[p] = geometry_utilities.PointsBoundingBox(cell2Ds_vertices[p]);
        cell2Ds_translation[p] = geometry_utilities.PolygonTranslation(cell2Ds_vertices[p]);
        cell2Ds_rotation_matrix[p] =
            geometry_utilities.PolygonRotationMatrix(cell2Ds_vertices[p], cell2Ds_normal[p], cell2Ds_translation[p]);
        cell2Ds_2D_vertices[p] = geometry_utilities.RotatePointsFrom3DTo2D(cell2Ds_vertices[p],
                                                                           cell2Ds_rotation_matrix[p].transpose(),
                                                                           cell2Ds_translation[p]);
    }

    Gedim::MeshUtilities::ComputeCell2DCell3DNeighbours(full_mesh);

    const auto cell1Ds_geometric_data = Gedim::MeshUtilities::FillMesh1DGeometricData(geometry_utilities, full_mesh);

    const auto cell3Ds_geometric_data = Gedim::MeshUtilities::FillMesh3DGeometricData(geometry_utilities, full_mesh);

    struct polyhedron_geometric_data final
    {
        std::vector<Eigen::MatrixXd> Edges_vertices;
        std::vector<Eigen::MatrixXd> Edges_bounding_box;
        Eigen::MatrixXd Edges_tangent;
        std::vector<Eigen::MatrixXd> Faces_vertices;
        std::vector<Eigen::MatrixXd> Faces_bounding_box;
        Eigen::Vector3d Internal_point;
        std::vector<Eigen::Vector3d> Faces_translation;
        std::vector<Eigen::Vector3d> Faces_normal;
        std::vector<bool> Faces_normal_direction;
        std::vector<Eigen::Matrix3d> Faces_rotation_matrix;
        std::vector<Eigen::MatrixXd> Faces_rotated_vertices;
        Eigen::MatrixXd Bounding_box;
    };

    const unsigned int num_holes = domains.size() - 1;

    std::vector<Gedim::GeometryUtilities::Polyhedron> polyhedrons(num_holes);
    for (unsigned int pol = 0; pol < num_holes; ++pol)
    {
        polyhedrons[pol].Vertices = domains[pol + 1].Vertices;
        polyhedrons[pol].Edges = domains[pol + 1].Edges;
        polyhedrons[pol].Faces = domains[pol + 1].Faces;
    }

    std::vector<std::list<Eigen::Vector3d>> cell3Ds_intersections_list(full_mesh.Cell3DTotalNumber());
    std::vector<polyhedron_geometric_data> polyhedrons_geometric_data(domains.size() - 1);

    for (unsigned int pol = 0; pol < num_holes; ++pol)
    {
        const auto &polyhedron = polyhedrons.at(pol);
        auto &polyhedron_geometric_data = polyhedrons_geometric_data[pol];

        polyhedron_geometric_data.Edges_vertices.resize(polyhedron.Edges.cols());
        polyhedron_geometric_data.Edges_bounding_box.resize(polyhedron.Edges.cols());
        for (unsigned int e = 0; e < polyhedron.Edges.cols(); e++)
        {
            auto &polyhedron_edges_vertices = polyhedron_geometric_data.Edges_vertices[e];
            auto &polyhedron_edges_bounding_box = polyhedron_geometric_data.Edges_bounding_box[e];

            polyhedron_edges_vertices.resize(3, 2);
            polyhedron_edges_vertices.col(0) << polyhedron.Vertices.col(polyhedron.Edges(0, e));
            polyhedron_edges_vertices.col(1) << polyhedron.Vertices.col(polyhedron.Edges(1, e));

            polyhedron_edges_bounding_box = geometry_utilities.PointsBoundingBox(polyhedron_edges_vertices);
        }

        polyhedron_geometric_data.Edges_tangent =
            geometry_utilities.PolyhedronEdgeTangents(polyhedron.Vertices, polyhedron.Edges);

        polyhedron_geometric_data.Faces_vertices =
            geometry_utilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);

        polyhedron_geometric_data.Faces_bounding_box.resize(polyhedron.Faces.size());
        for (unsigned int f = 0; f < polyhedron.Faces.size(); f++)
        {
            auto &polyhdron_faces_vertices = polyhedron_geometric_data.Faces_vertices[f];
            auto &polyhedron_faces_bounding_box = polyhedron_geometric_data.Faces_bounding_box[f];
            polyhedron_faces_bounding_box = geometry_utilities.PointsBoundingBox(polyhdron_faces_vertices);
        }

        polyhedron_geometric_data.Internal_point = geometry_utilities.PolyhedronBarycenter(polyhedron.Vertices);
        polyhedron_geometric_data.Faces_translation =
            geometry_utilities.PolyhedronFaceTranslations(polyhedron_geometric_data.Faces_vertices);
        polyhedron_geometric_data.Faces_normal = geometry_utilities.PolyhedronFaceNormals(polyhedron_geometric_data.Faces_vertices);
        polyhedron_geometric_data.Faces_normal_direction =
            geometry_utilities.PolyhedronFaceNormalDirections(polyhedron_geometric_data.Faces_vertices,
                                                              polyhedron_geometric_data.Internal_point,
                                                              polyhedron_geometric_data.Faces_normal);
        polyhedron_geometric_data.Faces_rotation_matrix =
            geometry_utilities.PolyhedronFaceRotationMatrices(polyhedron_geometric_data.Faces_vertices,
                                                              polyhedron_geometric_data.Faces_normal,
                                                              polyhedron_geometric_data.Faces_translation);
        polyhedron_geometric_data.Faces_rotated_vertices =
            geometry_utilities.PolyhedronFaceRotatedVertices(polyhedron_geometric_data.Faces_vertices,
                                                             polyhedron_geometric_data.Faces_translation,
                                                             polyhedron_geometric_data.Faces_rotation_matrix);
        polyhedron_geometric_data.Bounding_box = geometry_utilities.PointsBoundingBox(polyhedron.Vertices);

        const auto intersections_mesh_polyedron =
            Gedim::MeshUtilities::Intersect_mesh_polyhedron(geometry_utilities,
                                                            polyhedron.Vertices,
                                                            polyhedron.Edges,
                                                            polyhedron_geometric_data.Edges_vertices,
                                                            polyhedron_geometric_data.Edges_tangent,
                                                            polyhedron_geometric_data.Edges_bounding_box,
                                                            polyhedron.Faces,
                                                            polyhedron_geometric_data.Faces_vertices,
                                                            polyhedron_geometric_data.Faces_rotated_vertices,
                                                            polyhedron_geometric_data.Faces_normal,
                                                            polyhedron_geometric_data.Faces_normal_direction,
                                                            polyhedron_geometric_data.Faces_translation,
                                                            polyhedron_geometric_data.Faces_rotation_matrix,
                                                            polyhedron_geometric_data.Faces_bounding_box,
                                                            polyhedron_geometric_data.Bounding_box,
                                                            full_mesh,
                                                            cell1Ds_geometric_data.Cell1DsBoundingBox,
                                                            cell1Ds_geometric_data.Cell1DsVertices,
                                                            cell1Ds_geometric_data.Cell1DsTangents,
                                                            cell2Ds_vertices,
                                                            cell2Ds_normal,
                                                            cell2Ds_2D_vertices,
                                                            cell2Ds_translation,
                                                            cell2Ds_rotation_matrix,
                                                            cell2Ds_bounding_box,
                                                            cell3Ds_geometric_data.Cell3DsBoundingBox,
                                                            cell3Ds_geometric_data.Cell3DsFaces,
                                                            cell3Ds_geometric_data.Cell3DsFaces3DVertices,
                                                            cell3Ds_geometric_data.Cell3DsFaces2DVertices,
                                                            cell3Ds_geometric_data.Cell3DsFacesNormals,
                                                            cell3Ds_geometric_data.Cell3DsFacesNormalDirections,
                                                            cell3Ds_geometric_data.Cell3DsFacesTranslations,
                                                            cell3Ds_geometric_data.Cell3DsFacesRotationMatrices);

        for (unsigned int c = 0; c < full_mesh.Cell3DTotalNumber(); ++c)
        {
            if (!full_mesh.Cell3DIsActive(c))
                continue;

            std::set<unsigned int> cell3D_intersections_index;
            auto &cell3D_intersections_list = cell3Ds_intersections_list.at(c);

            for (unsigned int e = 0; e < full_mesh.Cell3DNumberEdges(c); ++e)
            {
                const unsigned int cell1D_index = full_mesh.Cell3DEdge(c, e);

                const auto cell1D_intersections =
                    intersections_mesh_polyedron.Mesh_intersections.Cell1Ds_intersections.find(cell1D_index);

                if (cell1D_intersections == intersections_mesh_polyedron.Mesh_intersections.Cell1Ds_intersections.end())
                    continue;

                if (cell1D_intersections->second.empty())
                    continue;

                const unsigned int num_intersections = cell1D_intersections->second.size();

                for (unsigned int i = 0; i < num_intersections; i++)
                {
                    const unsigned int intersection_index = cell1D_intersections->second[i];

                    if (cell3D_intersections_index.find(intersection_index) != cell3D_intersections_index.end())
                        continue;

                    cell3D_intersections_index.insert(intersection_index);
                    cell3D_intersections_list.push_back(intersections_mesh_polyedron.Intersections_Coordinates.col(intersection_index));
                }
            }

            for (unsigned int f = 0; f < full_mesh.Cell3DNumberFaces(c); ++f)
            {
                const unsigned int cell2D_index = full_mesh.Cell3DFace(c, f);

                const auto cell2D_intersections =
                    intersections_mesh_polyedron.Mesh_intersections.Cell2Ds_intersections.find(cell2D_index);

                if (cell2D_intersections == intersections_mesh_polyedron.Mesh_intersections.Cell2Ds_intersections.end())
                    continue;

                if (cell2D_intersections->second.empty())
                    continue;

                const unsigned int num_intersections = cell2D_intersections->second.size();

                for (unsigned int i = 0; i < num_intersections; i++)
                {
                    const unsigned int intersection_index = cell2D_intersections->second[i];

                    if (cell3D_intersections_index.find(intersection_index) != cell3D_intersections_index.end())
                        continue;

                    cell3D_intersections_index.insert(intersection_index);
                    cell3D_intersections_list.push_back(intersections_mesh_polyedron.Intersections_Coordinates.col(intersection_index));
                }
            }

            const auto cell3D_intersections = intersections_mesh_polyedron.Mesh_intersections.Cell3Ds_intersections.find(c);

            if (cell3D_intersections == intersections_mesh_polyedron.Mesh_intersections.Cell3Ds_intersections.end())
                continue;

            if (cell3D_intersections->second.empty())
                continue;

            const unsigned int num_intersections = cell3D_intersections->second.size();

            for (unsigned int i = 0; i < num_intersections; i++)
            {
                const unsigned int intersection_index = cell3D_intersections->second[i];

                if (cell3D_intersections_index.find(intersection_index) != cell3D_intersections_index.end())
                    continue;

                cell3D_intersections_index.insert(intersection_index);
                cell3D_intersections_list.push_back(intersections_mesh_polyedron.Intersections_Coordinates.col(intersection_index));
            }
        }
    }

    std::vector<Eigen::MatrixXd> cell3Ds_constrained_points(full_mesh.Cell3DTotalNumber());
    for (unsigned int c = 0; c < full_mesh.Cell3DTotalNumber(); ++c)
    {
        if (!full_mesh.Cell3DIsActive(c))
            continue;

        cell3Ds_constrained_points[c] = Gedim::Eigen_Utilities::CoordinatesToMatrix(cell3Ds_intersections_list.at(c));
    }

    auto empty_function = [&](const Eigen::MatrixXd &points) {
        Eigen::VectorXi marks(points.cols());

        for (unsigned int p = 0; p < points.cols(); p++)
        {
            const Eigen::Vector3d point = points.col(p);

            marks[p] = 0;

            for (unsigned int pol = 0; pol < polyhedrons.size(); ++pol)
            {
                const auto &polyhedron = polyhedrons.at(pol);
                const auto &polyhedron_geometric_data = polyhedrons_geometric_data.at(pol);

                if (!geometry_utilities.IsPointInBoundingBox(point, polyhedron_geometric_data.Bounding_box))
                {
                    continue;
                }

                const auto point_polyhedron_position =
                    geometry_utilities.PointPolyhedronPosition(point,
                                                               polyhedron.Faces,
                                                               polyhedron_geometric_data.Faces_vertices,
                                                               polyhedron_geometric_data.Faces_rotated_vertices,
                                                               polyhedron_geometric_data.Faces_normal,
                                                               polyhedron_geometric_data.Faces_normal_direction,
                                                               polyhedron_geometric_data.Faces_translation,
                                                               polyhedron_geometric_data.Faces_rotation_matrix);

                switch (point_polyhedron_position.Type)
                {
                case Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::Outside:
                    continue;
                case Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace:
                case Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge:
                case Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex:
                    continue;
                case Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::Inside:
                    break;
                default:
                    throw std::runtime_error("Unexpected point intersection");
                }

                marks[p] = 1;
                break;
            }
        }

        return marks;
    };

    const auto agglomeration_info =
        AgglomerateMeshByFunction(geometry_utilities, graphUtilities, full_mesh_data, empty_function, 0.0, cell3Ds_constrained_points);

    Gedim::MeshUtilities::CheckMesh3DConfiguration checkMeshDataConfig;
    checkMeshDataConfig.Cell2D_CheckConvexity = false;
    checkMeshDataConfig.Cell3D_CheckConvexity = false;
    checkMeshDataConfig.Cell2D_CheckMeasure = false;
    checkMeshDataConfig.Cell3D_CheckMeasure = false;

    Gedim::MeshUtilities::ComputeCell2DCell3DNeighbours(full_mesh);

    const auto mesh_geometric_data = Gedim::MeshUtilities::FillMesh3DGeometricData(geometry_utilities,
                                                                                   full_mesh,
                                                                                   agglomeration_info.Cell3Ds_tetra_vertices,
                                                                                   agglomeration_info.Cell2Ds_triangles_vertices);

    Gedim::MeshUtilities::ComputeCell2DCell3DNeighbours(full_mesh);
    Gedim::MeshUtilities::ComputeCell1DCell3DNeighbours(full_mesh);
    Gedim::MeshUtilities::ComputeCell0DCell3DNeighbours(full_mesh);

    const unsigned int numCell3D = full_mesh.Cell3DTotalNumber();

    for (unsigned int c = 0; c < numCell3D; c++)
        if (agglomeration_info.Cell3Ds_mark[c] == 1)
            full_mesh.Cell3DSetState(c, false);

    const unsigned int markerObstacles = domains[0].Vertices.cols() + domains[0].Edges.cols() + domains[0].Faces.size() + 1;

    for (unsigned int c = 0; c < numCell3D; c++)
    {
        if (agglomeration_info.Cell3Ds_mark[c] == 1)
        {
            const std::vector<unsigned int> faces = full_mesh.Cell3DFaces(c);

            for (unsigned int f : faces)
            {
                if (!full_mesh.Cell2DIsActive(f))
                    continue;

                unsigned int count = 0;
                const std::vector<unsigned int> neighbours = full_mesh.Cell2DNeighbourCell3Ds(f);

                for (unsigned int n : neighbours)
                    if (agglomeration_info.Cell3Ds_mark[n] == 1)
                    {
                        count++;
                    }

                if (count == neighbours.size())
                    full_mesh.Cell2DSetState(f, false);
                else if (count == 1)
                    full_mesh.Cell2DSetMarker(f, markerObstacles);
                else
                    throw std::runtime_error("Not valid configuration");
            }
        }
    }

    for (unsigned int c = 0; c < numCell3D; c++)
    {
        if (agglomeration_info.Cell3Ds_mark[c] == 1)
        {
            const std::vector<unsigned int> edges = full_mesh.Cell3DEdges(c);

            for (unsigned int e : edges)
            {
                if (!full_mesh.Cell1DIsActive(e))
                    continue;

                unsigned int count = 0;
                const std::vector<unsigned int> neighbours = full_mesh.Cell1DNeighbourCell3Ds(e);

                for (unsigned int n : neighbours)
                    if (agglomeration_info.Cell3Ds_mark[n] == 1)
                    {
                        count++;
                    }

                if (count == neighbours.size())
                    full_mesh.Cell1DSetState(e, false);
                else
                    full_mesh.Cell1DSetMarker(e, markerObstacles + 1);
            }
        }
    }

    std::vector<unsigned int> tempNearObstacleIds;
    for (unsigned int c = 0; c < numCell3D; c++)
    {
        if (agglomeration_info.Cell3Ds_mark[c] == 1)
        {
            const std::vector<unsigned int> vertices = full_mesh.Cell3DVertices(c);

            for (unsigned int v : vertices)
            {
                if (!full_mesh.Cell0DIsActive(v))
                    continue;

                unsigned int count = 0;
                const std::vector<unsigned int> neighbours = full_mesh.Cell0DNeighbourCell3Ds(v);

                for (unsigned int n : neighbours)
                    if (agglomeration_info.Cell3Ds_mark[n] == 1)
                    {
                        count++;
                    }

                if (count == neighbours.size())
                    full_mesh.Cell0DSetState(v, false);
                else
                {
                    for (unsigned int n : neighbours)
                        if (agglomeration_info.Cell3Ds_mark[n] != 1)
                        {
                            tempNearObstacleIds.push_back(n);
                        }

                    full_mesh.Cell0DSetMarker(v, markerObstacles + 1);
                }
            }
        }
    }

    const auto filter_active_mesh = Gedim::MeshUtilities::FilterActiveMesh(full_mesh);
    const auto extract_mesh_data = Gedim::MeshUtilities::ExtractMesh3D(filter_active_mesh.Cell0Ds,
                                                                       filter_active_mesh.Cell1Ds,
                                                                       filter_active_mesh.Cell2Ds,
                                                                       filter_active_mesh.Cell3Ds,
                                                                       full_mesh,
                                                                       mesh);
}
// ***************************************************************************
} // namespace MeshUtilities
} // namespace External
} // namespace Gedim
