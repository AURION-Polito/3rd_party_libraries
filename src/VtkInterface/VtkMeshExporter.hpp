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

#ifndef __VtkMeshExpertor_H
#define __VtkMeshExpertor_H

#include "3rd_party_libraries_Macro.hpp"

#include <vector>

#include "Eigen/Eigen"
#include "GeometryUtilities.hpp"
#include "IMeshDAO.hpp"

#if ENABLE_VTK == 1
#include "vtkAppendFilter.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCellTypes.h"
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkWeakPointer.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#endif // ENABLE_VTK

namespace Gedim
{
namespace External
{
namespace MeshUtilities
{
/// \brief Export Mesh To VTU
/// \param mesh the mesh
/// \param exportFolder the folder in which the mesh is exported
void ExportMeshToVTU(const Gedim::IMeshDAO &mesh, const std::string &exportFolder, const std::string &fileName, const bool &separateFile = false);

/// \brief Export Cell2D To VTU
/// \param mesh the mesh
/// \param cell2DIndex the cell2D index
/// \param cell2DVertices the cell2D vertices
/// \param cell2DTriangulations the cell2D triangulation
/// \param cell2DArea the cell2D area
/// \param cell2DCentroid the cell2D centroid
/// \param exportFolder the folder in which to export
void ExportCell2DToVTU(const Gedim::IMeshDAO &mesh,
                       const unsigned int &cell2DIndex,
                       const Eigen::MatrixXd &cell2DVertices,
                       const std::vector<Eigen::Matrix3d> &cell2DTriangulations,
                       const double &cell2DArea,
                       const Eigen::Vector3d &cell2DCentroid,
                       const std::string &exportFolder);

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
                       const std::string &exportFolder);
} // namespace MeshUtilities
} // namespace External
} // namespace Gedim

namespace Gedim
{
namespace External
{
namespace GeometryUtilities
{

/// \brief Export Polyhedron To VTU
/// \param polyhedronVertices the polyhedron vertices
/// \param polyhedronEdges the polyhedron edges
/// \param polyhedronFaces the polyhedron faces
/// \param exportFolder the folder in which to export
void ExportPolyhedronToVTU(const Eigen::MatrixXd &polyhedronVertices,
                           const Eigen::MatrixXi &polyhedronEdges,
                           const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                           const std::string &exportFolder);

inline void ExportPolyhedronToVTU(const Gedim::GeometryUtilities::Polyhedron &polyhedron, const std::string &exportFolder)
{
    Gedim::External::GeometryUtilities::ExportPolyhedronToVTU(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, exportFolder);
}

/// \brief Export Polyhedron To VTU
/// \param polyhedronVertices the polyhedron vertices
/// \param polyhedronEdges the polyhedron edges
/// \param polyhedronFaces the polyhedron faces
/// \param exportFolder the folder in which to export
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
                           const std::string &exportFolder);

void ExportPolygonToVTU(const unsigned int &index,
                        const Eigen::MatrixXd &polygon,
                        const std::vector<Eigen::Matrix3d> &polygon_triangles,
                        const double &polygon_volume,
                        const Eigen::Vector3d &polygon_centroid,
                        const Eigen::MatrixXd &polygon_edges_centroid,
                        const Eigen::MatrixXd &polygon_edges_normal,
                        const std::vector<bool> &polygon_edges_normal_direction,
                        const std::string &exportFolder);

void ExportPlanteToVTU(const Gedim::GeometryUtilities &geometryUtilities,
                       const Eigen::Vector3d &plane_origin,
                       const Eigen::Vector3d &plane_normal,
                       const Eigen::Vector3d &plane_translation,
                       const Eigen::Matrix3d &plane_rotation_matrix,
                       const std::string &exportFolder);

} // namespace GeometryUtilities
} // namespace External
} // namespace Gedim

#endif //__VtkMeshExpertor_H
