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

#ifndef __TriangleInterface_H
#define __TriangleInterface_H

#include "3rd_party_libraries_Macro.hpp"
#include "Eigen/Eigen"
#include "GeometryUtilities.hpp"

#if ENABLE_TRIANGLE == 1
#include "mytriangle.hpp"
#endif

#include "IMeshDAO.hpp"

namespace Gedim
{
class TriangleInterface final
{
  private:
#if ENABLE_TRIANGLE == 1
    void CreateTriangleInput(const Eigen::MatrixXd &polygonVertices,
                             struct triangulateio &triangleInput,
                             const Eigen::MatrixXd &constrainedPoints = Eigen::MatrixXd(),
                             const Eigen::MatrixXi &constrainedSegments = Eigen::MatrixXi()) const;
    void CreateTriangleOutput(const double &maxTriangleArea,
                              struct triangulateio &triangleInput,
                              struct triangulateio &triangleOutput,
                              const std::string &triangleOptions = "-QDzpqnea") const;

    void ConvertTriangleMeshToGedimMesh(struct triangulateio &triangleOutput, IMeshDAO &mesh) const;

    void DeleteTriangleStructure(struct triangulateio &triangleInput, struct triangulateio &triangleOutput) const;

    void ExportTriangleMesh(const struct triangulateio &triangleInput,
                            const struct triangulateio &triangleOutput,
                            const std::string &nameFolder,
                            const std::string &nameFile) const;
#endif
  public:
    TriangleInterface();
    ~TriangleInterface();

    void ExportMesh(const Eigen::MatrixXd &polygonVertices,
                    const double &maxTriangleArea,
                    const std::string &folder_name,
                    const std::string &file_name,
                    const std::string &triangleOptions = "-QDzpqnea") const;

    void CreateMesh(const Eigen::MatrixXd &polygonVertices,
                    const double &maxTriangleArea,
                    Gedim::IMeshDAO &mesh,
                    const std::string &triangleOptions = "-QDzpqnea") const;
};

namespace External
{
namespace MeshUtilities
{
/// \brief Create triangular mesh on 2D polygon
/// \param polygonVertices the 2D polygon vertices, size 3xnumVertices
/// \param maxTriangleArea the maximum triangular area
/// \param options mesh options, see https://www.cs.cmu.edu/~quake/triangle.switch.html
/// \note markers on border are set as { 1, 2, 3, 4, ..., numVertices } for cell0Ds and { 5, 6, 7, 8, ..., 2 *
/// numVertices } for cell1Ds \note use triangle library
void CreateTriangularMesh(const Eigen::MatrixXd &polygonVertices,
                          const double &maxTriangleArea,
                          Gedim::IMeshDAO &mesh,
                          const std::string &options = "-QDzpqnea");

void CreateQuadrilateralMeshFromTriangularMesh(const Gedim::GeometryUtilities &geometryUtilities,
                                               const Eigen::MatrixXd &polygonVertices,
                                               const double &maxTriangleArea,
                                               Gedim::IMeshDAO &mesh);

} // namespace MeshUtilities
} // namespace External

} // namespace Gedim

#endif // __TriangleInterface_H
