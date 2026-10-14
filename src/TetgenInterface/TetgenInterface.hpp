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

#ifndef __TetgenInterface_H
#define __TetgenInterface_H

#include "3rd_party_libraries_Macro.hpp"
#include "IMeshDAO.hpp"

#if ENABLE_TETGEN == 1
#include "tetgen.h"
#endif

#include "Eigen/Eigen"

namespace Gedim
{
/// \brief The Tetgen Interface
/// \see http://wias-berlin.de/software/tetgen/files/tetcall.cxx
class TetgenInterface final
{
  public:
    struct Region
    {
        int id;                   // unique id of the region
        Eigen::Vector3d centroid; // internal point of the region
        double max_volume;        // default -1.0
    };

  private:
#if ENABLE_TETGEN == 1
    void CreateTetgenInput(const Eigen::MatrixXd &polyhedronVertices,
                           const Eigen::MatrixXi &polyhedronEdges,
                           const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                           tetgenio &tetgenInput,
                           const Eigen::MatrixXd &constrainedPoints = Eigen::MatrixXd(),
                           const std::vector<Eigen::VectorXi> &constrainedFaces = std::vector<Eigen::VectorXi>()) const;
    void CreateTetgenInput(const Eigen::MatrixXd &points,
                           const std::vector<std::vector<unsigned int>> &facets,
                           const std::vector<TetgenInterface::Region> &regions,
                           const std::vector<Eigen::Vector3d> &holes,
                           tetgenio &tetgenInput) const;
    void CreateDelaunayInput(const Eigen::MatrixXd &points, const std::vector<unsigned int> &points_marker, tetgenio &tetgenInput) const;
    void CreateTetgenOutput(const double &maxTetrahedronArea,
                            tetgenio &tetgenInput,
                            tetgenio &tetgenOutput,
                            const std::string &tetgenOptions = "Qpqfezna") const;
    void CreateTetgenOutput(tetgenio &tetgenInput, tetgenio &tetgenOutput, const std::string &tetgenOptions) const;
    void ConvertTetgenOutputToMeshDAO(const tetgenio &tetgenOutput, IMeshDAO &mesh) const;

    void DeleteTetgenStructure(tetgenio &tetgenInput, tetgenio &tetgenOutput) const;
    void ExportTetgenOutput(const std::string &nameFolder, const std::string &nameFile, tetgenio &tetgenOutput) const;
#endif

  public:
    TetgenInterface();
    ~TetgenInterface();

    void CreateDelaunay(const Eigen::MatrixXd &points, const std::vector<unsigned int> &points_marker, Gedim::IMeshDAO &mesh) const;

    void CreateMesh(const Eigen::MatrixXd &polyhedronVertices,
                    const Eigen::MatrixXi &polyhedronEdges,
                    const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                    const double &maxTetrahedronVolume,
                    Gedim::IMeshDAO &mesh,
                    const std::string &tetgenOptions = "Qpqfezna") const;

    void CreateMesh(const Eigen::MatrixXd &points,
                    const std::vector<std::vector<unsigned int>> &facets,
                    const double &maxTetrahedronVolume,
                    Gedim::IMeshDAO &mesh,
                    const std::string &tetgenOptions = "Qpqfezna") const;

    void CreateMesh(const Eigen::MatrixXd &points,
                    const std::vector<std::vector<unsigned int>> &facets,
                    const std::vector<Region> &regions,
                    Gedim::IMeshDAO &mesh,
                    const std::string &tetgenOptions = "QpqfeznaA") const;

    void CreateMesh(const Eigen::MatrixXd &points,
                    const std::vector<std::vector<unsigned int>> &facets,
                    const std::vector<TetgenInterface::Region> &regions,
                    const std::vector<Eigen::Vector3d> &holes,
                    IMeshDAO &mesh,
                    const std::string &tetgenOptions = "QpqfeznaA") const;
};

namespace External
{
namespace MeshUtilities
{

/// \brief Create tetrahedral mesh on 3D polyhedron
/// \param polyhedronVertices the polyhedron vertices, size 3 x numVertices
/// \param polyhedronEdges the polyhedron edges, size 2 x numEdges
/// \param polyhedronFaces the polyhedron face vertices and edges, size numFaces x 2 x numVertices
/// \param maxTetrahedronVolume the maximum tetrahedron area
/// \param options mesh options, see https://wias-berlin.de/software/tetgen/1.5/doc/manual/manual005.html#cmd-q
/// \note markers on border are set as { 1, 2, 3, 4, ..., numVertices } for cell0Ds and { 5, 6, 7, 8, ..., 2 *
/// numVertices } for cell1Ds \note use tetgen library
void CreateTetrahedralMesh(const Eigen::MatrixXd &polyhedronVertices,
                           const Eigen::MatrixXi &polyhedronEdges,
                           const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                           const double &maxTetrahedronVolume,
                           Gedim::IMeshDAO &mesh,
                           const std::string &options = "Qpqfezna");
void CreateTetrahedralMesh(const Eigen::MatrixXd &points,
                           const std::vector<std::vector<unsigned int>> &facets,
                           const double &maxTetrahedronVolume,
                           Gedim::IMeshDAO &mesh,
                           const std::string &options = "Qpqfezna");

void CreateDelaunayMesh3D(const Eigen::MatrixXd &points, const std::vector<unsigned int> &points_marker, Gedim::IMeshDAO &mesh);

} // namespace MeshUtilities

} // namespace External

} // namespace Gedim

#endif // __TetgenInterface_H
