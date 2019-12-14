// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the eWoms project.

  eWoms is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  eWoms is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with eWoms.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \brief A test for numerical integration using the vertex-centered finite
 * volume geometries.
 */
#include "config.h"

#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/common/mcmgmapper.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dune/common/version.hh>
#include <dune/geometry/quadraturerules.hh>

#include <ewoms/numerics/discretizations/vcfv/vcfvstencil.hh>

#include <ewoms/common/unused.hh>

#if HAVE_DUNE_ALUGRID
#define EWOMS_NO_ALUGRID_UNUSED
#else
#define EWOMS_NO_ALUGRID_UNUSED  EWOMS_UNUSED
#endif

const unsigned dim = 3;
typedef double Scalar;
typedef Ewoms::QuadrialteralQuadratureGeometry<Scalar, dim> QuadratureGeom;
typedef QuadratureGeom::LocalPosition LocalPosition;
typedef QuadratureGeom::GlobalPosition GlobalPosition;

// function prototypes
GlobalPosition::field_type f(const GlobalPosition &pos);
void testIdenityMapping();
template <class Grid>
void writeTetrahedronSubControlVolumes(const Grid &grid);
void testTetrahedron();
template <class Grid>
void writeCubeSubControlVolumes(const Grid &grid);
void testCube();
void testQuadrature();

GlobalPosition::field_type f(const GlobalPosition &pos)
{
    GlobalPosition::field_type result = 1;
    for (unsigned i = 0; i < GlobalPosition::dimension; ++i)
        result *= pos[i];
    return result;
}

void testIdenityMapping()
{
    QuadratureGeom foo;

    Scalar corners[][3] = { { 0, 0, 0 },
                            { 1, 0, 0 },
                            { 0, 1, 0 },
                            { 1, 1, 0 },
                            { 0, 0, 1 },
                            { 1, 0, 1 },
                            { 0, 1, 1 },
                            { 1, 1, 1 } };
    foo.setCorners(corners, 8);

    std::cout << "testing identity mapping...\n";
    unsigned n = 100;
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            for (unsigned k = 0; k < n; ++k) {
                LocalPosition localPos;

                localPos[0] = Scalar(i) / (n - 1);
                localPos[1] = Scalar(j) / (n - 1);
                localPos[2] = Scalar(k) / (n - 1);

                GlobalPosition globalPos = foo.global(localPos);

                GlobalPosition diff(localPos);
                diff -= globalPos;
                assert(diff.two_norm() < 1e-10);
            }
        }
    }
}

template <class Grid>
void writeTetrahedronSubControlVolumes(const Grid& EWOMS_NO_ALUGRID_UNUSED grid)
{
#if HAVE_DUNE_ALUGRID
    typedef typename Grid::LeafGridView GridView;

    typedef Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming> Grid2;
    typedef typename Grid2::LeafGridView GridView2;
    typedef Dune::GridFactory<Grid2> GridFactory2;

    GridFactory2 gf2;
    const auto &gridView = grid.leafView();
    typedef Ewoms::VcfvStencil<Scalar, GridView> Stencil;
    typedef typename Stencil :: Mapper Mapper;

#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
    Mapper mapper(gridView, Dune::mcmgVertexLayout());
#else
    Mapper mapper(gridView);
#endif
    Stencil stencil(gridView, mapper);

    auto eIt = gridView.template begin<0>();
    const auto &eEndIt = gridView.template end<0>();
    for (; eIt != eEndIt; ++eIt) {
        stencil.update(*eIt);
        for (unsigned scvIdx = 0; scvIdx < stencil.numDof(); ++scvIdx) {
            const auto &scvLocalGeom = stencil.subControlVolume(scvIdx).localGeometry();

            for (unsigned i = 0; i < scvLocalGeom.numCorners; ++i) {
                GlobalPosition pos(
                    eIt->geometry().global(scvLocalGeom.corner(i)));
                gf2.insertVertex(pos);
            }
        }
    }

    unsigned cornerOffset = 0;
    eIt = gridView.template begin<0>();
    for (; eIt != eEndIt; ++eIt) {
        stencil.update(*eIt);
        for (unsigned scvIdx = 0; scvIdx < stencil.numDof(); ++scvIdx) {
            const auto &scvLocalGeom = stencil.subControlVolume(scvIdx).localGeometry();

            std::vector<unsigned> vertexIndices;
            for (unsigned i = 0; i < scvLocalGeom.numCorners; ++i) {
                vertexIndices.push_back(cornerOffset);
                ++cornerOffset;
            }

            gf2.insertElement(Dune::GeometryType(/*topologyId=*/(1 << dim) - 1, dim),
                              vertexIndices);
        }
    }

    std::unique_ptr<Grid2> grid2(gf2.createGrid());
    typedef Dune::VTKWriter<GridView2> VtkWriter;
    VtkWriter writer(grid2->leafView(), Dune::VTK::conforming);
    writer.write("tetrahedron-scvs", Dune::VTK::ascii);
#endif // HAVE_DUNE_ALUGRID
}

void testTetrahedron()
{
#if HAVE_DUNE_ALUGRID
    typedef Dune::ALUGrid<dim, dim, Dune::simplex, Dune::nonconforming> Grid;
    typedef Dune::GridFactory<Grid> GridFactory;
    GridFactory gf;
    Scalar corners[][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

    for (unsigned i = 0; i < sizeof(corners) / sizeof(corners[0]); ++i) {
        GlobalPosition pos;
        for (unsigned j = 0; j < dim; ++j)
            pos[j] = corners[i][j];
        gf.insertVertex(pos);
    }
    std::vector<unsigned int> v = { 0, 1, 2, 3 };
    // in Dune >= 2.6 topologyIds seem to be opaque integers. WTF!?
    gf.insertElement(Dune::GeometryType(/*topologyId=*/0, dim), v);
    std::unique_ptr<Grid> gridPtr(gf.createGrid());

    // write the sub-control volumes to a VTK file.
    writeTetrahedronSubControlVolumes(*gridPtr);
#endif // HAVE_DUNE_ALUGRID
}

template <class Grid>
void writeCubeSubControlVolumes(const Grid& EWOMS_NO_ALUGRID_UNUSED grid)
{
#if HAVE_DUNE_ALUGRID
    typedef typename Grid::LeafGridView GridView;

    typedef Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming> Grid2;
    typedef typename Grid2::LeafGridView GridView2;
    typedef Dune::GridFactory<Grid2> GridFactory2;
    typedef Ewoms::VcfvStencil<Scalar, GridView> Stencil;

    GridFactory2 gf2;
    const auto &gridView = grid.leafView();

#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView> VertexMapper;
    VertexMapper vertexMapper(gridView, Dune::mcmgVertexLayout());
#else
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGVertexLayout> VertexMapper;
    VertexMapper vertexMapper(gridView);
#endif
    Stencil stencil(gridView, vertexMapper);
    auto eIt = gridView.template begin<0>();
    const auto &eEndIt = gridView.template end<0>();
    for (; eIt != eEndIt; ++eIt) {
        stencil.update(*eIt);
        for (unsigned scvIdx = 0; scvIdx < stencil.numDof(); ++scvIdx) {
            const auto &scvLocalGeom = stencil.subControlVolume(scvIdx).localGeometry();

            for (unsigned i = 0; i < scvLocalGeom.numCorners; ++i) {
                GlobalPosition pos(
                    eIt->geometry().global(scvLocalGeom.corner(i)));
                gf2.insertVertex(pos);
            }
        }
    }

    unsigned cornerOffset = 0;
    eIt = gridView.template begin<0>();
    for (; eIt != eEndIt; ++eIt) {
        stencil.update(*eIt);
        for (unsigned scvIdx = 0; scvIdx < stencil.numDof(); ++scvIdx) {
            const auto &scvLocalGeom = stencil.subControlVolume(scvIdx).localGeometry();

            std::vector<unsigned int> vertexIndices;
            for (unsigned i = 0; i < scvLocalGeom.numCorners; ++i) {
                vertexIndices.push_back(cornerOffset);
                ++cornerOffset;
            }

            gf2.insertElement(Dune::GeometryType(/*topologyId=*/(1 << dim) - 1, dim),
                              vertexIndices);
        }
    }

    std::unique_ptr<Grid2> grid2(gf2.createGrid());
    typedef Dune::VTKWriter<GridView2> VtkWriter;
    VtkWriter writer(grid2->leafView(), Dune::VTK::conforming);
    writer.write("cube-scvs", Dune::VTK::ascii);
#endif // HAVE_DUNE_ALUGRID
}

void testCube()
{
#if HAVE_DUNE_ALUGRID
    typedef Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming> Grid;
    typedef Dune::GridFactory<Grid> GridFactory;
    GridFactory gf;
    Scalar corners[][3] = { { 0, 0, 0 },
                            { 1, 0, 0 },
                            { 0, 2, 0 },
                            { 3, 3, 0 },
                            { 0, 0, 4 },
                            { 5, 0, 5 },
                            { 0, 6, 6 },
                            { 7, 7, 7 }, };

    for (unsigned i = 0; i < sizeof(corners) / sizeof(corners[0]); ++i) {
        GlobalPosition pos;
        for (unsigned j = 0; j < dim; ++j)
            pos[j] = corners[i][j];
        gf.insertVertex(pos);
    }
    std::vector<unsigned int> v = { 0, 1, 2, 3, 4, 5, 6, 7 };
    // in Dune >= 2.6 topologyIds seem to be opaque integers. WTF!?
    gf.insertElement(Dune::GeometryType((1 << dim) - 1, dim), v);
    std::unique_ptr<Grid> gridPtr(gf.createGrid());

    // write the sub-control volumes to a VTK file.
    writeCubeSubControlVolumes(*gridPtr);
#endif // HAVE_DUNE_ALUGRID
}

void testQuadrature()
{
    std::cout << "testing SCV quadrature...\n";

    std::array<int, dim> cellRes;

    std::fill(cellRes.begin(), cellRes.end(), 10);

    GlobalPosition upperRight(1.0);

    typedef Dune::YaspGrid<dim> Grid;
    typedef Grid::LeafGridView GridView;
    Grid grid(upperRight, cellRes);

    // compute approximate integral
    auto gridView = grid.leafGridView();
    auto eIt = gridView.begin<0>();
    const auto eEndIt = gridView.end<0>();
    Scalar result = 0;
    // instanciate a stencil
    typedef Ewoms::VcfvStencil<Scalar, GridView> Stencil;
    typedef typename Stencil :: Mapper Mapper;

#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
    Mapper mapper(gridView, Dune::mcmgVertexLayout());
#else
    Mapper mapper(gridView);
#endif
    Stencil stencil(gridView, mapper);
    for (; eIt != eEndIt; ++eIt) {
        const auto &elemGeom = eIt->geometry();

        stencil.update(*eIt);

        // loop over all sub-control volumes
        for (unsigned scvIdx = 0; scvIdx < stencil.numDof(); ++scvIdx) {
            const auto &scvLocalGeom = stencil.subControlVolume(scvIdx).localGeometry();

            Dune::GeometryType geomType = scvLocalGeom.type();
            static const unsigned quadratureOrder = 2;
            const auto &rule
                = Dune::QuadratureRules<Scalar, dim>::rule(geomType,
                                                           quadratureOrder);

            // integrate f over the sub-control volume
            for (auto it = rule.begin(); it != rule.end(); ++it) {
                auto posScvLocal = it->position();
                auto posElemLocal = scvLocalGeom.global(posScvLocal);
                auto posGlobal = elemGeom.global(posScvLocal);

                Scalar fval = f(posGlobal);
                Scalar weight = it->weight();
                Scalar detjac = scvLocalGeom.integrationElement(posScvLocal)
                                * elemGeom.integrationElement(posElemLocal);

                result += fval * weight * detjac;
            }
        }
    }

    std::cout << "result: " << result
              << " (expected value: " << 1.0 / (1 << dim) << ")\n";
}

int main(int argc, char **argv)
{
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);

    testIdenityMapping();
// test the quadrature in a tetrahedron. since the CLang compiler
// prior to version 3.2 generates incorrect code here, we do not
// do it if the compiler is clang 3.1 or older.
#if !__clang__ || __clang_major__ > 3                                          \
    || (__clang_major__ == 3 && __clang_minor__ >= 3)
    testTetrahedron();
    testCube();
#endif
    testQuadrature();

    return 0;
}
