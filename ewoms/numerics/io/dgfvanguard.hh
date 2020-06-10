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
 * \copydoc Ewoms::DgfVanguard
 */
#ifndef EWOMS_DGF_GRID_VANGUARD_HH
#define EWOMS_DGF_GRID_VANGUARD_HH

#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <ewoms/numerics/models/discretefracture/fracturemapper.hh>

#include <ewoms/numerics/io/basevanguard.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <type_traits>
#include <string>

BEGIN_PROPERTIES

NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridFile);
NEW_PROP_TAG(Vanguard);
NEW_PROP_TAG(GridGlobalRefinements);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Simulator);

END_PROPERTIES

namespace Ewoms {

/*!
 * \brief Provides a simulator vanguard which creates a grid by parsing a Dune Grid
 *        Format (DGF) file.
 */
template <class TypeTag>
class DgfVanguard : public BaseVanguard<TypeTag>
{
    typedef BaseVanguard<TypeTag> ParentType;
    typedef GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef Ewoms::FractureMapper<TypeTag> FractureMapper;

    typedef std::unique_ptr< Grid > GridPointer;

public:
    /*!
     * \brief Register all run-time parameters for the DGF simulator vanguard.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, std::string, GridFile,
                             "The file name of the DGF file to load");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, GridGlobalRefinements,
                             "The number of global refinements of the grid "
                             "executed after it was loaded");
    }

    /*!
     * \brief Load the grid from the file.
     */
    DgfVanguard(Simulator& simulator)
        : ParentType(simulator)
    {
        const std::string dgfFileName = EWOMS_GET_PARAM(TypeTag, std::string, GridFile);
        unsigned numRefinments = EWOMS_GET_PARAM(TypeTag, unsigned, GridGlobalRefinements);

        {
            // create DGF GridPtr from a dgf file
            Dune::GridPtr< Grid > dgfPointer( dgfFileName );

            // this is only implemented for 2d currently
            addFractures_( dgfPointer );

            // store pointer to dune grid
            gridPtr_.reset( dgfPointer.release() );
        }

        if (numRefinments > 0)
            gridPtr_->globalRefine(static_cast<int>(numRefinments));

        this->finalizeInit_();
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    Grid& grid()
    { return *gridPtr_; }

    /*!
     * \brief Returns a reference to the grid.
     */
    const Grid& grid() const
    { return *gridPtr_; }

    /*!
     * \brief Distributes the grid on all processes of a parallel
     *        computation.
     *
     * This grid manager plays nice and also distributes the data of
     * the DGF...
     */
    void loadBalance()
    { gridPtr_->loadBalance(); }

    /*!
     * \brief Returns the fracture mapper
     *
     * The fracture mapper determines the topology of the fractures.
     */
    FractureMapper& fractureMapper()
    { return fractureMapper_; }

    /*!
     * \brief Returns the fracture mapper
     *
     * The fracture mapper determines the topology of the fractures.
     */
    const FractureMapper& fractureMapper() const
    { return fractureMapper_; }

protected:
    void addFractures_(Dune::GridPtr<Grid>& dgfPointer)
    {
        typedef typename Grid::LevelGridView LevelGridView;

        // check if fractures are available (only 2d currently)
        if (dgfPointer.nofParameters(static_cast<int>(Grid::dimension)) == 0)
            return;

        LevelGridView gridView = dgfPointer->levelGridView(/*level=*/0);
        const unsigned edgeCodim = Grid::dimension - 1;

#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<LevelGridView> VertexMapper;
        VertexMapper vertexMapper(gridView, Dune::mcmgVertexLayout());
#else
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<LevelGridView, Dune::MCMGVertexLayout> VertexMapper;
        VertexMapper vertexMapper(gridView);
#endif

        // first create a map of the dune to ART vertex indices
        auto eIt = gridView.template begin</*codim=*/0>();
        const auto eEndIt = gridView.template end</*codim=*/0>();
        for (; eIt != eEndIt; ++eIt) {
            const auto& element = *eIt;
            const auto& refElem =
                Dune::ReferenceElements<Scalar, Grid::dimension>::general(element.type());

            const int edges = refElem.size( edgeCodim );
            for (int edge = 0; edge < edges; ++edge) {
                const int vertices = refElem.size(edge, edgeCodim, Grid::dimension);
                std::vector<unsigned> vertexIndices;
                vertexIndices.reserve(Grid::dimension);
                for (int vx = 0; vx < vertices; ++vx) {
                    // get local vertex number from edge
                    const int localVx = refElem.subEntity(edge, edgeCodim, vx, Grid::dimension);

                    // get vertex
                    const auto vertex = element.template subEntity<Grid::dimension>(localVx);

                    // if vertex has parameter 1 insert as a fracture vertex
                    if (dgfPointer.parameters( vertex )[ 0 ] > 0)
                        vertexIndices.push_back(
                            static_cast<unsigned>(vertexMapper.subIndex(element,
                                                                        static_cast<int>(localVx),
                                                                        Grid::dimension)));
                }
                // if 2 vertices have been found with flag 1 insert a fracture edge
                if (static_cast<int>(vertexIndices.size()) == Grid::dimension)
                    fractureMapper_.addFractureEdge(vertexIndices[0], vertexIndices[1]);
            }
        }
    }

private:
    GridPointer    gridPtr_;
    FractureMapper fractureMapper_;
};

} // namespace Ewoms

#endif
