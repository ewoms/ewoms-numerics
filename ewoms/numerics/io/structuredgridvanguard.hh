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
 * \copydoc Ewoms::StructuredGridVanguard
 */
#ifndef EWOMS_STRUCTURED_GRID_VANGUARD_HH
#define EWOMS_STRUCTURED_GRID_VANGUARD_HH

#include <ewoms/numerics/io/basevanguard.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#endif

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <vector>
#include <memory>

namespace Ewoms {

template <class TypeTag>
class StructuredGridVanguard;

} // namespace Ewoms

BEGIN_PROPERTIES

NEW_TYPE_TAG(StructuredGridVanguard);

// declare the properties required by the for the structured grid simulator vanguard
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(Scalar);

NEW_PROP_TAG(DomainSizeX);
NEW_PROP_TAG(DomainSizeY);
NEW_PROP_TAG(DomainSizeZ);

NEW_PROP_TAG(CellsX);
NEW_PROP_TAG(CellsY);
NEW_PROP_TAG(CellsZ);

NEW_PROP_TAG(GridGlobalRefinements);

// GRIDDIM is only set by the finger problem
#ifndef GRIDDIM
static const int dim = 2;
#else
static const int dim = GRIDDIM;
#endif

// set the Grid and Vanguard properties
#if HAVE_DUNE_ALUGRID
SET_TYPE_PROP(StructuredGridVanguard, Grid, Dune::ALUGrid< dim, dim, Dune::cube, Dune::nonconforming >);
#else
SET_TYPE_PROP(StructuredGridVanguard, Grid, Dune::YaspGrid< dim >);
#endif

SET_TYPE_PROP(StructuredGridVanguard, Vanguard, Ewoms::StructuredGridVanguard<TypeTag>);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup TestProblems
 *
 * \brief Helper class for grid instantiation of the lens problem.
 */
template <class TypeTag>
class StructuredGridVanguard : public BaseVanguard<TypeTag>
{
    using ParentType = BaseVanguard<TypeTag>;
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using Grid = GET_PROP_TYPE(TypeTag, Grid);

    using GridPointer = std::unique_ptr<Grid>;

    static const int dim = Grid::dimension;

public:
    /*!
     * \brief Register all run-time parameters for the structured grid simulator vanguard.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, GridGlobalRefinements,
                             "The number of global refinements of the grid "
                             "executed after it was loaded");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, DomainSizeX,
                             "The size of the domain in x direction");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, CellsX,
                             "The number of intervalls in x direction");
        if (dim > 1) {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, DomainSizeY,
                                 "The size of the domain in y direction");
            EWOMS_REGISTER_PARAM(TypeTag, unsigned, CellsY,
                                 "The number of intervalls in y direction");
        }
        if (dim > 2) {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, DomainSizeZ,
                                 "The size of the domain in z direction");
            EWOMS_REGISTER_PARAM(TypeTag, unsigned, CellsZ,
                                 "The number of intervalls in z direction");
        }
    }

    /*!
     * \brief Create the grid for the lens problem
     */
    StructuredGridVanguard(Simulator& simulator)
        : ParentType(simulator)
    {
        Dune::FieldVector<int, dim> cellRes;

        using GridScalar = double;
        Dune::FieldVector<GridScalar, dim> upperRight;
        Dune::FieldVector<GridScalar, dim> lowerLeft( 0 );

        upperRight[0] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeX);
        upperRight[1] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeY);

        cellRes[0] = EWOMS_GET_PARAM(TypeTag, unsigned, CellsX);
        cellRes[1] = EWOMS_GET_PARAM(TypeTag, unsigned, CellsY);
        if (dim == 3) {
            upperRight[2] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeZ);
            cellRes[2] = EWOMS_GET_PARAM(TypeTag, unsigned, CellsZ);
        }

        std::stringstream dgffile;
        dgffile << "DGF" << std::endl;
        dgffile << "INTERVAL" << std::endl;
        dgffile << lowerLeft  << std::endl;
        dgffile << upperRight << std::endl;
        dgffile << cellRes    << std::endl;
        dgffile << "#" << std::endl;
        dgffile << "GridParameter" << std::endl;
        dgffile << "overlap 1" << std::endl;
        dgffile << "#" << std::endl;
        dgffile << "Simplex" << std::endl;
        dgffile << "#" << std::endl;

        // use DGF parser to create a grid from interval block
        gridPtr_.reset(Dune::GridPtr<Grid>(dgffile).release());

        unsigned numRefinements = EWOMS_GET_PARAM(TypeTag, unsigned, GridGlobalRefinements);
        gridPtr_->globalRefine(static_cast<int>(numRefinements));

        this->finalizeInit_();
    }

    /*!
     * \brief Return a reference to the grid object.
     */
    Grid& grid()
    { return *gridPtr_; }

    /*!
     * \brief Return a constant reference to the grid object.
     */
    const Grid& grid() const
    { return *gridPtr_; }

private:
    GridPointer gridPtr_;
};

} // namespace Ewoms

#endif
