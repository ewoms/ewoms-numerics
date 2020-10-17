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
*/
/*!
 * \file
 * \copydoc Ewoms::SimplexGridVanguard
 */
#ifndef EWOMS_SIMPLEX_GRID_VANGUARD_HH
#define EWOMS_SIMPLEX_GRID_VANGUARD_HH

#include <ewoms/numerics/utils/basicproperties.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/common/fvector.hh>

#include <memory>

BEGIN_PROPERTIES

NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Grid);

NEW_PROP_TAG(DomainSizeX);
NEW_PROP_TAG(DomainSizeY);
NEW_PROP_TAG(DomainSizeZ);

NEW_PROP_TAG(CellsX);
NEW_PROP_TAG(CellsY);
NEW_PROP_TAG(CellsZ);

NEW_PROP_TAG(GridGlobalRefinements);

END_PROPERTIES

namespace Ewoms {
/*!
 * \brief Provides a simulator vanguard which a creates regular grid made of simplices.
 */
template <class TypeTag>
class SimplexGridVanguard
{
    using ParentType = BaseVanguard<TypeTag>;
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using Grid = GET_PROP_TYPE(TypeTag, Grid);

    using GridPointer = std::unique_ptr<Grid>;
    using CoordScalar = typename Grid::ctype;
    enum { dimWorld = Grid::dimensionworld };
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

public:
    /*!
     * \brief Register all run-time parameters for the grid manager.
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
        if (dimWorld > 1) {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, DomainSizeY,
                                 "The size of the domain in y direction");
            EWOMS_REGISTER_PARAM(TypeTag, unsigned, CellsY,
                                 "The number of intervalls in y direction");
        }
        if (dimWorld > 2) {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, DomainSizeZ,
                                 "The size of the domain in z direction");
            EWOMS_REGISTER_PARAM(TypeTag, unsigned, CellsZ,
                                 "The number of intervalls in z direction");
        }
    }

    /*!
     * \brief Create the Grid
     */
    SimplexGridVanguard(Simulator& simulator)
        : ParentType(simulator)
    {
        Dune::array<unsigned, dimWorld> cellRes;
        GlobalPosition upperRight;
        GlobalPosition lowerLeft;

        lowerLeft[0] = 0.0;
        upperRight[0] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeX);
        cellRes[0] = EWOMS_GET_PARAM(TypeTag, unsigned, CellsX);
        if (dimWorld > 1) {
            lowerLeft[1] = 0.0;
            upperRight[1] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeY);
            cellRes[1] = EWOMS_GET_PARAM(TypeTag, unsigned, CellsY);
        }
        if (dimWorld > 2) {
            lowerLeft[2] = 0.0;
            upperRight[2] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeZ);
            cellRes[2] = EWOMS_GET_PARAM(TypeTag, unsigned, CellsZ);
        }

        simplexGrid_ = Dune::StructuredGridFactory<Grid>::createSimplexGrid(lowerLeft,
                                                                            upperRight,
                                                                            cellRes);

        unsigned numRefinments = EWOMS_GET_PARAM(TypeTag, unsigned, GridGlobalRefinements);
        simplexGrid_->globalRefine(numRefinments);

        this->finalizeInit_();
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    Grid& grid()
    { return simplexGrid_; }

    /*!
     * \brief Returns a reference to the grid.
     */
    const Grid& grid() const
    { return *simplexGrid_; }

private:
    GridPointer simplexGrid_;
};
} // namespace Ewoms

#endif
