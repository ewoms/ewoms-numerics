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
 *
 * \copydoc Ewoms::DiscreteFracturePrimaryVariables
 */
#ifndef EWOMS_DISCRETE_FRACTURE_PRIMARY_VARIABLES_HH
#define EWOMS_DISCRETE_FRACTURE_PRIMARY_VARIABLES_HH

#include "discretefractureproperties.hh"

#include <ewoms/numerics/models/immiscible/immiscibleprimaryvariables.hh>

namespace Ewoms {
/*!
 * \ingroup DiscreteFractureModel
 *
 * \brief Represents the primary variables used by the discrete fracture
 *        multi-phase model.
 */
template <class TypeTag>
class DiscreteFracturePrimaryVariables
    : public ImmisciblePrimaryVariables<TypeTag>
{
    using ParentType = ImmisciblePrimaryVariables<TypeTag>;

    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using MaterialLaw = GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = GET_PROP_TYPE(TypeTag, MaterialLawParams);

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

public:
    /*!
     * \brief Default constructor
     */
    DiscreteFracturePrimaryVariables() : ParentType()
    {}

    /*!
     * \brief Constructor with assignment from scalar
     *
     * \param value The scalar value to which all entries of the vector will be set.
     */
    DiscreteFracturePrimaryVariables(Scalar value) : ParentType(value)
    {}

    /*!
     * \brief Copy constructor
     *
     * \param value The primary variables that will be duplicated.
     */
    DiscreteFracturePrimaryVariables(const DiscreteFracturePrimaryVariables& value) = default;
    DiscreteFracturePrimaryVariables& operator=(const DiscreteFracturePrimaryVariables& value) = default;

    /*!
     * \brief Directly retrieve the primary variables from an
     *        arbitrary fluid state of the fractures.
     *
     * \param fractureFluidState The fluid state of the fractures
     *                           which should be represented by the
     *                           primary variables. The temperatures,
     *                           pressures and compositions of all
     *                           phases must be defined.
     * \param matParams The parameters for the capillary-pressure law
     *                  which apply for the fracture.
     */
    template <class FluidState>
    void assignNaiveFromFracture(const FluidState& fractureFluidState,
                                 const MaterialLawParams& matParams)
    {
        FluidState matrixFluidState;
        fractureToMatrixFluidState_(matrixFluidState, fractureFluidState,
                                    matParams);

        ParentType::assignNaive(matrixFluidState);
    }

private:
    template <class FluidState>
    void fractureToMatrixFluidState_(FluidState& matrixFluidState,
                                     const FluidState& fractureFluidState,
                                     const MaterialLawParams& matParams) const
    {
        // start with the same fluid state as in the fracture
        matrixFluidState.assign(fractureFluidState);

        // the condition for the equilibrium is that the pressures are
        // the same in the fracture and in the matrix. This means that
        // we have to find saturations for the matrix which result in
        // the same pressures as in the fracture. this can be done by
        // inverting the capillary pressure-saturation curve.
        Scalar saturations[numPhases];
        MaterialLaw::saturations(saturations, matParams, matrixFluidState);

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            matrixFluidState.setSaturation(phaseIdx, saturations[phaseIdx]);
    }
};

} // namespace Ewoms

#endif
