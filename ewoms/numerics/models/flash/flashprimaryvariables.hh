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
 * \copydoc Ewoms::FlashPrimaryVariables
 */
#ifndef EWOMS_FLASH_PRIMARY_VARIABLES_HH
#define EWOMS_FLASH_PRIMARY_VARIABLES_HH

#include "flashindices.hh"
#include "flashproperties.hh"

#include <ewoms/numerics/discretizations/common/fvbaseprimaryvariables.hh>
#include <ewoms/numerics/common/energymodule.hh>

#include <ewoms/material/constraintsolvers/ncpflash.hh>
#include <ewoms/material/fluidstates/compositionalfluidstate.hh>
#include <ewoms/common/valgrind.hh>
#include <ewoms/common/unused.hh>

#include <dune/common/fvector.hh>

#include <iostream>

namespace Ewoms {

/*!
 * \ingroup FlashModel
 *
 * \brief Represents the primary variables used by the compositional
 *        flow model based on flash calculations.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
template <class TypeTag>
class FlashPrimaryVariables : public FvBasePrimaryVariables<TypeTag>
{
    using ParentType = FvBasePrimaryVariables<TypeTag>;

    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using MaterialLawParams = GET_PROP_TYPE(TypeTag, MaterialLawParams);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);

    // primary variable indices
    enum { cTot0Idx = Indices::cTot0Idx };

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    using EnergyModule = Ewoms::EnergyModule<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy)>;

public:
    FlashPrimaryVariables() : ParentType()
    { Ewoms::Valgrind::SetDefined(*this); }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
     */
    FlashPrimaryVariables(Scalar value) : ParentType(value)
    {
        Ewoms::Valgrind::CheckDefined(value);
        Ewoms::Valgrind::SetDefined(*this);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const
     * ImmisciblePrimaryVariables& )
     */
    FlashPrimaryVariables(const FlashPrimaryVariables& value) = default;
    FlashPrimaryVariables& operator=(const FlashPrimaryVariables& value) = default;

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignMassConservative
     */
    template <class FluidState>
    void assignMassConservative(const FluidState& fluidState,
                                const MaterialLawParams& matParams EWOMS_UNUSED,
                                bool isInEquilibrium EWOMS_UNUSED= false)
    {
        // there is no difference between naive and mass conservative
        // assignment in the flash model. (we only need the total
        // concentrations.)
        assignNaive(fluidState);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignNaive
     */
    template <class FluidState>
    void assignNaive(const FluidState& fluidState)
    {
        // reset everything
        (*this) = 0.0;

        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(*this, fluidState);

        // determine the phase presence.
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                this->operator[](cTot0Idx + compIdx) +=
                    fluidState.molarity(phaseIdx, compIdx) * fluidState.saturation(phaseIdx);
            }
        }
    }

    /*!
     * \brief Prints the names of the primary variables and their values.
     *
     * \param os The \c std::ostream which should be used for the output.
     */
    void print(std::ostream& os = std::cout) const
    {
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            os << "(c_tot," << FluidSystem::componentName(compIdx) << " = "
               << this->operator[](cTot0Idx + compIdx);
        }
        os << ")" << std::flush;
    }
};

} // namespace Ewoms

#endif
