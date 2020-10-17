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
 * \copydoc Ewoms::RichardsPrimaryVariables
 */
#ifndef EWOMS_RICHARDS_PRIMARY_VARIABLES_HH
#define EWOMS_RICHARDS_PRIMARY_VARIABLES_HH

#include "richardsproperties.hh"

#include <ewoms/numerics/discretizations/common/fvbaseprimaryvariables.hh>

#include <ewoms/material/constraintsolvers/immiscibleflash.hh>
#include <ewoms/material/fluidstates/immisciblefluidstate.hh>
#include <ewoms/common/valgrind.hh>
#include <ewoms/common/unused.hh>

#include <dune/common/fvector.hh>

namespace Ewoms {

/*!
 * \ingroup RichardsModel
 *
 * \brief Represents the primary variables used in the Richards model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
template <class TypeTag>
class RichardsPrimaryVariables : public FvBasePrimaryVariables<TypeTag>
{
    using ParentType = FvBasePrimaryVariables<TypeTag>;
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
    using MaterialLaw = GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = GET_PROP_TYPE(TypeTag, MaterialLawParams);
    using EnergyModule = GET_PROP_TYPE(TypeTag, IntensiveQuantities);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);

    // primary variable indices
    enum { pressureWIdx = Indices::pressureWIdx };

    enum { liquidPhaseIdx = GET_PROP_VALUE(TypeTag, LiquidPhaseIndex) };
    enum { gasPhaseIdx = GET_PROP_VALUE(TypeTag, GasPhaseIndex) };

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;
    using PhaseVector = Dune::FieldVector<Scalar, numPhases>;
    using ImmiscibleFlash = Ewoms::ImmiscibleFlash<Scalar, FluidSystem>;

public:
    RichardsPrimaryVariables() : ParentType()
    { Ewoms::Valgrind::SetUndefined(*this); }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
     */
    RichardsPrimaryVariables(Scalar value) : ParentType(value)
    {}

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const
     * ImmisciblePrimaryVariables& )
     */
    RichardsPrimaryVariables(const RichardsPrimaryVariables& value) = default;
    RichardsPrimaryVariables& operator=(const RichardsPrimaryVariables& value) = default;

    /*!
     * \brief Set the primary variables with the wetting phase
     *        pressure, saturation and temperature.
     *
     * \param T The temperature [K]
     * \param pw The pressure of the wetting phase [Pa]
     * \param Sw The saturation of the wetting phase []
     * \param matParams The capillary pressure law parameters
     */
    void assignImmiscibleFromWetting(Scalar T, Scalar pw, Scalar Sw,
                                     const MaterialLawParams& matParams)
    {
        Ewoms::ImmiscibleFluidState<Scalar, FluidSystem> fs;

        fs.setTemperature(T);
        fs.setSaturation(liquidPhaseIdx, Sw);
        fs.setSaturation(gasPhaseIdx, 1 - Sw);

        // set phase pressures
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(liquidPhaseIdx, pw);
        fs.setPressure(gasPhaseIdx, pw + (pC[gasPhaseIdx] - pC[liquidPhaseIdx]));

        assignNaive(fs);
    }

    /*!
     * \brief Set the primary variables with the non-wetting phase
     *        pressure, saturation and temperature.
     *
     * \param T The temperature [K]
     * \param pn The pressure of the non-wetting phase [Pa]
     * \param Sn The saturation of the non-wetting phase []
     * \param matParams The capillary pressure law parameters
     */
    void assignImmiscibleFromNonWetting(Scalar T, Scalar pn, Scalar Sn,
                                        const MaterialLawParams& matParams)
    {
        Ewoms::ImmiscibleFluidState<Scalar, FluidSystem> fs;

        fs.setTemperature(T);
        fs.setSaturation(liquidPhaseIdx, 1 - Sn);
        fs.setSaturation(gasPhaseIdx, Sn);

        // set phase pressures
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(gasPhaseIdx, pn);
        fs.setPressure(gasPhaseIdx, pn + (pC[liquidPhaseIdx] - pC[gasPhaseIdx]));

        assignNaive(fs);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignMassConservative
     */
    template <class FluidState>
    void assignMassConservative(const FluidState& fluidState,
                                const MaterialLawParams& matParams,
                                bool isInEquilibrium EWOMS_UNUSED= false)
    {
        ComponentVector globalMolarities(0.0);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                globalMolarities[compIdx] +=
                    fluidState.molarity(phaseIdx, compIdx) * fluidState.saturation(phaseIdx);
            }
        }

        Ewoms::ImmiscibleFluidState<Scalar, FluidSystem> fsFlash;
        fsFlash.assign(fluidState);
        typename FluidSystem::ParameterCache paramCache;
        ImmiscibleFlash::template solve<MaterialLaw>(fsFlash, paramCache,
                                                     matParams,
                                                     globalMolarities);

        assignNaive(fsFlash);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignNaive
     */
    template <class FluidState>
    void assignNaive(const FluidState& fluidState)
    {
        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(*this, fluidState);

        (*this)[pressureWIdx] = fluidState.pressure(liquidPhaseIdx);
    }
};

} // namespace Ewoms

#endif
