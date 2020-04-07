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
 *
 * \copydoc Ewoms::NcpPrimaryVariables
 */
#ifndef EWOMS_NCP_PRIMARY_VARIABLES_HH
#define EWOMS_NCP_PRIMARY_VARIABLES_HH

#include "ncpproperties.hh"

#include <ewoms/numerics/discretizations/common/fvbaseprimaryvariables.hh>
#include <ewoms/numerics/common/energymodule.hh>

#include <ewoms/material/constraintsolvers/ncpflash.hh>
#include <ewoms/material/fluidstates/compositionalfluidstate.hh>
#include <ewoms/common/densead/math.hh>

#include <dune/common/fvector.hh>

namespace Ewoms {

/*!
 * \ingroup NcpModel
 *
 * \brief Represents the primary variables used by the compositional
 *        multi-phase NCP model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
template <class TypeTag>
class NcpPrimaryVariables : public FvBasePrimaryVariables<TypeTag>
{
    typedef FvBasePrimaryVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };
    enum { fugacity0Idx = Indices::fugacity0Idx };

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef Ewoms::EnergyModule<TypeTag, enableEnergy> EnergyModule;

    typedef Ewoms::NcpFlash<Scalar, FluidSystem> NcpFlash;

public:
    NcpPrimaryVariables() : ParentType()
    {}

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
     */
    NcpPrimaryVariables(Scalar value) : ParentType(value)
    {}

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const
     * ImmisciblePrimaryVariables& )
     */
    NcpPrimaryVariables(const NcpPrimaryVariables& value) = default;
    NcpPrimaryVariables& operator=(const NcpPrimaryVariables& value) = default;

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignMassConservative
     */
    template <class FluidState>
    void assignMassConservative(const FluidState& fluidState,
                                const MaterialLawParams& matParams,
                                bool isInEquilibrium = false)
    {
        typedef Ewoms::MathToolbox<typename FluidState::Scalar> FsToolbox;

#ifndef NDEBUG
        // make sure the temperature is the same in all fluid phases
        for (unsigned phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx) {
            assert(fluidState.temperature(0) == fluidState.temperature(phaseIdx));
        }
#endif // NDEBUG

        // for the equilibrium case, we don't need complicated
        // computations.
        if (isInEquilibrium) {
            assignNaive(fluidState);
            return;
        }

        // use a flash calculation to calculate a fluid state in
        // thermodynamic equilibrium
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        Ewoms::CompositionalFluidState<Scalar, FluidSystem> fsFlash;

        // use the externally given fluid state as initial value for
        // the flash calculation
        fsFlash.assign(fluidState);

        // calculate the phase densities
        paramCache.updateAll(fsFlash);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar rho = FluidSystem::density(fsFlash, paramCache, phaseIdx);
            fsFlash.setDensity(phaseIdx, rho);
        }

        // calculate the "global molarities"
        ComponentVector globalMolarities(0.0);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                globalMolarities[compIdx] +=
                    FsToolbox::value(fsFlash.saturation(phaseIdx))
                    * FsToolbox::value(fsFlash.molarity(phaseIdx, compIdx));
            }
        }

        // run the flash calculation
        NcpFlash::template solve<MaterialLaw>(fsFlash, matParams, paramCache, globalMolarities);

        // use the result to assign the primary variables
        assignNaive(fsFlash);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignNaive
     */
    template <class FluidState>
    void assignNaive(const FluidState& fluidState, unsigned refPhaseIdx = 0)
    {
        typedef Ewoms::MathToolbox<typename FluidState::Scalar> FsToolbox;

        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(*this, fluidState);

        // assign fugacities.
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updatePhase(fluidState, refPhaseIdx);
        Scalar pRef = FsToolbox::value(fluidState.pressure(refPhaseIdx));
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            // we always compute the fugacities because they are quite exotic quantities
            // and this easily forgotten to be specified
            Scalar fugCoeff =
                FluidSystem::template fugacityCoefficient<FluidState, Scalar>(fluidState,
                                                                              paramCache,
                                                                              refPhaseIdx,
                                                                              compIdx);
            (*this)[fugacity0Idx + compIdx] =
                fugCoeff*fluidState.moleFraction(refPhaseIdx, compIdx)*pRef;
        }

        // assign pressure of first phase
        (*this)[pressure0Idx] = FsToolbox::value(fluidState.pressure(/*phaseIdx=*/0));

        // assign first M - 1 saturations
        for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            (*this)[saturation0Idx + phaseIdx] = FsToolbox::value(fluidState.saturation(phaseIdx));
    }
};

} // namespace Ewoms

#endif
