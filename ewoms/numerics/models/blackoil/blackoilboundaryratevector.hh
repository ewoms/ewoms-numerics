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
 * \copydoc Ewoms::BlackOilBoundaryRateVector
 */
#ifndef EWOMS_BLACK_OIL_BOUNDARY_RATE_VECTOR_HH
#define EWOMS_BLACK_OIL_BOUNDARY_RATE_VECTOR_HH

#include <ewoms/common/valgrind.hh>
#include <ewoms/material/constraintsolvers/ncpflash.hh>

#include "blackoilintensivequantities.hh"
#include "blackoilenergymodules.hh"

namespace Ewoms {

/*!
 * \ingroup BlackOilModel
 *
 * \brief Implements a boundary vector for the fully implicit black-oil model.
 */
template <class TypeTag>
class BlackOilBoundaryRateVector : public GET_PROP_TYPE(TypeTag, RateVector)
{
    using ParentType = GET_PROP_TYPE(TypeTag, RateVector);
    using ExtensiveQuantities = GET_PROP_TYPE(TypeTag, ExtensiveQuantities);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
    using LocalResidual = GET_PROP_TYPE(TypeTag, LocalResidual);
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);
    using RateVector = GET_PROP_TYPE(TypeTag, RateVector);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { enableSolvent = GET_PROP_VALUE(TypeTag, EnableSolvent) };
    enum { enablePolymer = GET_PROP_VALUE(TypeTag, EnablePolymer) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { contiEnergyEqIdx = Indices::contiEnergyEqIdx };

    static constexpr bool blackoilConserveSurfaceVolume = GET_PROP_VALUE(TypeTag, BlackoilConserveSurfaceVolume);

    using EnergyModule = Ewoms::BlackOilEnergyModule<TypeTag, enableEnergy>;

public:
    /*!
     * \brief Default constructor
     */
    BlackOilBoundaryRateVector() : ParentType()
    {}

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(Scalar)
     */
    BlackOilBoundaryRateVector(Scalar value) : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(const ImmiscibleBoundaryRateVector& )
     */
    BlackOilBoundaryRateVector(const BlackOilBoundaryRateVector& value) = default;
    BlackOilBoundaryRateVector& operator=(const BlackOilBoundaryRateVector& value) = default;

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setFreeFlow
     */
    template <class Context, class FluidState>
    void setFreeFlow(const Context& context,
                     unsigned bfIdx,
                     unsigned timeIdx,
                     const FluidState& fluidState)
    {
        ExtensiveQuantities extQuants;
        extQuants.updateBoundary(context, bfIdx, timeIdx, fluidState);
        const auto& insideIntQuants = context.intensiveQuantities(bfIdx, timeIdx);
        unsigned focusDofIdx = context.focusDofIndex();
        unsigned interiorDofIdx = context.interiorScvIndex(bfIdx, timeIdx);

        ////////
        // advective fluxes of all components in all phases
        ////////
        (*this) = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const auto& pBoundary = fluidState.pressure(phaseIdx);
            const Evaluation& pInside = insideIntQuants.fluidState().pressure(phaseIdx);

            RateVector tmp;

            // mass conservation
            if (pBoundary < pInside)
                // outflux
                LocalResidual::template evalPhaseFluxes_<Evaluation>(tmp,
                                                                     phaseIdx,
                                                                     insideIntQuants.pvtRegionIndex(),
                                                                     extQuants,
                                                                     insideIntQuants.fluidState());
            else if (pBoundary > pInside) {
                typedef typename std::conditional<std::is_same<typename FluidState::Scalar, Evaluation>::value,
                                                  Evaluation, Scalar>::type RhsEval;
                // influx
                LocalResidual::template evalPhaseFluxes_<RhsEval>(tmp,
                                                                  phaseIdx,
                                                                  insideIntQuants.pvtRegionIndex(),
                                                                  extQuants,
                                                                  fluidState);
            }

            for (unsigned i = 0; i < tmp.size(); ++i)
                (*this)[i] += tmp[i];

            // energy conservation
            if (enableEnergy) {
                Evaluation density;
                Evaluation specificEnthalpy;
                if (pBoundary > pInside) {
                    if (focusDofIdx == interiorDofIdx) {
                        density = fluidState.density(phaseIdx);
                        specificEnthalpy = fluidState.enthalpy(phaseIdx);
                    }
                    else {
                        density = Ewoms::getValue(fluidState.density(phaseIdx));
                        specificEnthalpy = Ewoms::getValue(fluidState.enthalpy(phaseIdx));
                    }
                }
                else if (focusDofIdx == interiorDofIdx) {
                    density = insideIntQuants.fluidState().density(phaseIdx);
                    specificEnthalpy = insideIntQuants.fluidState().enthalpy(phaseIdx);
                }
                else {
                    density = Ewoms::getValue(insideIntQuants.fluidState().density(phaseIdx));
                    specificEnthalpy = Ewoms::getValue(insideIntQuants.fluidState().enthalpy(phaseIdx));
                }

                Evaluation enthalpyRate = density*extQuants.volumeFlux(phaseIdx)*specificEnthalpy;
                EnergyModule::addToEnthalpyRate(*this, enthalpyRate);
            }
        }

        if (enableSolvent) {
            (*this)[Indices::contiSolventEqIdx] = extQuants.solventVolumeFlux();
            if (blackoilConserveSurfaceVolume)
                (*this)[Indices::contiSolventEqIdx] *= insideIntQuants.solventInverseFormationVolumeFactor();
            else
                (*this)[Indices::contiSolventEqIdx] *= insideIntQuants.solventDensity();

        }

        if (enablePolymer) {
            (*this)[Indices::contiPolymerEqIdx] = extQuants.volumeFlux(FluidSystem::waterPhaseIdx) * insideIntQuants.polymerConcentration();
        }

        // make sure that the right mass conservation quantities are used
        LocalResidual::adaptMassConservationQuantities_(*this, insideIntQuants.pvtRegionIndex());

        // heat conduction
        if (enableEnergy)
            EnergyModule::addToEnthalpyRate(*this, extQuants.energyFlux());

#ifndef NDEBUG
        for (unsigned i = 0; i < numEq; ++i) {
            Ewoms::Valgrind::CheckDefined((*this)[i]);
        }
        Ewoms::Valgrind::CheckDefined(*this);
#endif
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setInFlow
     */
    template <class Context, class FluidState>
    void setInFlow(const Context& context,
                   unsigned bfIdx,
                   unsigned timeIdx,
                   const FluidState& fluidState)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fluidState);

        // we only allow fluxes in the direction opposite to the outer
        // unit normal
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            Scalar& val = this->operator[](eqIdx);
            val = std::min<Scalar>(0.0, val);
        }
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setOutFlow
     */
    template <class Context, class FluidState>
    void setOutFlow(const Context& context,
                    unsigned bfIdx,
                    unsigned timeIdx,
                    const FluidState& fluidState)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fluidState);

        // we only allow fluxes in the same direction as the outer
        // unit normal
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            Scalar& val = this->operator[](eqIdx);
            val = std::max( Scalar(0), val);
        }
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setNoFlow
     */
    void setNoFlow()
    { (*this) = Scalar(0); }

    /*!
     * \copydoc Specify an energy flux that corresponds to the thermal conduction from
     *          the domain boundary
     *
     * This means that a "thermal flow" boundary is a no-flow condition for mass and thermal
     * conduction for energy.
     */
    template <class Context, class FluidState>
    void setThermalFlow(const Context& context,
                        unsigned bfIdx,
                        unsigned timeIdx,
                        const FluidState& boundaryFluidState)
    {
        // set the mass no-flow condition
        setNoFlow();

        if (!enableEnergy)
            // if we do not conserve energy there is nothing we should do in addition
            return;

        ExtensiveQuantities extQuants;
        extQuants.updateBoundary(context, bfIdx, timeIdx, boundaryFluidState);

        (*this)[contiEnergyEqIdx] += extQuants.energyFlux();

#ifndef NDEBUG
        for (unsigned i = 0; i < numEq; ++i)
            Ewoms::Valgrind::CheckDefined((*this)[i]);
        Ewoms::Valgrind::CheckDefined(*this);
#endif
    }
};

} // namespace Ewoms

#endif
