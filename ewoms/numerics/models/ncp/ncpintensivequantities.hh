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
 * \copydoc Ewoms::NcpIntensiveQuantities
 */
#ifndef EWOMS_NCP_INTENSIVE_QUANTITIES_HH
#define EWOMS_NCP_INTENSIVE_QUANTITIES_HH

#include "ncpproperties.hh"

#include <ewoms/numerics/common/energymodule.hh>
#include <ewoms/numerics/common/diffusionmodule.hh>

#include <ewoms/material/constraintsolvers/ncpflash.hh>
#include <ewoms/material/fluidstates/compositionalfluidstate.hh>
#include <ewoms/material/constraintsolvers/compositionfromfugacities.hh>
#include <ewoms/common/valgrind.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Ewoms {
/*!
 * \ingroup NcpModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the compositional multi-phase NCP model.
 */
template <class TypeTag>
class NcpIntensiveQuantities
    : public GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities)
    , public DiffusionIntensiveQuantities<TypeTag, GET_PROP_VALUE(TypeTag, EnableDiffusion) >
    , public EnergyIntensiveQuantities<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy) >
    , public GET_PROP_TYPE(TypeTag, FluxModule)::FluxIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluxModule) FluxModule;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    enum { fugacity0Idx = Indices::fugacity0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { dimWorld = GridView::dimensionworld };

    typedef Ewoms::CompositionFromFugacities<Scalar, FluidSystem, Evaluation>
        CompositionFromFugacitiesSolver;
    typedef Ewoms::CompositionalFluidState<Evaluation, FluidSystem,
                                         /*storeEnthalpy=*/enableEnergy> FluidState;
    typedef Dune::FieldVector<Evaluation, numComponents> ComponentVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Ewoms::DiffusionIntensiveQuantities<TypeTag, enableDiffusion> DiffusionIntensiveQuantities;
    typedef Ewoms::EnergyIntensiveQuantities<TypeTag, enableEnergy> EnergyIntensiveQuantities;
    typedef typename FluxModule::FluxIntensiveQuantities FluxIntensiveQuantities;

public:
    NcpIntensiveQuantities()
    {}

    NcpIntensiveQuantities(const NcpIntensiveQuantities& other) = default;

    NcpIntensiveQuantities& operator=(const NcpIntensiveQuantities& other)  = default;

    /*!
     * \brief IntensiveQuantities::update
     */
    void update(const ElementContext& elemCtx,
                unsigned dofIdx,
                unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);
        ParentType::checkDefined();

        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);

        // set the phase saturations
        Evaluation sumSat = 0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
            const Evaluation& val = priVars.makeEvaluation(saturation0Idx + phaseIdx, timeIdx);
            fluidState_.setSaturation(phaseIdx, val);
            sumSat += val;
        }
        fluidState_.setSaturation(numPhases - 1, 1.0 - sumSat);
        Ewoms::Valgrind::CheckDefined(sumSat);

        // set the fluid phase temperature
        EnergyIntensiveQuantities::updateTemperatures_(fluidState_, elemCtx, dofIdx, timeIdx);

        // retrieve capillary pressure parameters
        const auto& problem = elemCtx.problem();
        const MaterialLawParams& materialParams =
            problem.materialLawParams(elemCtx, dofIdx, timeIdx);
        // calculate capillary pressures
        Evaluation capPress[numPhases];
        MaterialLaw::capillaryPressures(capPress, materialParams, fluidState_);
        // add to the pressure of the first fluid phase
        const Evaluation& pressure0 = priVars.makeEvaluation(pressure0Idx, timeIdx);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            fluidState_.setPressure(phaseIdx, pressure0 + (capPress[phaseIdx] - capPress[0]));

        ComponentVector fug;
        // retrieve component fugacities
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            fug[compIdx] = priVars.makeEvaluation(fugacity0Idx + compIdx, timeIdx);

        // calculate phase compositions
        const auto *hint = elemCtx.thermodynamicHint(dofIdx, timeIdx);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // initial guess
            if (hint) {
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    // use the hint for the initial mole fraction!
                    const Evaluation& moleFracIJ = hint->fluidState().moleFraction(phaseIdx, compIdx);
                    fluidState_.setMoleFraction(phaseIdx, compIdx, moleFracIJ);
                }
            }
            else // !hint
                CompositionFromFugacitiesSolver::guessInitial(fluidState_, phaseIdx, fug);

            // calculate the phase composition from the component
            // fugacities
            CompositionFromFugacitiesSolver::solve(fluidState_, paramCache, phaseIdx, fug);
        }

        // porosity
        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);
        Ewoms::Valgrind::CheckDefined(porosity_);

        // relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_, materialParams, fluidState_);

        // dynamic viscosities
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // viscosities
            const Evaluation& mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            fluidState_.setViscosity(phaseIdx, mu);

            mobility_[phaseIdx] = relativePermeability_[phaseIdx]/mu;
        }

        // intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, dofIdx, timeIdx);

        // update the quantities specific for the velocity model
        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

        // energy related quantities
        EnergyIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        // update the diffusion specific quantities of the intensive quantities
        DiffusionIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        checkDefined();
    }

    /*!
     * \brief ImmiscibleIntensiveQuantities::fluidState
     */
    const FluidState& fluidState() const
    { return fluidState_; }

    /*!
     * \brief ImmiscibleIntensiveQuantities::intrinsicPermeability
     */
    const DimMatrix& intrinsicPermeability() const
    { return intrinsicPerm_; }

    /*!
     * \brief ImmiscibleIntensiveQuantities::relativePermeability
     */
    const Evaluation& relativePermeability(unsigned phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \brief ImmiscibleIntensiveQuantities::mobility
     */
    const Evaluation& mobility(unsigned phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \brief ImmiscibleIntensiveQuantities::porosity
     */
    const Evaluation& porosity() const
    { return porosity_; }

    /*!
     * \brief IntensiveQuantities::checkDefined
     */
    void checkDefined() const
    {
#if !defined NDEBUG && HAVE_VALGRIND
        ParentType::checkDefined();

        Ewoms::Valgrind::CheckDefined(porosity_);
        Ewoms::Valgrind::CheckDefined(relativePermeability_);

        fluidState_.checkDefined();
#endif
    }

private:
    DimMatrix intrinsicPerm_;
    FluidState fluidState_;
    Evaluation porosity_;
    Evaluation relativePermeability_[numPhases];
    Evaluation mobility_[numPhases];
};

} // namespace Ewoms

#endif
