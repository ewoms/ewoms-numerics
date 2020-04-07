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
 * \copydoc Ewoms::PvsIntensiveQuantities
 */
#ifndef EWOMS_PVS_INTENSIVE_QUANTITIES_HH
#define EWOMS_PVS_INTENSIVE_QUANTITIES_HH

#include "pvsproperties.hh"

#include <ewoms/numerics/common/energymodule.hh>
#include <ewoms/numerics/common/diffusionmodule.hh>

#include <ewoms/material/constraintsolvers/computefromreferencephase.hh>
#include <ewoms/material/constraintsolvers/misciblemultiphasecomposition.hh>
#include <ewoms/material/fluidstates/compositionalfluidstate.hh>
#include <ewoms/common/valgrind.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <iostream>

namespace Ewoms {
/*!
 * \ingroup PvsModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the compositional multi-phase primary
 *        variable switching model.
 */
template <class TypeTag>
class PvsIntensiveQuantities
    : public GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities)
    , public DiffusionIntensiveQuantities<TypeTag, GET_PROP_VALUE(TypeTag, EnableDiffusion) >
    , public EnergyIntensiveQuantities<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy) >
    , public GET_PROP_TYPE(TypeTag, FluxModule)::FluxIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluxModule) FluxModule;

    enum { switch0Idx = Indices::switch0Idx };
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { dimWorld = GridView::dimensionworld };

    typedef Ewoms::MiscibleMultiPhaseComposition<Scalar, FluidSystem, Evaluation> MiscibleMultiPhaseComposition;
    typedef Ewoms::ComputeFromReferencePhase<Scalar, FluidSystem, Evaluation> ComputeFromReferencePhase;

    typedef Dune::FieldVector<Evaluation, numPhases> EvalPhaseVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

    typedef typename FluxModule::FluxIntensiveQuantities FluxIntensiveQuantities;
    typedef Ewoms::DiffusionIntensiveQuantities<TypeTag, enableDiffusion> DiffusionIntensiveQuantities;
    typedef Ewoms::EnergyIntensiveQuantities<TypeTag, enableEnergy> EnergyIntensiveQuantities;

public:
    //! The type of the object returned by the fluidState() method
    typedef Ewoms::CompositionalFluidState<Evaluation, FluidSystem> FluidState;

    PvsIntensiveQuantities()
    { }

    PvsIntensiveQuantities(const PvsIntensiveQuantities& other) = default;

    PvsIntensiveQuantities& operator=(const PvsIntensiveQuantities& other) = default;

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::update
     */
    void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);
        EnergyIntensiveQuantities::updateTemperatures_(fluidState_, elemCtx, dofIdx, timeIdx);

        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const auto& problem = elemCtx.problem();

        /////////////
        // set the saturations
        /////////////
        Evaluation sumSat = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fluidState_.setSaturation(phaseIdx, priVars.explicitSaturationValue(phaseIdx, timeIdx));
            Ewoms::Valgrind::CheckDefined(fluidState_.saturation(phaseIdx));
            sumSat += fluidState_.saturation(phaseIdx);
        }
        Ewoms::Valgrind::CheckDefined(priVars.implicitSaturationIdx());
        Ewoms::Valgrind::CheckDefined(sumSat);
        fluidState_.setSaturation(priVars.implicitSaturationIdx(), 1.0 - sumSat);

        /////////////
        // set the pressures of the fluid phases
        /////////////

        // calculate capillary pressure
        const MaterialLawParams& materialParams =
            problem.materialLawParams(elemCtx, dofIdx, timeIdx);
        EvalPhaseVector pC;
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);

        // set the absolute phase pressures in the fluid state
        const Evaluation& p0 = priVars.makeEvaluation(pressure0Idx, timeIdx);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            fluidState_.setPressure(phaseIdx, p0 + (pC[phaseIdx] - pC[0]));

        /////////////
        // calculate the phase compositions
        /////////////

        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        unsigned lowestPresentPhaseIdx = priVars.lowestPresentPhaseIdx();
        unsigned numNonPresentPhases = 0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!priVars.phaseIsPresent(phaseIdx))
                ++numNonPresentPhases;
        }

        // now comes the tricky part: calculate phase compositions
        if (numNonPresentPhases == numPhases - 1) {
            // only one phase is present, i.e. the primary variables
            // contain the complete composition of the phase
            Evaluation sumx = 0.0;
            for (unsigned compIdx = 1; compIdx < numComponents; ++compIdx) {
                const Evaluation& x = priVars.makeEvaluation(switch0Idx + compIdx - 1, timeIdx);
                fluidState_.setMoleFraction(lowestPresentPhaseIdx, compIdx, x);
                sumx += x;
            }

            // set the mole fraction of the first component
            fluidState_.setMoleFraction(lowestPresentPhaseIdx, 0, 1 - sumx);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState_, paramCache,
                                             lowestPresentPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setEnthalpy=*/false);
        }
        else {
            // create the auxiliary constraints
            unsigned numAuxConstraints = numComponents + numNonPresentPhases - numPhases;
            Ewoms::MMPCAuxConstraint<Evaluation> auxConstraints[numComponents];

            unsigned auxIdx = 0;
            unsigned switchIdx = 0;
            for (; switchIdx < numPhases - 1; ++switchIdx) {
                unsigned compIdx = switchIdx + 1;
                unsigned switchPhaseIdx = switchIdx;
                if (switchIdx >= lowestPresentPhaseIdx)
                    switchPhaseIdx += 1;

                if (!priVars.phaseIsPresent(switchPhaseIdx)) {
                    auxConstraints[auxIdx].set(lowestPresentPhaseIdx, compIdx,
                                               priVars.makeEvaluation(switch0Idx + switchIdx, timeIdx));
                    ++auxIdx;
                }
            }

            for (; auxIdx < numAuxConstraints; ++auxIdx, ++switchIdx) {
                unsigned compIdx = numPhases - numNonPresentPhases + auxIdx;
                auxConstraints[auxIdx].set(lowestPresentPhaseIdx, compIdx,
                                           priVars.makeEvaluation(switch0Idx + switchIdx, timeIdx));
            }

            // both phases are present, i.e. phase compositions are a result of the the
            // gas <-> liquid equilibrium. This is the job of the
            // "MiscibleMultiPhaseComposition" constraint solver
            MiscibleMultiPhaseComposition::solve(fluidState_, paramCache,
                                                 priVars.phasePresence(),
                                                 auxConstraints,
                                                 numAuxConstraints,
                                                 /*setViscosity=*/true,
                                                 /*setEnthalpy=*/false);
        }

#ifndef NDEBUG
        // make valgrind happy and set the enthalpies to NaN
        if (!enableEnergy) {
            Scalar myNan = std::numeric_limits<Scalar>::quiet_NaN();
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fluidState_.setEnthalpy(phaseIdx, myNan);
        }
#endif

        /////////////
        // calculate the remaining quantities
        /////////////

        // calculate relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_,
                                            materialParams, fluidState_);
        Ewoms::Valgrind::CheckDefined(relativePermeability_);

        // mobilities
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            mobility_[phaseIdx] =
                relativePermeability_[phaseIdx] / fluidState().viscosity(phaseIdx);

        // porosity
        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);
        Ewoms::Valgrind::CheckDefined(porosity_);

        // intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, dofIdx, timeIdx);

        // update the quantities specific for the velocity model
        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

        // energy related quantities
        EnergyIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        // update the diffusion specific quantities of the intensive quantities
        DiffusionIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        fluidState_.checkDefined();
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::fluidState
     */
    const FluidState& fluidState() const
    { return fluidState_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::intrinsicPermeability
     */
    const DimMatrix& intrinsicPermeability() const
    { return intrinsicPerm_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::relativePermeability
     */
    const Evaluation& relativePermeability(unsigned phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::mobility
     */
    const Evaluation& mobility(unsigned phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::porosity
     */
    const Evaluation& porosity() const
    { return porosity_; }

private:
    FluidState fluidState_;
    Evaluation porosity_;
    DimMatrix intrinsicPerm_;
    Evaluation relativePermeability_[numPhases];
    Evaluation mobility_[numPhases];
};

} // namespace Ewoms

#endif
