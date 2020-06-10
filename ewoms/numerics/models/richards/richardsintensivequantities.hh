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
 * \copydoc Ewoms::RichardsIntensiveQuantities
 */
#ifndef EWOMS_RICHARDS_INTENSIVE_QUANTITIES_HH
#define EWOMS_RICHARDS_INTENSIVE_QUANTITIES_HH

#include "richardsproperties.hh"

#include <ewoms/material/fluidstates/immisciblefluidstate.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Ewoms {

/*!
 * \ingroup RichardsModel
 * \ingroup IntensiveQuantities
 *
 * \brief Intensive quantities required by the Richards model.
 */
template <class TypeTag>
class RichardsIntensiveQuantities
    : public GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities)
    , public GET_PROP_TYPE(TypeTag, FluxModule)::FluxIntensiveQuantities
{
    typedef GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities) ParentType;
    typedef GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef GET_PROP_TYPE(TypeTag, FluxModule) FluxModule;

    typedef GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { pressureWIdx = Indices::pressureWIdx };
    enum { numPhases = FluidSystem::numPhases };
    enum { liquidPhaseIdx = GET_PROP_VALUE(TypeTag, LiquidPhaseIndex) };
    enum { gasPhaseIdx = GET_PROP_VALUE(TypeTag, GasPhaseIndex) };
    enum { dimWorld = GridView::dimensionworld };

    typedef typename FluxModule::FluxIntensiveQuantities FluxIntensiveQuantities;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<Scalar, numPhases> ScalarPhaseVector;
    typedef Dune::FieldVector<Evaluation, numPhases> PhaseVector;
    typedef Ewoms::MathToolbox<Evaluation> Toolbox;

public:
    //! The type returned by the fluidState() method
    typedef Ewoms::ImmiscibleFluidState<Evaluation, FluidSystem> FluidState;

    RichardsIntensiveQuantities()
    {}

    RichardsIntensiveQuantities(const RichardsIntensiveQuantities& other) = default;

    RichardsIntensiveQuantities& operator=(const RichardsIntensiveQuantities& other) = default;

    /*!
     * \copydoc IntensiveQuantities::update
     */
    void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);

        const auto& T = elemCtx.problem().temperature(elemCtx, dofIdx, timeIdx);
        fluidState_.setTemperature(T);

        // material law parameters
        const auto& problem = elemCtx.problem();
        const typename MaterialLaw::Params& materialParams =
            problem.materialLawParams(elemCtx, dofIdx, timeIdx);
        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);

        /////////
        // calculate the pressures
        /////////

        // first, we have to find the minimum capillary pressure (i.e. Sw = 0)
        fluidState_.setSaturation(liquidPhaseIdx, 1.0);
        fluidState_.setSaturation(gasPhaseIdx, 0.0);
        ScalarPhaseVector pC;
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);

        // non-wetting pressure can be larger than the
        // reference pressure if the medium is fully
        // saturated by the wetting phase
        const Evaluation& pW = priVars.makeEvaluation(pressureWIdx, timeIdx);
        Evaluation pN =
            Toolbox::max(elemCtx.problem().referencePressure(elemCtx, dofIdx, /*timeIdx=*/0),
                         pW + (pC[gasPhaseIdx] - pC[liquidPhaseIdx]));

        /////////
        // calculate the saturations
        /////////
        fluidState_.setPressure(liquidPhaseIdx, pW);
        fluidState_.setPressure(gasPhaseIdx, pN);

        PhaseVector sat;
        MaterialLaw::saturations(sat, materialParams, fluidState_);
        fluidState_.setSaturation(liquidPhaseIdx, sat[liquidPhaseIdx]);
        fluidState_.setSaturation(gasPhaseIdx, sat[gasPhaseIdx]);

        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        paramCache.updateAll(fluidState_);

        // compute and set the wetting phase viscosity
        const Evaluation& mu = FluidSystem::viscosity(fluidState_, paramCache, liquidPhaseIdx);
        fluidState_.setViscosity(liquidPhaseIdx, mu);
        fluidState_.setViscosity(gasPhaseIdx, 1e-20);

        // compute and set the wetting phase density
        const Evaluation& rho = FluidSystem::density(fluidState_, paramCache, liquidPhaseIdx);
        fluidState_.setDensity(liquidPhaseIdx, rho);
        fluidState_.setDensity(gasPhaseIdx, 1e-20);

        // relperms
        MaterialLaw::relativePermeabilities(relativePermeability_, materialParams, fluidState_);

        // mobilities
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            mobility_[phaseIdx] = relativePermeability_[phaseIdx]/fluidState_.viscosity(phaseIdx);

        // porosity
        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);

        // intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, dofIdx, timeIdx);

        // update the quantities specific for the velocity model
        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::fluidState
     */
    const FluidState& fluidState() const
    { return fluidState_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::porosity
     */
    const Evaluation& porosity() const
    { return porosity_; }

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

private:
    FluidState fluidState_;
    DimMatrix intrinsicPerm_;
    Evaluation relativePermeability_[numPhases];
    Evaluation mobility_[numPhases];
    Evaluation porosity_;
};

} // namespace Ewoms

#endif
