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
 * \copydoc Ewoms::InfiltrationProblem
 */
#ifndef EWOMS_INFILTRATION_PROBLEM_HH
#define EWOMS_INFILTRATION_PROBLEM_HH

#include <ewoms/numerics/models/pvs/pvsproperties.hh>

#include <ewoms/material/fluidstates/compositionalfluidstate.hh>
#include <ewoms/material/fluidsystems/h2oairmesitylenefluidsystem.hh>
#include <ewoms/material/fluidmatrixinteractions/threephaseparkervangenuchten.hh>
#include <ewoms/material/fluidmatrixinteractions/materialtraits.hh>
#include <ewoms/material/constraintsolvers/computefromreferencephase.hh>
#include <ewoms/common/valgrind.hh>
#include <ewoms/common/unused.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <string>

namespace Ewoms {
template <class TypeTag>
class InfiltrationProblem;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(InfiltrationBaseProblem);

// Set the grid type
SET_TYPE_PROP(InfiltrationBaseProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(InfiltrationBaseProblem, Problem,
              Ewoms::InfiltrationProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(
    InfiltrationBaseProblem, FluidSystem,
    Ewoms::H2OAirMesityleneFluidSystem<GET_PROP_TYPE(TypeTag, Scalar)>);

// Enable gravity?
SET_BOOL_PROP(InfiltrationBaseProblem, EnableGravity, true);

// Write newton convergence?
SET_BOOL_PROP(InfiltrationBaseProblem, NewtonWriteConvergence, false);

// -1 backward differences, 0: central differences, +1: forward differences
SET_INT_PROP(InfiltrationBaseProblem, NumericDifferenceMethod, 1);

// Set the material Law
SET_PROP(InfiltrationBaseProblem, MaterialLaw)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);

    typedef Ewoms::ThreePhaseMaterialTraits<
        Scalar,
        /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
        /*nonWettingPhaseIdx=*/FluidSystem::naplPhaseIdx,
        /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx> Traits;

public:
    using type = Ewoms::ThreePhaseParkerVanGenuchten<Traits>;
};

// The default for the end time of the simulation
SET_SCALAR_PROP(InfiltrationBaseProblem, EndTime, 6e3);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(InfiltrationBaseProblem, InitialTimeStepSize, 60);

// The default DGF file to load
SET_STRING_PROP(InfiltrationBaseProblem, GridFile,
                "./data/infiltration_50x3.dgf");

END_PROPERTIES

namespace Ewoms {
/*!
 * \ingroup TestProblems
 * \brief Isothermal NAPL infiltration problem where LNAPL
 *        contaminates the unsaturated and the saturated groundwater
 *        zone.
 *
 * The 2D domain of this test problem is 500 m long and 10 m deep,
 * where the lower part represents a slightly inclined groundwater
 * table, and the upper part is the vadose zone.  A LNAPL (Non-Aqueous
 * Phase Liquid which is lighter than water) infiltrates (modelled
 * with a Neumann boundary condition) into the vadose zone. Upon
 * reaching the water table, it spreads (since lighter than water) and
 * migrates on top of the water table in the direction of the slope.
 * On its way through the vadose zone, it leaves a trace of residually
 * trapped immobile NAPL, which can in the following dissolve and
 * evaporate slowly, and eventually be transported by advection and
 * diffusion.
 *
 * Left and right boundaries are constant hydraulic head boundaries
 * (Dirichlet), Top and bottom are Neumann boundaries, all no-flow
 * except for the small infiltration zone in the upper left part.
 */
template <class TypeTag>
class InfiltrationProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    using ParentType = GET_PROP_TYPE(TypeTag, BaseProblem);

    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = GET_PROP_TYPE(TypeTag, GridView);
    using MaterialLaw = GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = GET_PROP_TYPE(TypeTag, MaterialLawParams);
    using PrimaryVariables = GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using EqVector = GET_PROP_TYPE(TypeTag, EqVector);
    using RateVector = GET_PROP_TYPE(TypeTag, RateVector);
    using BoundaryRateVector = GET_PROP_TYPE(TypeTag, BoundaryRateVector);
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
    using Model = GET_PROP_TYPE(TypeTag, Model);

    // copy some indices for convenience
    using Indices = GET_PROP_TYPE(TypeTag, Indices);
    enum {
        // equation indices
        conti0EqIdx = Indices::conti0EqIdx,

        // number of phases/components
        numPhases = FluidSystem::numPhases,

        // component indices
        NAPLIdx = FluidSystem::NAPLIdx,
        H2OIdx = FluidSystem::H2OIdx,
        airIdx = FluidSystem::airIdx,

        // phase indices
        waterPhaseIdx = FluidSystem::waterPhaseIdx,
        gasPhaseIdx = FluidSystem::gasPhaseIdx,
        naplPhaseIdx = FluidSystem::naplPhaseIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    InfiltrationProblem(Simulator& simulator)
        : ParentType(simulator)
        , eps_(1e-6)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        temperature_ = 273.15 + 10.0; // -> 10 degrees Celsius
        FluidSystem::init(/*tempMin=*/temperature_ - 1,
                          /*tempMax=*/temperature_ + 1,
                          /*nTemp=*/3,
                          /*pressMin=*/0.8 * 1e5,
                          /*pressMax=*/3 * 1e5,
                          /*nPress=*/200);

        // intrinsic permeabilities
        fineK_ = this->toDimMatrix_(1e-11);
        coarseK_ = this->toDimMatrix_(1e-11);

        // porosities
        porosity_ = 0.40;

        // residual saturations
        materialParams_.setSwr(0.12);
        materialParams_.setSwrx(0.12);
        materialParams_.setSnr(0.07);
        materialParams_.setSgr(0.03);

        // parameters for the three-phase van Genuchten law
        materialParams_.setVgAlpha(0.0005);
        materialParams_.setVgN(4.);
        materialParams_.setkrRegardsSnr(false);

        materialParams_.finalize();
        materialParams_.checkDefined();
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::shouldWriteRestartFile
     *
     * This problem writes a restart file after every time step.
     */
    bool shouldWriteRestartFile() const
    { return true; }

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << "infiltration_" << Model::name();
        return oss.str();
    }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        this->model().checkConservativeness();

        // Calculate storage terms
        EqVector storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage: " << storage << std::endl << std::flush;
        }
#endif // NDEBUG
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context EWOMS_UNUSED,
                       unsigned spaceIdx EWOMS_UNUSED,
                       unsigned timeIdx EWOMS_UNUSED) const
    { return temperature_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix&
    intrinsicPermeability(const Context& context,
                          unsigned spaceIdx,
                          unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context EWOMS_UNUSED,
                    unsigned spaceIdx EWOMS_UNUSED,
                    unsigned timeIdx EWOMS_UNUSED) const
    { return porosity_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams&
    materialLawParams(const Context& context EWOMS_UNUSED,
                      unsigned spaceIdx EWOMS_UNUSED,
                      unsigned timeIdx EWOMS_UNUSED) const
    { return materialParams_; }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& context,
                  unsigned spaceIdx,
                  unsigned timeIdx) const
    {
        const auto& pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos) || onRightBoundary_(pos)) {
            Ewoms::CompositionalFluidState<Scalar, FluidSystem> fs;

            initialFluidState_(fs, context, spaceIdx, timeIdx);

            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else if (onInlet_(pos)) {
            RateVector molarRate(0.0);
            molarRate[conti0EqIdx + NAPLIdx] = -0.001;

            values.setMolarRate(molarRate);
            Ewoms::Valgrind::CheckDefined(values);
        }
        else
            values.setNoFlow();
    }

    //! \}

    /*!
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values,
                 const Context& context,
                 unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        Ewoms::CompositionalFluidState<Scalar, FluidSystem> fs;

        initialFluidState_(fs, context, spaceIdx, timeIdx);

        const auto& matParams = materialLawParams(context, spaceIdx, timeIdx);
        values.assignMassConservative(fs, matParams, /*inEquilibrium=*/true);
        Ewoms::Valgrind::CheckDefined(values);
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& context EWOMS_UNUSED,
                unsigned spaceIdx EWOMS_UNUSED,
                unsigned timeIdx EWOMS_UNUSED) const
    { rate = Scalar(0.0); }

    //! \}

private:
    bool onLeftBoundary_(const GlobalPosition& pos) const
    { return pos[0] < eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition& pos) const
    { return pos[1] < eps_; }

    bool onUpperBoundary_(const GlobalPosition& pos) const
    { return pos[1] > this->boundingBoxMax()[1] - eps_; }

    bool onInlet_(const GlobalPosition& pos) const
    { return onUpperBoundary_(pos) && 50 < pos[0] && pos[0] < 75; }

    template <class FluidState, class Context>
    void initialFluidState_(FluidState& fs, const Context& context,
                            unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition pos = context.pos(spaceIdx, timeIdx);
        Scalar y = pos[1];
        Scalar x = pos[0];

        Scalar densityW = 1000.0;
        Scalar pc = 9.81 * densityW * (y - (5 - 5e-4 * x));
        if (pc < 0.0)
            pc = 0.0;

        // set pressures
        const auto& matParams = materialLawParams(context, spaceIdx, timeIdx);
        Scalar Sw = matParams.Swr();
        Scalar Swr = matParams.Swr();
        Scalar Sgr = matParams.Sgr();
        if (Sw < Swr)
            Sw = Swr;
        if (Sw > 1 - Sgr)
            Sw = 1 - Sgr;
        Scalar Sg = 1 - Sw;

        Ewoms::Valgrind::CheckDefined(Sw);
        Ewoms::Valgrind::CheckDefined(Sg);

        fs.setSaturation(waterPhaseIdx, Sw);
        fs.setSaturation(gasPhaseIdx, Sg);
        fs.setSaturation(naplPhaseIdx, 0);

        // set temperature of all phases
        fs.setTemperature(temperature_);

        // compute pressures
        Scalar pcAll[numPhases];
        Scalar pg = 1e5;
        if (onLeftBoundary_(pos))
            pg += 10e3;
        MaterialLaw::capillaryPressures(pcAll, matParams, fs);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            fs.setPressure(phaseIdx, pg + (pcAll[phaseIdx] - pcAll[gasPhaseIdx]));

        // set composition of gas phase
        fs.setMoleFraction(gasPhaseIdx, H2OIdx, 1e-6);
        fs.setMoleFraction(gasPhaseIdx, airIdx,
                           1 - fs.moleFraction(gasPhaseIdx, H2OIdx));
        fs.setMoleFraction(gasPhaseIdx, NAPLIdx, 0);

        using CFRP = Ewoms::ComputeFromReferencePhase<Scalar, FluidSystem>;
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        CFRP::solve(fs, paramCache, gasPhaseIdx,
                    /*setViscosity=*/true,
                    /*setEnthalpy=*/false);

        fs.setMoleFraction(waterPhaseIdx, H2OIdx,
                           1 - fs.moleFraction(waterPhaseIdx, H2OIdx));
    }

    bool isFineMaterial_(const GlobalPosition& pos) const
    {  return 70. <= pos[0] && pos[0] <= 85. && 7.0 <= pos[1] && pos[1] <= 7.50; }

    DimMatrix fineK_;
    DimMatrix coarseK_;

    Scalar porosity_;

    MaterialLawParams materialParams_;

    Scalar temperature_;
    Scalar eps_;
};
} // namespace Ewoms

#endif
