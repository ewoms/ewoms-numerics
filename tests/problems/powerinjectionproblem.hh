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
 * \copydoc Ewoms::PowerInjectionProblem
 */
#ifndef EWOMS_POWER_INJECTION_PROBLEM_HH
#define EWOMS_POWER_INJECTION_PROBLEM_HH

#include <ewoms/numerics/models/immiscible/immisciblemodel.hh>
#include <ewoms/numerics/io/cubegridvanguard.hh>

#include <ewoms/material/fluidmatrixinteractions/regularizedvangenuchten.hh>
#include <ewoms/material/fluidmatrixinteractions/linearmaterial.hh>
#include <ewoms/material/fluidmatrixinteractions/efftoabslaw.hh>
#include <ewoms/material/fluidmatrixinteractions/materialtraits.hh>
#include <ewoms/material/fluidsystems/twophaseimmisciblefluidsystem.hh>
#include <ewoms/material/fluidstates/immisciblefluidstate.hh>
#include <ewoms/material/components/simpleh2o.hh>
#include <ewoms/material/components/air.hh>
#include <ewoms/common/unused.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <string>
#include <type_traits>
#include <iostream>

namespace Ewoms {
template <class TypeTag>
class PowerInjectionProblem;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(PowerInjectionBaseProblem);

// Set the grid implementation to be used
SET_TYPE_PROP(PowerInjectionBaseProblem, Grid, Dune::YaspGrid</*dim=*/1>);

// set the Vanguard property
SET_TYPE_PROP(PowerInjectionBaseProblem, Vanguard,
              Ewoms::CubeGridVanguard<TypeTag>);

// Set the problem property
SET_TYPE_PROP(PowerInjectionBaseProblem, Problem,
              Ewoms::PowerInjectionProblem<TypeTag>);

// Set the wetting phase
SET_PROP(PowerInjectionBaseProblem, WettingPhase)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);

public:
    using type = Ewoms::LiquidPhase<Scalar, Ewoms::SimpleH2O<Scalar> >;
};

// Set the non-wetting phase
SET_PROP(PowerInjectionBaseProblem, NonwettingPhase)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);

public:
    using type = Ewoms::GasPhase<Scalar, Ewoms::Air<Scalar> >;
};

// Set the material Law
SET_PROP(PowerInjectionBaseProblem, MaterialLaw)
{
private:
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
    enum { wettingPhaseIdx = FluidSystem::wettingPhaseIdx };
    enum { nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx };

    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    typedef Ewoms::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::wettingPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::nonWettingPhaseIdx>
        Traits;

    // define the material law which is parameterized by effective
    // saturations
    using EffectiveLaw = Ewoms::RegularizedVanGenuchten<Traits>;

public:
    // define the material law parameterized by absolute saturations
    using type = Ewoms::EffToAbsLaw<EffectiveLaw>;
};

// Write out the filter velocities for this problem
SET_BOOL_PROP(PowerInjectionBaseProblem, VtkWriteFilterVelocities, true);

// Disable gravity
SET_BOOL_PROP(PowerInjectionBaseProblem, EnableGravity, false);

// define the properties specific for the power injection problem
SET_SCALAR_PROP(PowerInjectionBaseProblem, DomainSizeX, 100.0);
SET_SCALAR_PROP(PowerInjectionBaseProblem, DomainSizeY, 1.0);
SET_SCALAR_PROP(PowerInjectionBaseProblem, DomainSizeZ, 1.0);

SET_INT_PROP(PowerInjectionBaseProblem, CellsX, 250);
SET_INT_PROP(PowerInjectionBaseProblem, CellsY, 1);
SET_INT_PROP(PowerInjectionBaseProblem, CellsZ, 1);

// The default for the end time of the simulation
SET_SCALAR_PROP(PowerInjectionBaseProblem, EndTime, 100);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(PowerInjectionBaseProblem, InitialTimeStepSize, 1e-3);

END_PROPERTIES

namespace Ewoms {
/*!
 * \ingroup TestProblems
 * \brief 1D Problem with very fast injection of gas on the left.
 *
 * The velocity model is chosen in the .cc file in this problem. The
 * spatial parameters are inspired by the ones given by
 *
 * V. Jambhekar: "Forchheimer Porous-media Flow models -- Numerical
 * Investigation and Comparison with Experimental Data", Master's
 * Thesis at Institute for Modelling Hydraulic and Environmental
 * Systems, University of Stuttgart, 2011
 */
template <class TypeTag>
class PowerInjectionProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    using ParentType = GET_PROP_TYPE(TypeTag, BaseProblem);

    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = GET_PROP_TYPE(TypeTag, GridView);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
    using WettingPhase = GET_PROP_TYPE(TypeTag, WettingPhase);
    using NonwettingPhase = GET_PROP_TYPE(TypeTag, NonwettingPhase);
    using PrimaryVariables = GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using EqVector = GET_PROP_TYPE(TypeTag, EqVector);
    using RateVector = GET_PROP_TYPE(TypeTag, RateVector);
    using BoundaryRateVector = GET_PROP_TYPE(TypeTag, BoundaryRateVector);
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);

    enum {
        // number of phases
        numPhases = FluidSystem::numPhases,

        // phase indices
        wettingPhaseIdx = FluidSystem::wettingPhaseIdx,
        nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx,

        // equation indices
        contiNEqIdx = Indices::conti0EqIdx + nonWettingPhaseIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using MaterialLaw = GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = GET_PROP_TYPE(TypeTag, MaterialLawParams);

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    PowerInjectionProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 3e-6;
        FluidSystem::init();

        temperature_ = 273.15 + 26.6;

        // parameters for the Van Genuchten law
        // alpha and n
        materialParams_.setVgAlpha(0.00045);
        materialParams_.setVgN(7.3);
        materialParams_.finalize();

        K_ = this->toDimMatrix_(5.73e-08); // [m^2]

        setupInitialFluidState_();
    }

    /*!
     * \name Auxiliary methods
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << "powerinjection_";
        if (std::is_same<GET_PROP_TYPE(TypeTag, FluxModule),
                         Ewoms::DarcyFluxModule<TypeTag> >::value)
            oss << "darcy";
        else
            oss << "forchheimer";

        if (std::is_same<GET_PROP_TYPE(TypeTag, LocalLinearizerSplice),
                         TTAG(AutoDiffLocalLinearizer)>::value)
            oss << "_" << "ad";
        else
            oss << "_" << "fd";

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
    //! \}

    /*!
     * \name Soil parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context EWOMS_UNUSED,
                                           unsigned spaceIdx EWOMS_UNUSED,
                                           unsigned timeIdx EWOMS_UNUSED) const
    { return K_; }

    /*!
     * \copydoc ForchheimerBaseProblem::ergunCoefficient
     */
    template <class Context>
    Scalar ergunCoefficient(const Context& context EWOMS_UNUSED,
                            unsigned spaceIdx EWOMS_UNUSED,
                            unsigned timeIdx EWOMS_UNUSED) const
    { return 0.3866; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context EWOMS_UNUSED,
                    unsigned spaceIdx EWOMS_UNUSED,
                    unsigned timeIdx EWOMS_UNUSED) const
    { return 0.558; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams&
    materialLawParams(const Context& context EWOMS_UNUSED,
                      unsigned spaceIdx EWOMS_UNUSED,
                      unsigned timeIdx EWOMS_UNUSED) const
    { return materialParams_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context EWOMS_UNUSED,
                       unsigned spaceIdx EWOMS_UNUSED,
                       unsigned timeIdx EWOMS_UNUSED) const
    { return temperature_; }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     *
     * This problem sets a very high injection rate of nitrogen on the
     * left and a free-flow boundary on the right.
     */
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& context,
                  unsigned spaceIdx,
                  unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos)) {
            RateVector massRate(0.0);
            massRate = 0.0;
            massRate[contiNEqIdx] = -1.00; // kg / (m^2 * s)

            // impose a forced flow boundary
            values.setMassRate(massRate);
        }
        else if (onRightBoundary_(pos))
            // free flow boundary with initial condition on the right
            values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidState_);
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
                 const Context& context EWOMS_UNUSED,
                 unsigned spaceIdx EWOMS_UNUSED,
                 unsigned timeIdx EWOMS_UNUSED) const
    {
        // assign the primary variables
        values.assignNaive(initialFluidState_);
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
    { return pos[0] < this->boundingBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    void setupInitialFluidState_()
    {
        initialFluidState_.setTemperature(temperature_);

        Scalar Sw = 1.0;
        initialFluidState_.setSaturation(wettingPhaseIdx, Sw);
        initialFluidState_.setSaturation(nonWettingPhaseIdx, 1 - Sw);

        Scalar p = 1e5;
        initialFluidState_.setPressure(wettingPhaseIdx, p);
        initialFluidState_.setPressure(nonWettingPhaseIdx, p);

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updateAll(initialFluidState_);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            initialFluidState_.setDensity(phaseIdx,
                                          FluidSystem::density(initialFluidState_, paramCache, phaseIdx));
            initialFluidState_.setViscosity(phaseIdx,
                                            FluidSystem::viscosity(initialFluidState_, paramCache, phaseIdx));
        }
    }

    DimMatrix K_;
    MaterialLawParams materialParams_;

    Ewoms::ImmiscibleFluidState<Scalar, FluidSystem> initialFluidState_;
    Scalar temperature_;
    Scalar eps_;
};

} // namespace Ewoms

#endif
