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
 * \copydoc Ewoms::FingerProblem
 */
#ifndef EWOMS_FINGER_PROBLEM_HH
#define EWOMS_FINGER_PROBLEM_HH

#include <ewoms/numerics/io/structuredgridvanguard.hh>

#include <ewoms/material/fluidmatrixinteractions/regularizedvangenuchten.hh>
#include <ewoms/material/fluidmatrixinteractions/linearmaterial.hh>
#include <ewoms/material/fluidmatrixinteractions/efftoabslaw.hh>
#include <ewoms/material/fluidmatrixinteractions/parkerlenhard.hh>
#include <ewoms/material/fluidmatrixinteractions/materialtraits.hh>

#include <ewoms/material/fluidsystems/twophaseimmisciblefluidsystem.hh>
#include <ewoms/material/fluidstates/immisciblefluidstate.hh>
#include <ewoms/material/components/simpleh2o.hh>
#include <ewoms/material/components/air.hh>

#include <ewoms/numerics/models/immiscible/immiscibleproperties.hh>
#include <ewoms/numerics/discretizations/common/restrictprolong.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <vector>
#include <string>

namespace Ewoms {
template <class TypeTag>
class FingerProblem;

} // namespace Ewoms

BEGIN_PROPERTIES

NEW_TYPE_TAG(FingerBaseProblem, INHERITS_FROM(StructuredGridVanguard));

#if HAVE_DUNE_ALUGRID
// use dune-alugrid if available
SET_TYPE_PROP(FingerBaseProblem,
              Grid,
              Dune::ALUGrid</*dim=*/2,
                            /*dimWorld=*/2,
                            Dune::cube,
                            Dune::nonconforming>);
#endif

// declare the properties used by the finger problem
NEW_PROP_TAG(InitialWaterSaturation);

// Set the problem property
SET_TYPE_PROP(FingerBaseProblem, Problem, Ewoms::FingerProblem<TypeTag>);

// Set the wetting phase
SET_PROP(FingerBaseProblem, WettingPhase)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);

public:
    using type = Ewoms::LiquidPhase<Scalar, Ewoms::SimpleH2O<Scalar> >;
};

// Set the non-wetting phase
SET_PROP(FingerBaseProblem, NonwettingPhase)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);

public:
    using type = Ewoms::GasPhase<Scalar, Ewoms::Air<Scalar> >;
};

// Set the material Law
SET_PROP(FingerBaseProblem, MaterialLaw)
{
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
    typedef Ewoms::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::wettingPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::nonWettingPhaseIdx> Traits;

    // use the parker-lenhard hysteresis law
    using ParkerLenhard = Ewoms::ParkerLenhard<Traits>;
    using type = ParkerLenhard;
};

// Write the solutions of individual newton iterations?
SET_BOOL_PROP(FingerBaseProblem, NewtonWriteConvergence, false);

// Use forward differences instead of central differences
SET_INT_PROP(FingerBaseProblem, NumericDifferenceMethod, +1);

// Enable constraints
SET_INT_PROP(FingerBaseProblem, EnableConstraints, true);

// Enable gravity
SET_BOOL_PROP(FingerBaseProblem, EnableGravity, true);

// define the properties specific for the finger problem
SET_SCALAR_PROP(FingerBaseProblem, DomainSizeX, 0.1);
SET_SCALAR_PROP(FingerBaseProblem, DomainSizeY, 0.3);
SET_SCALAR_PROP(FingerBaseProblem, DomainSizeZ, 0.1);

SET_SCALAR_PROP(FingerBaseProblem, InitialWaterSaturation, 0.01);

SET_INT_PROP(FingerBaseProblem, CellsX, 20);
SET_INT_PROP(FingerBaseProblem, CellsY, 70);
SET_INT_PROP(FingerBaseProblem, CellsZ, 1);

// The default for the end time of the simulation
SET_SCALAR_PROP(FingerBaseProblem, EndTime, 215);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(FingerBaseProblem, InitialTimeStepSize, 10);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup TestProblems
 *
 * \brief Two-phase problem featuring some gravity-driven saturation
 *        fingers.
 *
 * The domain of this problem is sized 10cm times 1m and is initially
 * dry. Water is then injected at three locations on the top of the
 * domain which leads to gravity fingering. The boundary conditions
 * used are no-flow for the left and right and top of the domain and
 * free-flow at the bottom. This problem uses the Parker-Lenhard
 * hystersis model which might lead to non-monotonic saturation in the
 * fingers if the right material parameters is chosen and the spatial
 * discretization is fine enough.
 */
template <class TypeTag>
class FingerProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    //!\cond SKIP_THIS
    using ParentType = GET_PROP_TYPE(TypeTag, BaseProblem);

    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = GET_PROP_TYPE(TypeTag, GridView);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
    using WettingPhase = GET_PROP_TYPE(TypeTag, WettingPhase);
    using NonwettingPhase = GET_PROP_TYPE(TypeTag, NonwettingPhase);
    using PrimaryVariables = GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using Constraints = GET_PROP_TYPE(TypeTag, Constraints);
    using Model = GET_PROP_TYPE(TypeTag, Model);

    enum {
        // number of phases
        numPhases = FluidSystem::numPhases,

        // phase indices
        wettingPhaseIdx = FluidSystem::wettingPhaseIdx,
        nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx,

        // equation indices
        contiWettingEqIdx = Indices::conti0EqIdx + wettingPhaseIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);
    using Stencil = GET_PROP_TYPE(TypeTag, Stencil) ;
    enum { codim = Stencil::Entity::codimension };
    using EqVector = GET_PROP_TYPE(TypeTag, EqVector);
    using RateVector = GET_PROP_TYPE(TypeTag, RateVector);
    using BoundaryRateVector = GET_PROP_TYPE(TypeTag, BoundaryRateVector);

    using ParkerLenhard = GET_PROP(TypeTag, MaterialLaw)::ParkerLenhard;
    using MaterialLaw = GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = GET_PROP_TYPE(TypeTag, MaterialLawParams);

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using Grid = typename GridView :: Grid;

    using MaterialLawParamsContainer = Dune::PersistentContainer< Grid, std::shared_ptr< MaterialLawParams > >  ;
    //!\endcond

public:
    using RestrictProlongOperator = CopyRestrictProlong< Grid, MaterialLawParamsContainer >;

    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    FingerProblem(Simulator& simulator)
        : ParentType(simulator),
          materialParams_( simulator.vanguard().grid(), codim )
    {
    }

    /*!
     * \name Auxiliary methods
     */
    //! \{

    /*!
     * \brief \copydoc FvBaseProblem::restrictProlongOperator
     */
    RestrictProlongOperator restrictProlongOperator()
    {
        return RestrictProlongOperator( materialParams_ );
    }

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    { return
            std::string("finger") +
            "_" + Model::name() +
            "_" + Model::discretizationName() +
            (this->model().enableGridAdaptation()?"_adaptive":"");
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, InitialWaterSaturation,
                             "The initial saturation in the domain [] of the "
                             "wetting phase");
    }

    /*!
     * \copydoc FvBaseProblem::finishInit()
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 3e-6;

        temperature_ = 273.15 + 20; // -> 20°C

        FluidSystem::init();

        // parameters for the Van Genuchten law of the main imbibition
        // and the main drainage curves.
        micParams_.setVgAlpha(0.0037);
        micParams_.setVgN(4.7);
        micParams_.finalize();

        mdcParams_.setVgAlpha(0.0037);
        mdcParams_.setVgN(4.7);
        mdcParams_.finalize();

        // initialize the material parameter objects of the individual
        // finite volumes, resize will resize the container to the number of elements
        materialParams_.resize();

        for (auto it = materialParams_.begin(),
                 end = materialParams_.end(); it != end; ++it ) {
            std::shared_ptr< MaterialLawParams >& materialParams = *it ;
            if( ! materialParams )
            {
                materialParams.reset( new MaterialLawParams() );
                materialParams->setMicParams(&micParams_);
                materialParams->setMdcParams(&mdcParams_);
                materialParams->setSwr(0.0);
                materialParams->setSnr(0.1);
                materialParams->finalize();
                ParkerLenhard::reset(*materialParams);
            }
        }

        K_ = this->toDimMatrix_(4.6e-10);

        setupInitialFluidState_();
    }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        // checkConservativeness() does not include the effect of constraints, so we
        // disable it for this problem...
        //this->model().checkConservativeness();

        // Calculate storage terms
        EqVector storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage: " << storage << std::endl << std::flush;
        }
#endif // NDEBUG

        // update the history of the hysteresis law
        ElementContext elemCtx(this->simulator());

        auto elemIt = this->gridView().template begin<0>();
        const auto& elemEndIt = this->gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            elemCtx.updateAll( elem );
            size_t numDofs = elemCtx.numDof(/*timeIdx=*/0);
            for (unsigned scvIdx = 0; scvIdx < numDofs; ++scvIdx)
            {
                MaterialLawParams& materialParam = materialLawParams( elemCtx, scvIdx, /*timeIdx=*/0 );
                const auto& fs = elemCtx.intensiveQuantities(scvIdx, /*timeIdx=*/0).fluidState();
                ParkerLenhard::update(materialParam, fs);
            }
        }
    }

    //! \}

    /*!
     * \name Soil parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& /*context*/, unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const
    { return temperature_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& /*context*/, unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const
    { return K_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& /*context*/, unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const
    { return 0.4; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    MaterialLawParams& materialLawParams(const Context& context,
                                         unsigned spaceIdx, unsigned timeIdx)
    {
        const auto& entity = context.stencil(timeIdx).entity(spaceIdx);
        assert(materialParams_[entity]);
        return *materialParams_[entity];
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        const auto& entity = context.stencil(timeIdx).entity( spaceIdx );
        assert(materialParams_[entity]);
        return *materialParams_[entity];
    }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector& values, const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos) || onRightBoundary_(pos) || onLowerBoundary_(pos))
            values.setNoFlow();
        else {
            assert(onUpperBoundary_(pos));

            values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidState_);
        }

        // override the value for the liquid phase by forced
        // imbibition of water on inlet boundary segments
        if (onInlet_(pos)) {
            values[contiWettingEqIdx] = -0.001; // [kg/(m^2 s)]
        }
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
    void initial(PrimaryVariables& values, const Context& /*context*/, unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const
    {
        // assign the primary variables
        values.assignNaive(initialFluidState_);
    }

    /*!
     * \copydoc FvBaseProblem::constraints
     */
    template <class Context>
    void constraints(Constraints& constraints, const Context& context,
                     unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (onUpperBoundary_(pos) && !onInlet_(pos)) {
            constraints.setActive(true);
            constraints.assignNaive(initialFluidState_);
        }
        else if (onLowerBoundary_(pos)) {
            constraints.setActive(true);
            constraints.assignNaive(initialFluidState_);
        }
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector& rate, const Context& /*context*/,
                unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const
    { rate = Scalar(0.0); }
    //! \}

private:
    bool onLeftBoundary_(const GlobalPosition& pos) const
    { return pos[0] < this->boundingBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition& pos) const
    { return pos[1] < this->boundingBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition& pos) const
    { return pos[1] > this->boundingBoxMax()[1] - eps_; }

    bool onInlet_(const GlobalPosition& pos) const
    {
        Scalar width = this->boundingBoxMax()[0] - this->boundingBoxMin()[0];
        Scalar lambda = (this->boundingBoxMax()[0] - pos[0]) / width;

        if (!onUpperBoundary_(pos))
            return false;

        Scalar xInject[] = { 0.25, 0.75 };
        Scalar injectLen[] = { 0.1, 0.1 };
        for (unsigned i = 0; i < sizeof(xInject) / sizeof(Scalar); ++i) {
            if (xInject[i] - injectLen[i] / 2 < lambda
                && lambda < xInject[i] + injectLen[i] / 2)
                return true;
        }
        return false;
    }

    void setupInitialFluidState_()
    {
        auto& fs = initialFluidState_;
        fs.setPressure(wettingPhaseIdx, /*pressure=*/1e5);

        Scalar Sw = EWOMS_GET_PARAM(TypeTag, Scalar, InitialWaterSaturation);
        fs.setSaturation(wettingPhaseIdx, Sw);
        fs.setSaturation(nonWettingPhaseIdx, 1 - Sw);

        fs.setTemperature(temperature_);

        // set the absolute pressures
        Scalar pn = 1e5;
        fs.setPressure(nonWettingPhaseIdx, pn);
        fs.setPressure(wettingPhaseIdx, pn);

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updateAll(fs);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            fs.setDensity(phaseIdx, FluidSystem::density(fs, paramCache, phaseIdx));
            fs.setViscosity(phaseIdx, FluidSystem::viscosity(fs, paramCache, phaseIdx));
        }

    }

    DimMatrix K_;

    typename MaterialLawParams::VanGenuchtenParams micParams_;
    typename MaterialLawParams::VanGenuchtenParams mdcParams_;

    MaterialLawParamsContainer materialParams_;

    Ewoms::ImmiscibleFluidState<Scalar, FluidSystem> initialFluidState_;

    Scalar temperature_;
    Scalar eps_;
};

} // namespace Ewoms

#endif
