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
 * \brief Contains the classes required to extend the black-oil model to include the effects of foam.
 */
#ifndef EWOMS_BLACK_OIL_FOAM_MODULE_HH
#define EWOMS_BLACK_OIL_FOAM_MODULE_HH

#include "blackoilproperties.hh"
//#include <ewoms/numerics/io/vtkblackoilfoammodule.hh>
#include <ewoms/numerics/common/quantitycallbacks.hh>

#include <ewoms/common/tabulated1dfunction.hh>
//#include <ewoms/common/intervaltabulated2dfunction.hh>

#if HAVE_ECL_INPUT
#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/foamadstable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/foammobtable.hh>
#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/opmlog/opmlog.hh>
#endif

#include <ewoms/common/valgrind.hh>
#include <ewoms/common/unused.hh>
#include <ewoms/common/exceptions.hh>

#include <dune/common/fvector.hh>

#include <string>
#include <math.h>

namespace Ewoms {
/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model to include the effects of foam.
 */
template <class TypeTag, bool enableFoamV = GET_PROP_VALUE(TypeTag, EnableFoam)>
class BlackOilFoamModule
{
    typedef GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef GET_PROP_TYPE(TypeTag, Model) Model;
    typedef GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef Ewoms::MathToolbox<Evaluation> Toolbox;

    typedef typename Ewoms::Tabulated1DFunction<Scalar> TabulatedFunction;

    static constexpr unsigned foamConcentrationIdx = Indices::foamConcentrationIdx;
    static constexpr unsigned contiFoamEqIdx = Indices::contiFoamEqIdx;
    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;

    static constexpr unsigned enableFoam = enableFoamV;
    static constexpr bool enableVtkOutput = GET_PROP_VALUE(TypeTag, EnableVtkOutput);

    static constexpr unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);

public:
    // a struct containing constants to calculate change to relative permeability,
    // based on model (1-9) in Table 1 of
    // Kun Ma, Guangwei Ren, Khalid Mateen, Danielle Morel, and Philippe Cordelier:
    // "Modeling techniques for foam flow in porous media", SPE Journal, 20(03):453–470, jun 2015.
    // The constants are provided by various deck keywords as shown in the comments below.
    struct FoamCoefficients {
        Scalar fm_min = 1e-20;   // FOAMFSC
        Scalar fm_mob = 1.0;     // FOAMFRM

        Scalar fm_surf = 1.0;    // FOAMFSC
        Scalar ep_surf = 1.0;    // FOAMFSC

        Scalar fm_oil = 1.0;     // FOAMFSO
        Scalar fl_oil = 0.0;     // FOAMFSO
        Scalar ep_oil = 0.0;     // FOAMFSO

        Scalar fm_cap = 1.0;     // FOAMFCN
        Scalar ep_cap = 0.0;     // FOAMFCN

        Scalar fm_dry = 1.0;     // FOAMFSW
        Scalar ep_dry = 0.0;     // FOAMFSW
    };

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize all internal data structures needed by the foam module
     */
    static void initFromEclState(const Ewoms::EclipseState& eclState)
    {
        // some sanity checks: if foam is enabled, the FOAM keyword must be
        // present, if foam is disabled the keyword must not be present.
        if (enableFoam && !eclState.runspec().phases().active(Phase::FOAM)) {
            throw std::runtime_error("Non-trivial foam treatment requested at compile time, but "
                                     "the deck does not contain the FOAM keyword");
        }
        else if (!enableFoam && eclState.runspec().phases().active(Phase::FOAM)) {
            throw std::runtime_error("Foam treatment disabled at compile time, but the deck "
                                     "contains the FOAM keyword");
        }

        if (!eclState.runspec().phases().active(Phase::FOAM)) {
            return; // foam treatment is supposed to be disabled
        }

        // Check that only implemented options are used.
        // We only support the default values of FOAMOPTS (GAS, TAB).
        if (eclState.getInitConfig().getFoamConfig().getTransportPhase() != Phase::GAS) {
            throw std::runtime_error("In FOAMOPTS, only GAS is allowed for the transport phase.");
        }
        if (eclState.getInitConfig().getFoamConfig().getMobilityModel() != FoamConfig::MobilityModel::TAB) {
            throw std::runtime_error("In FOAMOPTS, only TAB is allowed for the gas mobility factor reduction model.");
        }

        const auto& tableManager = eclState.getTableManager();
        const unsigned int numSatRegions = tableManager.getTabdims().getNumSatTables();
        setNumSatRegions(numSatRegions);
        const unsigned int numPvtRegions = tableManager.getTabdims().getNumPVTTables();
        setNumPvtRegions(numPvtRegions);

        // Get and check FOAMROCK data.
        const Ewoms::FoamConfig& foamConf = eclState.getInitConfig().getFoamConfig();
        if (numSatRegions != foamConf.size()) {
            throw std::runtime_error("Inconsistent sizes, number of saturation regions differ from the number of elements "
                                     "in FoamConfig, which typically corresponds to the number of records in FOAMROCK.");
        }

        // Get and check FOAMADS data.
        const auto& foamadsTables = tableManager.getFoamadsTables();
        if (foamadsTables.empty()) {
            throw std::runtime_error("FOAMADS must be specified in FOAM runs");
        }
        if (numSatRegions != foamadsTables.size()) {
            throw std::runtime_error("Inconsistent sizes, number of saturation regions differ from the "
                                     "number of FOAMADS tables.");
        }

        // Set data that vary with saturation region.
        for (std::size_t satReg = 0; satReg < numSatRegions; ++satReg) {
            const auto& rec = foamConf.getRecord(satReg);
            foamCoefficients_[satReg] = FoamCoefficients();
            foamCoefficients_[satReg].fm_min = rec.minimumSurfactantConcentration();
            foamCoefficients_[satReg].fm_surf = rec.referenceSurfactantConcentration();
            foamCoefficients_[satReg].ep_surf = rec.exponent();
            foamRockDensity_[satReg] = rec.rockDensity();
            foamAllowDesorption_[satReg] = rec.allowDesorption();
            const auto& foamadsTable = foamadsTables.template getTable<Ewoms::FoamadsTable>(satReg);
            const auto& conc = foamadsTable.getFoamConcentrationColumn();
            const auto& ads = foamadsTable.getAdsorbedFoamColumn();
            adsorbedFoamTable_[satReg].setXYContainers(conc, ads);
        }

        // Get and check FOAMMOB data.
        const auto& foammobTables = tableManager.getFoammobTables();
        if (foammobTables.empty()) {
            // When in the future adding support for the functional
            // model, FOAMMOB will not be required anymore (functional
            // family of keywords can be used instead, FOAMFSC etc.).
            throw std::runtime_error("FOAMMOB must be specified in FOAM runs");
        }
        if (numPvtRegions != foammobTables.size()) {
            throw std::runtime_error("Inconsistent sizes, number of PVT regions differ from the "
                                     "number of FOAMMOB tables.");
        }

        // Set data that vary with PVT region.
        for (std::size_t pvtReg = 0; pvtReg < numPvtRegions; ++pvtReg) {
            const auto& foammobTable = foammobTables.template getTable<Ewoms::FoammobTable>(pvtReg);
            const auto& conc = foammobTable.getFoamConcentrationColumn();
            const auto& mobMult = foammobTable.getMobilityMultiplierColumn();
            gasMobilityMultiplierTable_[pvtReg].setXYContainers(conc, mobMult);
        }
    }
#endif

    /*!
     * \brief Specify the number of saturation regions.
     */
    static void setNumSatRegions(unsigned numRegions)
    {
        foamCoefficients_.resize(numRegions);
        foamRockDensity_.resize(numRegions);
        foamAllowDesorption_.resize(numRegions);
        adsorbedFoamTable_.resize(numRegions);
    }

    /*!
     * \brief Specify the number of PVT regions.
     */
    static void setNumPvtRegions(unsigned numRegions)
    {
        gasMobilityMultiplierTable_.resize(numRegions);
    }

    /*!
     * \brief Register all run-time parameters for the black-oil foam module.
     */
    static void registerParameters()
    {
        if (!enableFoam)
            // foam has been disabled at compile time
            return;

        //Ewoms::VtkBlackOilFoamModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all foam specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model EWOMS_UNUSED,
                                      Simulator& simulator EWOMS_UNUSED)
    {
        if (!enableFoam)
            // foam have been disabled at compile time
            return;

#if HAVE_EWOMS_ECLIO
        if (enableVtkOutput) {
            Ewoms::OpmLog::warning("VTK output requested, currently unsupported by the foam module.");
        }
#endif
        //model.addOutputModule(new Ewoms::VtkBlackOilFoamModule<TypeTag>(simulator));
    }

    static bool primaryVarApplies(unsigned pvIdx)
    {
        if (!enableFoam) {
            return false;
        } else {
            return pvIdx == foamConcentrationIdx;
        }
    }

    static std::string primaryVarName(unsigned pvIdx EWOMS_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));
        return "foam_concentration";
    }

    static Scalar primaryVarWeight(unsigned pvIdx EWOMS_OPTIM_UNUSED)
    {
       assert(primaryVarApplies(pvIdx));

       // TODO: it may be beneficial to chose this differently.
       return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if (!enableFoam)
            return false;

        return eqIdx == contiFoamEqIdx;

    }

    static std::string eqName(unsigned eqIdx EWOMS_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        return "conti^foam";
    }

    static Scalar eqWeight(unsigned eqIdx EWOMS_OPTIM_UNUSED)
    {
       assert(eqApplies(eqIdx));

       // TODO: it may be beneficial to chose this differently.
       return static_cast<Scalar>(1.0);
    }

    // must be called after water storage is computed
    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if (!enableFoam)
            return;

        const auto& fs = intQuants.fluidState();

        LhsEval surfaceVolumeFreeGas =
            Toolbox::template decay<LhsEval>(fs.saturation(gasPhaseIdx))
            * Toolbox::template decay<LhsEval>(fs.invB(gasPhaseIdx))
            * Toolbox::template decay<LhsEval>(intQuants.porosity());

        // Avoid singular matrix if no gas is present.
        surfaceVolumeFreeGas = Ewoms::max(surfaceVolumeFreeGas, 1e-10);

        // Foam/surfactant in gas phase.
        const LhsEval gasFoam = surfaceVolumeFreeGas
            * Toolbox::template decay<LhsEval>(intQuants.foamConcentration());

        // Adsorbed foam/surfactant.
        const LhsEval adsorbedFoam =
            Toolbox::template decay<LhsEval>(1.0 - intQuants.porosity())
            * Toolbox::template decay<LhsEval>(intQuants.foamRockDensity())
            * Toolbox::template decay<LhsEval>(intQuants.foamAdsorbed());

        LhsEval accumulationFoam = gasFoam + adsorbedFoam;
        storage[contiFoamEqIdx] += accumulationFoam;
    }

    static void computeFlux(RateVector& flux,
                            const ElementContext& elemCtx,
                            unsigned scvfIdx,
                            unsigned timeIdx)

    {
        if (!enableFoam) {
            return;
        }

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        const unsigned upIdx = extQuants.upstreamIndex(FluidSystem::gasPhaseIdx);
        const unsigned inIdx = extQuants.interiorIndex();
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);

        // The effect of the gas mobility reduction factor is
        // incorporated in the mobility, so the oil (if vaporized oil
        // is active) and gas fluxes do not need modification here.
        if (upIdx == inIdx) {
            flux[contiFoamEqIdx] =
                extQuants.volumeFlux(gasPhaseIdx)
                *up.fluidState().invB(gasPhaseIdx)
                *up.foamConcentration();
        } else {
            flux[contiFoamEqIdx] =
                extQuants.volumeFlux(gasPhaseIdx)
                *Ewoms::decay<Scalar>(up.fluidState().invB(gasPhaseIdx))
                *Ewoms::decay<Scalar>(up.foamConcentration());
        }

    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables& oldPv EWOMS_UNUSED,
                                     const EqVector& delta EWOMS_UNUSED)
    {
        // do not consider the change of foam primary variables for convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    template <class DofEntity>
    static void serializeEntity(const Model& model, std::ostream& outstream, const DofEntity& dof)
    {
        if (!enableFoam)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);
        const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
        outstream << priVars[foamConcentrationIdx];
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
        if (!enableFoam)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);
        PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
        PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

        instream >> priVars0[foamConcentrationIdx];

        // set the primary variables for the beginning of the current time step.
        priVars1[foamConcentrationIdx] = priVars0[foamConcentrationIdx];
    }

    static const Scalar foamRockDensity(const ElementContext& elemCtx,
                                        unsigned scvIdx,
                                        unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return foamRockDensity_[satnumRegionIdx];
    }

    static bool foamAllowDesorption(const ElementContext& elemCtx,
                                    unsigned scvIdx,
                                    unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return foamAllowDesorption_[satnumRegionIdx];
    }

    static const TabulatedFunction& adsorbedFoamTable(const ElementContext& elemCtx,
                                                      unsigned scvIdx,
                                                      unsigned timeIdx)
    {
       unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
       return adsorbedFoamTable_[satnumRegionIdx];
    }

    static const TabulatedFunction& gasMobilityMultiplierTable(const ElementContext& elemCtx,
                                                               unsigned scvIdx,
                                                               unsigned timeIdx)
    {
       unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
       return gasMobilityMultiplierTable_[pvtnumRegionIdx];
    }

    static const FoamCoefficients& foamCoefficients(const ElementContext& elemCtx,
                                                    const unsigned scvIdx,
                                                    const unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return foamCoefficients_[satnumRegionIdx];
    }

private:
    static std::vector<Scalar> foamRockDensity_;
    static std::vector<bool> foamAllowDesorption_;
    static std::vector<FoamCoefficients> foamCoefficients_;
    static std::vector<TabulatedFunction> adsorbedFoamTable_;
    static std::vector<TabulatedFunction> gasMobilityMultiplierTable_;
};

template <class TypeTag, bool enableFoam>
std::vector<typename BlackOilFoamModule<TypeTag, enableFoam>::Scalar>
BlackOilFoamModule<TypeTag, enableFoam>::foamRockDensity_;

template <class TypeTag, bool enableFoam>
std::vector<bool>
BlackOilFoamModule<TypeTag, enableFoam>::foamAllowDesorption_;

template <class TypeTag, bool enableFoam>
std::vector<typename BlackOilFoamModule<TypeTag, enableFoam>::FoamCoefficients>
BlackOilFoamModule<TypeTag, enableFoam>::foamCoefficients_;

template <class TypeTag, bool enableFoam>
std::vector<typename BlackOilFoamModule<TypeTag, enableFoam>::TabulatedFunction>
BlackOilFoamModule<TypeTag, enableFoam>::adsorbedFoamTable_;

template <class TypeTag, bool enableFoam>
std::vector<typename BlackOilFoamModule<TypeTag, enableFoam>::TabulatedFunction>
BlackOilFoamModule<TypeTag, enableFoam>::gasMobilityMultiplierTable_;

/*!
 * \ingroup BlackOil
 * \class Ewoms::BlackOilFoamIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        polymers extension of the black-oil model.
 */
template <class TypeTag, bool enableFoam = GET_PROP_VALUE(TypeTag, EnableFoam)>
class BlackOilFoamIntensiveQuantities
{
    typedef GET_PROP_TYPE(TypeTag, IntensiveQuantities) Implementation;

    typedef GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef BlackOilFoamModule<TypeTag> FoamModule;

    static constexpr int foamConcentrationIdx = Indices::foamConcentrationIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr unsigned oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;

public:

    /*!
     * \brief Update the intensive properties needed to handle polymers from the
     *        primary variables
     *
     */
    void foamPropertiesUpdate_(const ElementContext& elemCtx,
                               unsigned dofIdx,
                               unsigned timeIdx)
    {
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        foamConcentration_ = priVars.makeEvaluation(foamConcentrationIdx, timeIdx);
        const auto& fs = asImp_().fluidState_;

        // Compute gas mobility reduction factor
        Evaluation mobilityReductionFactor = 1.0;
        if (false) {
            // The functional model is used.
            // TODO: allow this model.
            // In order to do this we must allow transport to be in the water phase, not just the gas phase.
            const auto& foamCoefficients = FoamModule::foamCoefficients(elemCtx, dofIdx, timeIdx);

            const Scalar fm_mob = foamCoefficients.fm_mob;

            const Scalar fm_surf = foamCoefficients.fm_surf;
            const Scalar ep_surf = foamCoefficients.ep_surf;

            const Scalar fm_oil = foamCoefficients.fm_oil;
            const Scalar fl_oil = foamCoefficients.fl_oil;
            const Scalar ep_oil = foamCoefficients.ep_oil;

            const Scalar fm_dry = foamCoefficients.fm_dry;
            const Scalar ep_dry = foamCoefficients.ep_dry;

            const Scalar fm_cap = foamCoefficients.fm_cap;
            const Scalar ep_cap = foamCoefficients.ep_cap;

            const Evaluation C_surf = foamConcentration_;
            const Evaluation Ca = 1e10; // TODO: replace with proper capillary number.
            const Evaluation S_o = fs.saturation(oilPhaseIdx);
            const Evaluation S_w = fs.saturation(waterPhaseIdx);

            Evaluation F1 = pow(C_surf/fm_surf, ep_surf);
            Evaluation F2 = pow((fm_oil-S_o)/(fm_oil-fl_oil), ep_oil);
            Evaluation F3 = pow(fm_cap/Ca, ep_cap);
            Evaluation F7 = 0.5 + atan(ep_dry*(S_w-fm_dry))/M_PI;

            mobilityReductionFactor = 1./(1. + fm_mob*F1*F2*F3*F7);
        } else {
            // The tabular model is used.
            // Note that the current implementation only includes the effect of foam concentration (FOAMMOB),
            // and not the optional pressure dependence (FOAMMOBP) or shear dependence (FOAMMOBS).
            const auto& gasMobilityMultiplier = FoamModule::gasMobilityMultiplierTable(elemCtx, dofIdx, timeIdx);
            mobilityReductionFactor = gasMobilityMultiplier.eval(foamConcentration_, /* extrapolate = */ true);
        }

        // adjust gas mobility
        asImp_().mobility_[gasPhaseIdx] *= mobilityReductionFactor;

        foamRockDensity_ = FoamModule::foamRockDensity(elemCtx, dofIdx, timeIdx);

        const auto& adsorbedFoamTable = FoamModule::adsorbedFoamTable(elemCtx, dofIdx, timeIdx);
        foamAdsorbed_ = adsorbedFoamTable.eval(foamConcentration_, /*extrapolate=*/true);
        if (!FoamModule::foamAllowDesorption(elemCtx, dofIdx, timeIdx)) {
            throw std::runtime_error("Foam module does not support the 'no desorption' option.");
        }
    }

    const Evaluation& foamConcentration() const
    { return foamConcentration_; }

    Scalar foamRockDensity() const
    { return foamRockDensity_; }

    const Evaluation& foamAdsorbed() const
    { return foamAdsorbed_; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation foamConcentration_;
    Scalar foamRockDensity_;
    Evaluation foamAdsorbed_;
};

template <class TypeTag>
class BlackOilFoamIntensiveQuantities<TypeTag, false>
{
    typedef GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    void foamPropertiesUpdate_(const ElementContext& elemCtx EWOMS_UNUSED,
                                  unsigned scvIdx EWOMS_UNUSED,
                                  unsigned timeIdx EWOMS_UNUSED)
    { }

    const Evaluation& foamConcentration() const
    { throw std::runtime_error("foamConcentration() called but foam is disabled"); }

    Scalar foamRockDensity() const
    { throw std::runtime_error("foamRockDensity() called but foam is disabled"); }

    Scalar foamAdsorbed() const
    { throw std::runtime_error("foamAdsorbed() called but foam is disabled"); }
};

} // namespace Ewoms

#endif
