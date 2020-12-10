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
 * \brief Contains the classes required to extend black-oil solvent model proposed by Sandve, Sævareid and Aavatsmark.
 *
 *  For details, refer:
 *  [*] T.H. Sandve, O. Sævareid and I. Aavatsmark: “Improved Extended Blackoil Formulation
 *  for CO2 EOR Simulations.” in ECMOR XVII – The 17th European Conference on the
 *  Mathematics of Oil Recovery,  September 2020.
 */
#ifndef EWOMS_BLACK_OIL_SSA_SOLVENT_HH
#define EWOMS_BLACK_OIL_SSA_SOLVENT_HH

#include "blackoilproperties.hh"

//#include <ewoms/numerics/io/vtkblackoilssasolventmodule.hh> //todo: missing ...

#include <ewoms/numerics/common/quantitycallbacks.hh>

#include <ewoms/common/tabulated1dfunction.hh>
#include <ewoms/common/uniformxtabulated2dfunction.hh>

#if HAVE_ECL_INPUT
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/ssfntable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/sof2table.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/msfntable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/pmisctable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/misctable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/sorwmistable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/sgcwmistable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tlpmixpatable.hh>
#endif

#include <ewoms/common/valgrind.hh>
#include <ewoms/common/unused.hh>
#include <ewoms/common/exceptions.hh>

#include <dune/common/fvector.hh>

#include <string>

namespace Ewoms {
/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model.
 */
template <class TypeTag, bool enableSsaSolventV = GET_PROP_VALUE(TypeTag, EnableSsaSolvent)>
class BlackOilSsaSolventModule
{
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);
    using PrimaryVariables = GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using IntensiveQuantities = GET_PROP_TYPE(TypeTag, IntensiveQuantities);
    using ExtensiveQuantities = GET_PROP_TYPE(TypeTag, ExtensiveQuantities);
    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
    using Model = GET_PROP_TYPE(TypeTag, Model);
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using EqVector = GET_PROP_TYPE(TypeTag, EqVector);
    using RateVector = GET_PROP_TYPE(TypeTag, RateVector);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);

    typedef Ewoms::MathToolbox<Evaluation> Toolbox;

    typedef typename Ewoms::Tabulated1DFunction<Scalar> TabulatedFunction;
    typedef typename Ewoms::UniformXTabulated2DFunction<Scalar> Tabulated2DFunction;

    static constexpr unsigned zFractionIdx = Indices::zFractionIdx;
    static constexpr unsigned contiZfracEqIdx = Indices::contiZfracEqIdx;
    static constexpr unsigned enableSsaSolvent = enableSsaSolventV;
    static constexpr unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr unsigned numPhases = FluidSystem::numPhases;
    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr unsigned oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr bool blackoilConserveSurfaceVolume = GET_PROP_VALUE(TypeTag, BlackoilConserveSurfaceVolume);

public:
#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize all internal data structures needed by the solvent module
     */
    static void initFromEclState(const Ewoms::EclipseState& eclState)
    {
        // some sanity checks: if extended BO is enabled, the PVTSOL keyword must be
        // present, if extended BO is disabled the keyword must not be present.
        if (enableSsaSolvent && !eclState.runspec().phases().active(Phase::ZFRACTION))
            throw std::runtime_error("Extended black oil treatment requested at compile "
                                     "time, but the deck does not contain the PVTSOL keyword");
        else if (!enableSsaSolvent && eclState.runspec().phases().active(Phase::ZFRACTION))
            throw std::runtime_error("Extended black oil treatment disabled at compile time, but the deck "
                                     "contains the PVTSOL keyword");

        if (!eclState.runspec().phases().active(Phase::ZFRACTION))
            return; // solvent treatment is supposed to be disabled

        // pvt properties from kw PVTSOL:

        const auto& tableManager = eclState.getTableManager();
        const auto& pvtsolTables = tableManager.getPvtsolTables();

        size_t numPvtRegions = pvtsolTables.size();

        BO_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        BG_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        RS_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        RV_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        X_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        Y_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        VISCO_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        VISCG_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});

        PBUB_RS_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        PBUB_RV_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});

        zLim_.resize(numPvtRegions);

        const bool extractCmpFromPvt = true; //<false>: Default values used in [*]
        oilCmp_.resize(numPvtRegions);
        gasCmp_.resize(numPvtRegions);

        for (unsigned regionIdx = 0; regionIdx < numPvtRegions; ++ regionIdx) {
          const auto& pvtsolTable = pvtsolTables[regionIdx];

          const auto& saturatedTable = pvtsolTable.getSaturatedTable();
          assert(saturatedTable.numRows() > 1);

          std::vector<Scalar> oilCmp(saturatedTable.numRows(), -4.0e-9); //Default values used in [*]
          std::vector<Scalar> gasCmp(saturatedTable.numRows(), -0.08);   //-------------"-------------
          zLim_[regionIdx] = 0.7;                                        //-------------"-------------
          std::vector<Scalar> zArg(saturatedTable.numRows(), 0.0);

          for (unsigned outerIdx = 0; outerIdx < saturatedTable.numRows(); ++ outerIdx) {
            Scalar ZCO2 = saturatedTable.get("ZCO2", outerIdx);

            zArg[outerIdx] = ZCO2;

            BO_[regionIdx].appendXPos(ZCO2);
            BG_[regionIdx].appendXPos(ZCO2);

            RS_[regionIdx].appendXPos(ZCO2);
            RV_[regionIdx].appendXPos(ZCO2);

            X_[regionIdx].appendXPos(ZCO2);
            Y_[regionIdx].appendXPos(ZCO2);

            VISCO_[regionIdx].appendXPos(ZCO2);
            VISCG_[regionIdx].appendXPos(ZCO2);

            PBUB_RS_[regionIdx].appendXPos(ZCO2);
            PBUB_RV_[regionIdx].appendXPos(ZCO2);

            const auto& underSaturatedTable = pvtsolTable.getUnderSaturatedTable(outerIdx);
            size_t numRows = underSaturatedTable.numRows();

            Scalar bo0=0.0;
            Scalar po0=0.0;
            for (unsigned innerIdx = 0; innerIdx < numRows; ++ innerIdx) {
              Scalar po = underSaturatedTable.get("P", innerIdx);
              Scalar bo = underSaturatedTable.get("B_O", innerIdx);
              Scalar bg = underSaturatedTable.get("B_G", innerIdx);
              Scalar rs = underSaturatedTable.get("RS", innerIdx)+innerIdx*1.0e-10;
              Scalar rv = underSaturatedTable.get("RV", innerIdx)+innerIdx*1.0e-10;
              Scalar xv = underSaturatedTable.get("XVOL", innerIdx);
              Scalar yv = underSaturatedTable.get("YVOL", innerIdx);
              Scalar mo = underSaturatedTable.get("MU_O", innerIdx);
              Scalar mg = underSaturatedTable.get("MU_G", innerIdx);

              if (bo0 > bo) { // This is undersaturated oil-phase for ZCO2 <= zLim ...
                              // Here we assume tabulated bo to decay beyond boiling point
                  if (extractCmpFromPvt) {
                      Scalar cmpFactor = (bo-bo0)/(po-po0);
                      oilCmp[outerIdx] = cmpFactor;
                      zLim_[regionIdx] = ZCO2;
                      //std::cout << "### cmpFactorOil: " << cmpFactor << "  zLim: " << zLim_[regionIdx] << std::endl;
                  }
                  break;
              } else if (bo0 == bo) { // This is undersaturated gas-phase for ZCO2 > zLim ...
                                      // Here we assume tabulated bo to be constant extrapolated beyond dew point
                  if (innerIdx+1 < numRows && ZCO2<1.0 && extractCmpFromPvt) {
                    Scalar rvNxt = underSaturatedTable.get("RV", innerIdx+1)+innerIdx*1.0e-10;
                    Scalar bgNxt = underSaturatedTable.get("B_G", innerIdx+1);
                    Scalar cmpFactor = (bgNxt-bg)/(rvNxt-rv);
                    gasCmp[outerIdx] = cmpFactor;
                    //std::cout << "### cmpFactorGas: " << cmpFactor << "  zLim: " << zLim_[regionIdx] << std::endl;
                  }

                  BO_[regionIdx].appendSamplePoint(outerIdx,po,bo);
                  BG_[regionIdx].appendSamplePoint(outerIdx,po,bg);
                  RS_[regionIdx].appendSamplePoint(outerIdx,po,rs);
                  RV_[regionIdx].appendSamplePoint(outerIdx,po,rv);
                  X_[regionIdx].appendSamplePoint(outerIdx,po,xv);
                  Y_[regionIdx].appendSamplePoint(outerIdx,po,yv);
                  VISCO_[regionIdx].appendSamplePoint(outerIdx,po,mo);
                  VISCG_[regionIdx].appendSamplePoint(outerIdx,po,mg);
                  break;
              }

              bo0=bo;
              po0=po;

              BO_[regionIdx].appendSamplePoint(outerIdx,po,bo);
              BG_[regionIdx].appendSamplePoint(outerIdx,po,bg);

              RS_[regionIdx].appendSamplePoint(outerIdx,po,rs);
              RV_[regionIdx].appendSamplePoint(outerIdx,po,rv);

              X_[regionIdx].appendSamplePoint(outerIdx,po,xv);
              Y_[regionIdx].appendSamplePoint(outerIdx,po,yv);

              VISCO_[regionIdx].appendSamplePoint(outerIdx,po,mo);
              VISCG_[regionIdx].appendSamplePoint(outerIdx,po,mg);

              // rs,rv -> pressure
              PBUB_RS_[regionIdx].appendSamplePoint(outerIdx, rs, po);
              PBUB_RV_[regionIdx].appendSamplePoint(outerIdx, rv, po);

            }
          }
          oilCmp_[regionIdx].setXYContainers(zArg, oilCmp, /*sortInput=*/false);
          gasCmp_[regionIdx].setXYContainers(zArg, gasCmp, /*sortInput=*/false);
        }

        // Reference density for pure z-component taken from kw SDENSITY
        const auto& sdensityTables = eclState.getTableManager().getSolventDensityTables();
        if (sdensityTables.size() == numPvtRegions) {
           zReferenceDensity_.resize(numPvtRegions);
           for (unsigned regionIdx = 0; regionIdx < numPvtRegions; ++ regionIdx) {
             Scalar rhoRefS = sdensityTables[regionIdx].getSolventDensityColumn().front();
             zReferenceDensity_[regionIdx]=rhoRefS;
           }
        }
        else
           throw std::runtime_error("SsaSolvent:  kw SDENSITY is missing or not aligned with NTPVT\n");
    }
#endif

    /*!
     * \brief Register all run-time parameters for the black-oil solvent module.
     */
    static void registerParameters()
    {
        if (!enableSsaSolvent)
            // SSA solvents have disabled at compile time
            return;

        //Ewoms::VtkBlackOilSsaSolventModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all solvent specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if (!enableSsaSolvent)
            // SSA solvents have disabled at compile time
            return;

        //model.addOutputModule(new Ewoms::VtkBlackOilSsaSolventModule<TypeTag>(simulator));
    }

    static bool primaryVarApplies(unsigned pvIdx)
    {
        if (!enableSsaSolvent)
            // SSA solvents have disabled at compile time
            return false;

        return pvIdx == zFractionIdx;
    }

    static std::string primaryVarName(unsigned pvIdx EWOMS_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        return "z_fraction";
    }

    static Scalar primaryVarWeight(unsigned pvIdx EWOMS_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if (!enableSsaSolvent)
            return false;

        return eqIdx == contiZfracEqIdx;
    }

    static std::string eqName(unsigned eqIdx EWOMS_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        return "conti^solvent";
    }

    static Scalar eqWeight(unsigned eqIdx EWOMS_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if (!enableSsaSolvent)
            return;

        if (blackoilConserveSurfaceVolume) {
            storage[contiZfracEqIdx] =
                    Toolbox::template decay<LhsEval>(intQuants.porosity())
                    * Toolbox::template decay<LhsEval>(intQuants.yVolume())
                    * Toolbox::template decay<LhsEval>(intQuants.fluidState().saturation(gasPhaseIdx))
                    * Toolbox::template decay<LhsEval>(intQuants.fluidState().invB(gasPhaseIdx));
            if (FluidSystem::enableDissolvedGas()) { // account for dissolved z in oil phase
                storage[contiZfracEqIdx] +=
                    Toolbox::template decay<LhsEval>(intQuants.porosity())
                    * Toolbox::template decay<LhsEval>(intQuants.xVolume())
                    * Toolbox::template decay<LhsEval>(intQuants.fluidState().Rs())
                    * Toolbox::template decay<LhsEval>(intQuants.fluidState().saturation(oilPhaseIdx))
                    * Toolbox::template decay<LhsEval>(intQuants.fluidState().invB(oilPhaseIdx));
            }
            // Reg. terms: Preliminary attempt to avoid singular behaviour when solvent is invading a pure water
            //             region. Results seems insensitive to the weighting factor.
            // TODO: Further investigations ...
            const Scalar regWghtFactor = 1.0e-6;
            storage[contiZfracEqIdx] += regWghtFactor*(1.0-Toolbox::template decay<LhsEval>(intQuants.zFraction()))
                                      + regWghtFactor*Toolbox::template decay<LhsEval>(intQuants.porosity())
                                                   * Toolbox::template decay<LhsEval>(intQuants.fluidState().saturation(gasPhaseIdx))
                                                   * Toolbox::template decay<LhsEval>(intQuants.fluidState().invB(gasPhaseIdx));
            storage[contiZfracEqIdx-1] += regWghtFactor*Toolbox::template decay<LhsEval>(intQuants.zFraction());
        }
        else {
            throw std::runtime_error("Only component conservation in terms of surface volumes is implemented. ");
        }

    }

    static void computeFlux(RateVector& flux,
                            const ElementContext& elemCtx,
                            unsigned scvfIdx,
                            unsigned timeIdx)

    {
        if (!enableSsaSolvent)
            return;

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        if (blackoilConserveSurfaceVolume) {
            unsigned inIdx = extQuants.interiorIndex();

            unsigned upIdxGas = static_cast<unsigned>(extQuants.upstreamIndex(gasPhaseIdx));
            const auto& upGas = elemCtx.intensiveQuantities(upIdxGas, timeIdx);
            const auto& fsGas = upGas.fluidState();
            if (upIdxGas == inIdx) {
                flux[contiZfracEqIdx] =
                    extQuants.volumeFlux(gasPhaseIdx)
                    * (upGas.yVolume())
                    * fsGas.invB(gasPhaseIdx);
            }
            else {
                flux[contiZfracEqIdx] =
                        extQuants.volumeFlux(gasPhaseIdx)
                        * (Ewoms::decay<Scalar>(upGas.yVolume()))
                        * Ewoms::decay<Scalar>(fsGas.invB(gasPhaseIdx));
            }
            if (FluidSystem::enableDissolvedGas()) { // account for dissolved z in oil phase
                unsigned upIdxOil = static_cast<unsigned>(extQuants.upstreamIndex(oilPhaseIdx));
                const auto& upOil = elemCtx.intensiveQuantities(upIdxOil, timeIdx);
                const auto& fsOil = upOil.fluidState();
                if (upIdxOil == inIdx) {
                    flux[contiZfracEqIdx] +=
                        extQuants.volumeFlux(oilPhaseIdx)
                        * upOil.xVolume()
                        * fsOil.Rs()
                        * fsOil.invB(oilPhaseIdx);
                }
                else {
                    flux[contiZfracEqIdx] +=
                        extQuants.volumeFlux(oilPhaseIdx)
                        * Ewoms::decay<Scalar>(upOil.xVolume())
                        * Ewoms::decay<Scalar>(fsOil.Rs())
                        * Ewoms::decay<Scalar>(fsOil.invB(oilPhaseIdx));
                }
            }
        }
        else {
            throw std::runtime_error("Only component conservation in terms of surface volumes is implemented. ");
        }
    }

    /*!
     * \brief Assign the solvent specific primary variables to a PrimaryVariables object
     */
    static void assignPrimaryVars(PrimaryVariables& priVars,
                                  Scalar zFraction)
    {
        if (!enableSsaSolvent)
            return;

        priVars[zFractionIdx] = zFraction;
    }

    /*!
     * \brief Do a Newton-Raphson update the primary variables of the solvents.
     */
    static void updatePrimaryVars(PrimaryVariables& newPv,
                                  const PrimaryVariables& oldPv,
                                  const EqVector& delta)
    {
        if (!enableSsaSolvent)
            return;

        // do a plain unchopped Newton update
        newPv[zFractionIdx] = oldPv[zFractionIdx] - delta[zFractionIdx];
    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables& oldPv EWOMS_UNUSED,
                                     const EqVector& delta EWOMS_UNUSED)
    {
        // do not consider consider the cange of solvent primary variables for
        // convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    /*!
     * \brief Return how much a residual is considered an error
     */
    static Scalar computeResidualError(const EqVector& resid)
    {
        // do not weight the residual of solvents when it comes to convergence
        return std::abs(Toolbox::scalarValue(resid[contiZfracEqIdx]));
    }

    template <class DofEntity>
    static void serializeEntity(const Model& model, std::ostream& outstream, const DofEntity& dof)
    {
        if (!enableSsaSolvent)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);

        const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
        outstream << priVars[zFractionIdx];
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
        if (!enableSsaSolvent)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);

        PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
        PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

        instream >> priVars0[zFractionIdx];

        // set the primary variables for the beginning of the current time step.
        priVars1 = priVars0[zFractionIdx];
    }

    template <typename Value>
    static Value xVolume(unsigned pvtRegionIdx, const Value& pressure, const Value& z) {
        return X_[pvtRegionIdx].eval(z, pressure);
    }

    template <typename Value>
    static Value yVolume(unsigned pvtRegionIdx, const Value& pressure, const Value& z) {
        return Y_[pvtRegionIdx].eval(z, pressure);
    }

    template <typename Value>
    static Value pbubRs(unsigned pvtRegionIdx, const Value& z, const Value& rs) {
        return PBUB_RS_[pvtRegionIdx].eval(z, rs);
    }

    template <typename Value>
    static Value pbubRv(unsigned pvtRegionIdx, const Value& z, const Value& rv) {
        return PBUB_RV_[pvtRegionIdx].eval(z, rv);
    }

    template <typename Value>
    static Value oilViscosity(unsigned pvtRegionIdx, const Value& pressure, const Value& z) {
        return VISCO_[pvtRegionIdx].eval(z, pressure);
    }

    template <typename Value>
    static Value gasViscosity(unsigned pvtRegionIdx, const Value& pressure, const Value& z) {
        return VISCG_[pvtRegionIdx].eval(z, pressure);
    }

    template <typename Value>
    static Value bo(unsigned pvtRegionIdx, const Value& pressure, const Value& z) {
        return BO_[pvtRegionIdx].eval(z, pressure);
    }

    template <typename Value>
    static Value bg(unsigned pvtRegionIdx, const Value& pressure, const Value& z) {
        return BG_[pvtRegionIdx].eval(z, pressure);
    }

    template <typename Value>
    static Value rs(unsigned pvtRegionIdx, const Value& pressure, const Value& z) {
        return RS_[pvtRegionIdx].eval(z, pressure);
    }

    template <typename Value>
    static Value rv(unsigned pvtRegionIdx, const Value& pressure, const Value& z) {
        return RV_[pvtRegionIdx].eval(z, pressure);
    }

    static Scalar referenceDensity(unsigned regionIdx) {
        return zReferenceDensity_[regionIdx];
    }

    static Scalar zLim(unsigned regionIdx) {
        return zLim_[regionIdx];
    }

    template <typename Value>
    static Value oilCmp(unsigned pvtRegionIdx, const Value& z) {
        return oilCmp_[pvtRegionIdx].eval(z);
    }

    template <typename Value>
    static Value gasCmp(unsigned pvtRegionIdx, const Value& z) {
        return gasCmp_[pvtRegionIdx].eval(z);
    }

private:
    static std::vector<Tabulated2DFunction> X_;
    static std::vector<Tabulated2DFunction> Y_;
    static std::vector<Tabulated2DFunction> PBUB_RS_;
    static std::vector<Tabulated2DFunction> PBUB_RV_;
    static std::vector<Tabulated2DFunction> VISCO_;
    static std::vector<Tabulated2DFunction> VISCG_;
    static std::vector<Tabulated2DFunction> BO_;
    static std::vector<Tabulated2DFunction> BG_;
    static std::vector<Tabulated2DFunction> RS_;
    static std::vector<Tabulated2DFunction> RV_;

    static std::vector<Scalar> zReferenceDensity_;

    static std::vector<Scalar> zLim_;
    static std::vector<TabulatedFunction> oilCmp_;
    static std::vector<TabulatedFunction> gasCmp_;
};

template <class TypeTag, bool enableSsaSolventV>
std::vector<typename BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::Tabulated2DFunction>
BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::X_;

template <class TypeTag, bool enableSsaSolventV>
std::vector<typename BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::Tabulated2DFunction>
BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::Y_;

template <class TypeTag, bool enableSsaSolventV>
std::vector<typename BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::Tabulated2DFunction>
BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::PBUB_RS_;

template <class TypeTag, bool enableSsaSolventV>
std::vector<typename BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::Tabulated2DFunction>
BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::PBUB_RV_;

template <class TypeTag, bool enableSsaSolventV>
std::vector<typename BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::Tabulated2DFunction>
BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::VISCO_;

template <class TypeTag, bool enableSsaSolventV>
std::vector<typename BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::Tabulated2DFunction>
BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::VISCG_;

template <class TypeTag, bool enableSsaSolventV>
std::vector<typename BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::Tabulated2DFunction>
BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::BO_;

template <class TypeTag, bool enableSsaSolventV>
std::vector<typename BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::Tabulated2DFunction>
BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::BG_;

template <class TypeTag, bool enableSsaSolventV>
std::vector<typename BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::Tabulated2DFunction>
BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::RS_;

template <class TypeTag, bool enableSsaSolventV>
std::vector<typename BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::Tabulated2DFunction>
BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::RV_;

template <class TypeTag, bool enableSsaSolventV>
std::vector<typename BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::Scalar>
BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::zReferenceDensity_;

template <class TypeTag, bool enableSsaSolventV>
std::vector<typename BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::Scalar>
BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::zLim_;

template <class TypeTag, bool enableSsaSolventV>
std::vector<typename BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::TabulatedFunction>
BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::oilCmp_;

template <class TypeTag, bool enableSsaSolventV>
std::vector<typename BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::TabulatedFunction>
BlackOilSsaSolventModule<TypeTag, enableSsaSolventV>::gasCmp_;

/*!
 * \ingroup BlackOil
 * \class Ewoms::BlackOilSsaSolventIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        solvents extension of the black-oil model.
 */
template <class TypeTag, bool enableSsaSolventV = GET_PROP_VALUE(TypeTag, EnableSsaSolvent)>
class BlackOilSsaSolventIntensiveQuantities
{
    using Implementation = GET_PROP_TYPE(TypeTag, IntensiveQuantities);

    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);
    using PrimaryVariables = GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
    using MaterialLaw = GET_PROP_TYPE(TypeTag, MaterialLaw);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);
    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);

    using SsaSolventModule = BlackOilSsaSolventModule<TypeTag>;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    static constexpr int zFractionIdx = Indices::zFractionIdx;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr double cutOff = 1e-12;

public:
    /*!
     * \brief Compute extended pvt properties from table lookups.
     *
     *  At this point the pressures of the fluid state are correct.
     */
    void zFractionUpdate_(const ElementContext& elemCtx,
                          unsigned dofIdx,
                          unsigned timeIdx)
    {
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        unsigned pvtRegionIdx = priVars.pvtRegionIndex();
        auto& fs = asImp_().fluidState_;

        zFraction_ = priVars.makeEvaluation(zFractionIdx, timeIdx);

        oilViscosity_ = SsaSolventModule::oilViscosity(pvtRegionIdx, fs.pressure(oilPhaseIdx), zFraction_);
        gasViscosity_ = SsaSolventModule::gasViscosity(pvtRegionIdx, fs.pressure(gasPhaseIdx), zFraction_);

        bo_ = SsaSolventModule::bo(pvtRegionIdx, fs.pressure(oilPhaseIdx), zFraction_);
        bg_ = SsaSolventModule::bg(pvtRegionIdx, fs.pressure(gasPhaseIdx), zFraction_);

        bz_ = SsaSolventModule::bg(pvtRegionIdx, fs.pressure(oilPhaseIdx), Evaluation{0.99});

        if (FluidSystem::enableDissolvedGas())
            rs_ = SsaSolventModule::rs(pvtRegionIdx, fs.pressure(oilPhaseIdx), zFraction_);
        else
            rs_ = 0.0;

        if (FluidSystem::enableVaporizedOil())
            rv_ = SsaSolventModule::rv(pvtRegionIdx, fs.pressure(gasPhaseIdx), zFraction_);
        else
            rv_ = 0.0;

        xVolume_ = SsaSolventModule::xVolume(pvtRegionIdx, fs.pressure(oilPhaseIdx), zFraction_);
        yVolume_ = SsaSolventModule::yVolume(pvtRegionIdx, fs.pressure(oilPhaseIdx), zFraction_);

        Evaluation pbub = fs.pressure(oilPhaseIdx);

        if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg) {
           static const Scalar thresholdWaterFilledCell = 1.0 - 1e-6;
           Scalar Sw = 0.0;
           if (Indices::waterEnabled)
              Sw = priVars.makeEvaluation(Indices::waterSaturationIdx, timeIdx).value();

           if (Sw >= thresholdWaterFilledCell)
              rs_ = 0.0;  // water only, zero rs_ ...
        }

        if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Rs) {
           rs_ = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
           const Evaluation zLim = SsaSolventModule::zLim(pvtRegionIdx);
           if (zFraction_ > zLim) {
             pbub = SsaSolventModule::pbubRs(pvtRegionIdx, zLim, rs_);
           } else {
             pbub = SsaSolventModule::pbubRs(pvtRegionIdx, zFraction_, rs_);
           }
           bo_ = SsaSolventModule::bo(pvtRegionIdx, pbub, zFraction_) + SsaSolventModule::oilCmp(pvtRegionIdx, zFraction_)*(fs.pressure(oilPhaseIdx)-pbub);

           xVolume_ = SsaSolventModule::xVolume(pvtRegionIdx, pbub, zFraction_);
        }

        if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_pg_Rv) {
           rv_ = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
           Evaluation rvsat = SsaSolventModule::rv(pvtRegionIdx, pbub, zFraction_);
           bg_ = SsaSolventModule::bg(pvtRegionIdx, pbub, zFraction_) + SsaSolventModule::gasCmp(pvtRegionIdx, zFraction_)*(rv_-rvsat);

           yVolume_ = SsaSolventModule::yVolume(pvtRegionIdx, pbub, zFraction_);
        }
    }

    /*!
     * \brief Re-compute face densities to account for zFraction dependency.
     *
     * At this point the pressures and saturations of the fluid state are correct.
     */
    void zPvtUpdate_()
    {
        const auto& iq = asImp_();
        auto& fs = asImp_().fluidState_;

        unsigned pvtRegionIdx = iq.pvtRegionIndex();
        zRefDensity_ = SsaSolventModule::referenceDensity(pvtRegionIdx);

        fs.setInvB(oilPhaseIdx, 1.0/bo_);
        fs.setInvB(gasPhaseIdx, 1.0/bg_);

        fs.setDensity(oilPhaseIdx,
                      fs.invB(oilPhaseIdx)
                      *(FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx)
                        + (1.0-xVolume_)*fs.Rs()*FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx)
                        + xVolume_*fs.Rs()*zRefDensity_ ));
        fs.setDensity(gasPhaseIdx,
                      fs.invB(gasPhaseIdx)
                      *(FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx)*(1.0-yVolume_)+yVolume_*zRefDensity_
                        + FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx)*fs.Rv()));
    }

    const Evaluation& zFraction() const
    { return zFraction_; }

    const Evaluation& xVolume() const
    { return xVolume_; }

    const Evaluation& yVolume() const
    { return yVolume_; }

    const Evaluation& oilViscosity() const
    { return oilViscosity_; }

    const Evaluation& gasViscosity() const
    { return gasViscosity_; }

    const Evaluation& bo() const
    { return bo_; }

    const Evaluation& bg() const
    { return bg_; }

    const Evaluation& rs() const
    { return rs_; }

    const Evaluation& rv() const
    { return rv_; }

    const Evaluation zPureInvFormationVolumeFactor() const
    { return 1.0/bz_; }

    const Scalar& zRefDensity() const
    { return zRefDensity_; }

private:

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    // Abstract "mass fraction" accounting for the solvent component. The relation between this
    // quantity and the actual mass fraction of solvent, is implicitly defined from the specific
    // pvt measurements as provided by kw PVTSOL.
    Evaluation zFraction_;

    // The solvent component is assumed gas at surface conditions
    Evaluation xVolume_; // Solvent volume fraction of Rs
    Evaluation yVolume_; // Solvent volume fraction of Sg/Bg

    // Standard black oil parameters modified for presence of solvent
    Evaluation oilViscosity_;
    Evaluation gasViscosity_;
    Evaluation bo_;
    Evaluation bg_;
    Evaluation rs_;
    Evaluation rv_;

    // Properties of pure solvent
    Evaluation bz_;
    Scalar zRefDensity_;
};

template <class TypeTag>
class BlackOilSsaSolventIntensiveQuantities<TypeTag, false>
{
    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);
    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);

public:

    void zPvtUpdate_()
    { }

    void zFractionUpdate_(const ElementContext& elemCtx EWOMS_UNUSED,
                          unsigned scvIdx EWOMS_UNUSED,
                          unsigned timeIdx EWOMS_UNUSED)
    { }

    const Evaluation& xVolume() const
    { throw std::runtime_error("xVolume() called but SSA solvents are disabled"); }

    const Evaluation& yVolume() const
    { throw std::runtime_error("yVolume() called but SSA solvents are disabled"); }

    const Evaluation& oilViscosity() const
    { throw std::runtime_error("oilViscosity() called but SSA solvents are disabled"); }

    const Evaluation& gasViscosity() const
    { throw std::runtime_error("gasViscosity() called but SSA solvents are disabled"); }

    const Evaluation& rs() const
    { throw std::runtime_error("rs() called but SSA solvents are disabled"); }

    const Evaluation& rv() const
    { throw std::runtime_error("rv() called but SSA solvents are disabled"); }

    const Evaluation& zPureInvFormationVolumeFactor() const
    { throw std::runtime_error("zPureInvFormationVolumeFactor() called but SSA solvents are disabled"); }

    const Evaluation& zFraction() const
    { throw std::runtime_error("zFraction() called but SSA solvents are disabled"); }

    const Evaluation& zInverseFormationVolumeFactor() const
     { throw std::runtime_error("zInverseFormationVolumeFactor() called but SSA solvents are disabled"); }

    const Scalar& zRefDensity() const
     { throw std::runtime_error("zRefDensity() called but SSA solvents are disabled"); }
};

/*!
 * \ingroup BlackOil
 * \class Ewoms::BlackOilSsaSolventExtensiveQuantities
 *
 * \brief Provides the solvent specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag, bool enableSsaSolventV = GET_PROP_VALUE(TypeTag, EnableSsaSolvent)>
class BlackOilSsaSolventExtensiveQuantities
{
    using Implementation = GET_PROP_TYPE(TypeTag, ExtensiveQuantities);

    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);
    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);
    using IntensiveQuantities = GET_PROP_TYPE(TypeTag, IntensiveQuantities);
    using ExtensiveQuantities = GET_PROP_TYPE(TypeTag, ExtensiveQuantities);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
    using GridView = GET_PROP_TYPE(TypeTag, GridView);

    using Toolbox = Ewoms::MathToolbox<Evaluation>;

    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int dimWorld = GridView::dimensionworld;

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldVector<Evaluation, dimWorld> DimEvalVector;

public:

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

};

template <class TypeTag>
class BlackOilSsaSolventExtensiveQuantities<TypeTag, false>
{
    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);
    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);

public:

};

} // namespace Ewoms

#endif
