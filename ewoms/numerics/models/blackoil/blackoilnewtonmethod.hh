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
 * \copydoc Ewoms::BlackOilNewtonMethod
 */
#ifndef EWOMS_BLACK_OIL_NEWTON_METHOD_HH
#define EWOMS_BLACK_OIL_NEWTON_METHOD_HH

#include "blackoilproperties.hh"

#include <ewoms/common/signum.hh>

#include <ewoms/common/unused.hh>

BEGIN_PROPERTIES

NEW_PROP_TAG(DpMaxRel);
NEW_PROP_TAG(DsMax);
NEW_PROP_TAG(PriVarOscilationThreshold);

SET_SCALAR_PROP(NewtonMethod, DpMaxRel, 0.3);
SET_SCALAR_PROP(NewtonMethod, DsMax, 0.2);
SET_SCALAR_PROP(NewtonMethod, PriVarOscilationThreshold, 1e-5);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup BlackOilModel
 *
 * \brief A newton solver which is specific to the black oil model.
 */
template <class TypeTag>
class BlackOilNewtonMethod : public GET_PROP_TYPE(TypeTag, DiscNewtonMethod)
{
    using ParentType = GET_PROP_TYPE(TypeTag, DiscNewtonMethod);
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using SolutionVector = GET_PROP_TYPE(TypeTag, SolutionVector);
    using GlobalEqVector = GET_PROP_TYPE(TypeTag, GlobalEqVector);
    using PrimaryVariables = GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using EqVector = GET_PROP_TYPE(TypeTag, EqVector);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);

    static const unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);

public:
    BlackOilNewtonMethod(Simulator& simulator) : ParentType(simulator)
    {
        priVarOscilationThreshold_ = EWOMS_GET_PARAM(TypeTag, Scalar, PriVarOscilationThreshold);
        dpMaxRel_ = EWOMS_GET_PARAM(TypeTag, Scalar, DpMaxRel);
        dsMax_ = EWOMS_GET_PARAM(TypeTag, Scalar, DsMax);
    }

    /*!
     * \copydoc NewtonMethod::finishInit()
     */
    void finishInit()
    {
        ParentType::finishInit();

        wasSwitched_.resize(this->model().numTotalDof());
        std::fill(wasSwitched_.begin(), wasSwitched_.end(), false);
    }

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, DpMaxRel, "Maximum relative change of pressure in a single iteration");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, DsMax, "Maximum absolute change of any saturation in a single iteration");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, PriVarOscilationThreshold,
                             "The threshold value for the primary variable switching conditions after its meaning has switched to hinder oscilations");
    }

    /*!
     * \brief Returns the number of degrees of freedom for which the
     *        interpretation has changed for the most recent iteration.
     */
    unsigned numPriVarsSwitched() const
    { return numPriVarsSwitched_; }

protected:
    friend NewtonMethod<TypeTag>;
    friend ParentType;

    /*!
     * \copydoc FvBaseNewtonMethod::beginIteration_
     */
    void beginIteration_()
    {
        numPriVarsSwitched_ = 0;
        ParentType::beginIteration_();
    }

    /*!
     * \copydoc FvBaseNewtonMethod::endIteration_
     */
    void endIteration_(SolutionVector& uCurrentIter,
                       const SolutionVector& uLastIter)
    {
#if HAVE_MPI
        // in the MPI enabled case we need to add up the number of DOF
        // for which the interpretation changed over all processes.
        int localSwitched = numPriVarsSwitched_;
        MPI_Allreduce(&localSwitched,
                      &numPriVarsSwitched_,
                      /*num=*/1,
                      MPI_INT,
                      MPI_SUM,
                      MPI_COMM_WORLD);
#endif // HAVE_MPI

        this->simulator_.model().newtonMethod().endIterMsg()
            << ", num switched=" << numPriVarsSwitched_;

        ParentType::endIteration_(uCurrentIter, uLastIter);
    }

public:
    void update_(SolutionVector& nextSolution,
                 const SolutionVector& currentSolution,
                 const GlobalEqVector& solutionUpdate,
                 const GlobalEqVector& currentResidual)
    {
        const auto& comm = this->simulator_.gridView().comm();

        int succeeded;
        try {
            ParentType::update_(nextSolution,
                                currentSolution,
                                solutionUpdate,
                                currentResidual);
            succeeded = 1;
        }
        catch (...) {
            std::cout << "Newton update threw an exception on rank "
                      << comm.rank() << "\n";
            succeeded = 0;
        }
        succeeded = comm.min(succeeded);

        if (!succeeded)
            throw Ewoms::NumericalIssue("A process did not succeed in adapting the primary variables");

        numPriVarsSwitched_ = comm.sum(numPriVarsSwitched_);
    }

protected:
    /*!
     * \copydoc FvBaseNewtonMethod::updatePrimaryVariables_
     */
    void updatePrimaryVariables_(unsigned globalDofIdx,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& currentResidual)
    {
        static constexpr bool enableSolvent = Indices::solventSaturationIdx >= 0;
        static constexpr bool enableSsaSolvent = Indices::zFractionIdx >= 0;
        static constexpr bool enablePolymer = Indices::polymerConcentrationIdx >= 0;
        static constexpr bool enablePolymerWeight = Indices::polymerMoleWeightIdx >= 0;
        static constexpr bool enableEnergy = Indices::temperatureIdx >= 0;
        static constexpr bool enableFoam = Indices::foamConcentrationIdx >= 0;
        static constexpr bool enableBrine = Indices::saltConcentrationIdx >= 0;

        currentValue.checkDefined();
        Ewoms::Valgrind::CheckDefined(update);
        Ewoms::Valgrind::CheckDefined(currentResidual);

        // saturation delta for each phase
        Scalar deltaSw = 0.0;
        Scalar deltaSo = 0.0;
        Scalar deltaSg = 0.0;
        Scalar deltaSs = 0.0;

        if (Indices::waterEnabled) {
            deltaSw = update[Indices::waterSaturationIdx];
            deltaSo = -deltaSw;
        }

        if (Indices::gasEnabled && currentValue.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg) {
            deltaSg = update[Indices::compositionSwitchIdx];
            deltaSo -= deltaSg;
        }

        if (enableSolvent) {
            deltaSs = update[Indices::solventSaturationIdx];
            deltaSo -= deltaSs;
        }

        // maximum saturation delta
        Scalar maxSatDelta = std::max(std::abs(deltaSg), std::abs(deltaSo));
        maxSatDelta = std::max(maxSatDelta, std::abs(deltaSw));
        maxSatDelta = std::max(maxSatDelta, std::abs(deltaSs));

        // scaling factor for saturation deltas to make sure that none of them exceeds
        // the specified threshold value.
        Scalar satAlpha = 1.0;
        if (maxSatDelta > dsMax_)
            satAlpha = dsMax_/maxSatDelta;

        for (int pvIdx = 0; pvIdx < int(numEq); ++pvIdx) {
            // calculate the update of the current primary variable. For the black-oil
            // model we limit the pressure delta relative to the pressure's current
            // absolute value (Default: 30%) and saturation deltas to an absolute change
            // (Default: 20%). Further, we ensure that the R factors, solvent
            // "saturation" and polymer concentration do not become negative after the
            // update.
            Scalar delta = update[pvIdx];

            // limit pressure delta
            if (pvIdx == Indices::pressureSwitchIdx) {
                if (std::abs(delta) > dpMaxRel_*currentValue[pvIdx])
                    delta = Ewoms::signum(delta)*dpMaxRel_*currentValue[pvIdx];
            }
            // water saturation delta
            else if (pvIdx == Indices::waterSaturationIdx)
                delta *= satAlpha;
            else if (pvIdx == Indices::compositionSwitchIdx) {
                // the switching primary variable for composition is tricky because the
                // "reasonable" value ranges it exhibits vary widely depending on its
                // interpretation since it can represent Sg, Rs or Rv. For now, we only
                // limit saturation deltas and ensure that the R factors do not become
                // negative.
                if (currentValue.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg)
                    delta *= satAlpha;
                else {
                    if (delta > currentValue[Indices::compositionSwitchIdx])
                        delta = currentValue[Indices::compositionSwitchIdx];
                }
            }
            else if (enableSolvent && pvIdx == Indices::solventSaturationIdx)
                // solvent saturation updates are also subject to the Appleyard chop
                delta *= satAlpha;
            else if (enableSsaSolvent && pvIdx == Indices::zFractionIdx) {
                // z fraction updates are also subject to the Appleyard chop
                if (delta > currentValue[Indices::zFractionIdx])
                        delta = currentValue[Indices::zFractionIdx];
                if (delta < currentValue[Indices::zFractionIdx]-1.0)
                        delta = currentValue[Indices::zFractionIdx]-1.0;
            }
            else if (enablePolymerWeight && pvIdx == Indices::polymerMoleWeightIdx) {
                const double sign = delta >= 0. ? 1. : -1.;
                // maximum change of polymer molecular weight, the unit is MDa.
                // applying this limit to stabilize the simulation. The value itself is still experimental.
                const double maxMolarWeightChange = 100.0;
                delta = sign * std::min(std::abs(delta), maxMolarWeightChange);
                delta *= satAlpha;
            }

            // do the actual update
            nextValue[pvIdx] = currentValue[pvIdx] - delta;

            // keep the solvent saturation between 0 and 1
            if (enableSolvent && pvIdx == Indices::solventSaturationIdx)
                nextValue[pvIdx] = std::min(std::max(nextValue[pvIdx], 0.0), 1.0);

            // keep the z fraction between 0 and 1
            if (enableSsaSolvent && pvIdx == Indices::zFractionIdx)
                nextValue[pvIdx] = std::min(std::max(nextValue[pvIdx], 0.0), 1.0);

            // keep the polymer concentration above 0
            if (enablePolymer && pvIdx == Indices::polymerConcentrationIdx)
                nextValue[pvIdx] = std::max(nextValue[pvIdx], 0.0);

            if (enablePolymerWeight && pvIdx == Indices::polymerMoleWeightIdx) {
                nextValue[pvIdx] = std::max(nextValue[pvIdx], 0.0);
                const double polymerConcentration = nextValue[Indices::polymerConcentrationIdx];
                if (polymerConcentration < 1.e-10)
                    nextValue[pvIdx] = 0.0;
            }

            // keep the foam concentration above 0
            if (enableFoam && pvIdx == Indices::foamConcentrationIdx)
                nextValue[pvIdx] = std::max(nextValue[pvIdx], 0.0);

            // keep the salt concentration above 0
            if (enableBrine && pvIdx == Indices::saltConcentrationIdx)
                nextValue[pvIdx] = std::max(nextValue[pvIdx], 0.0);

            // keep the temperature above 100 and below 1000 Kelvin
            if (enableEnergy && pvIdx == Indices::temperatureIdx)
                nextValue[pvIdx] = std::max(std::min(nextValue[pvIdx], 1000.0), 100.0);
        }

        // switch the new primary variables to something which is physically meaningful.
        // use a threshold value after a switch to make it harder to switch back
        // immediately.
        if (wasSwitched_[globalDofIdx])
            wasSwitched_[globalDofIdx] = nextValue.adaptPrimaryVariables(this->problem(), globalDofIdx, priVarOscilationThreshold_);
        else
            wasSwitched_[globalDofIdx] = nextValue.adaptPrimaryVariables(this->problem(), globalDofIdx);

        if (wasSwitched_[globalDofIdx])
            ++ numPriVarsSwitched_;

        nextValue.checkDefined();
    }

private:
    int numPriVarsSwitched_;

    Scalar priVarOscilationThreshold_;
    Scalar dpMaxRel_;
    Scalar dsMax_;

    // keep track of cells where the primary variable meaning has changed
    // to detect and hinder oscillations
    std::vector<bool> wasSwitched_;
};
} // namespace Ewoms

#endif
