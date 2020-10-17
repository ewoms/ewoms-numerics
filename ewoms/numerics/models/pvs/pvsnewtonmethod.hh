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
 * \copydoc Ewoms::PvsNewtonMethod
 */
#ifndef EWOMS_PVS_NEWTON_METHOD_HH
#define EWOMS_PVS_NEWTON_METHOD_HH

#include "pvsproperties.hh"

namespace Ewoms {

/*!
 * \ingroup PvsModel
 *
 * \brief A newton solver which is specific to the compositional
 *        multi-phase PVS model.
 */
template <class TypeTag>
class PvsNewtonMethod : public GET_PROP_TYPE(TypeTag, DiscNewtonMethod)
{
    using ParentType = GET_PROP_TYPE(TypeTag, DiscNewtonMethod);
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using SolutionVector = GET_PROP_TYPE(TypeTag, SolutionVector);
    using PrimaryVariables = GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using EqVector = GET_PROP_TYPE(TypeTag, EqVector);
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);

    enum { numPhases = FluidSystem::numPhases };

    // primary variable indices
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { switch0Idx = Indices::switch0Idx };

public:
    PvsNewtonMethod(Simulator& simulator) : ParentType(simulator)
    {}

protected:
    friend NewtonMethod<TypeTag>;
    friend ParentType;

    /*!
     * \copydoc FvBaseNewtonMethod::updatePrimaryVariables_
     */
    void updatePrimaryVariables_(unsigned globalDofIdx EWOMS_UNUSED,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& currentResidual EWOMS_UNUSED)
    {
        // normal Newton-Raphson update
        nextValue = currentValue;
        nextValue -= update;

        ////
        // put crash barriers along the update path
        ////
        // saturations: limit the change of any saturation to at most 20%
        Scalar sumSatDelta = 0.0;
        Scalar maxSatDelta = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
            if (!currentValue.phaseIsPresent(phaseIdx))
                continue;

            maxSatDelta = std::max(std::abs(update[switch0Idx + phaseIdx]),
                                   maxSatDelta);
            sumSatDelta += update[switch0Idx + phaseIdx];
        }
        maxSatDelta = std::max(std::abs(- sumSatDelta), maxSatDelta);

        if (maxSatDelta > 0.2) {
            Scalar alpha = 0.2/maxSatDelta;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
                if (!currentValue.phaseIsPresent(phaseIdx))
                    continue;

                nextValue[switch0Idx + phaseIdx] =
                    currentValue[switch0Idx + phaseIdx]
                    - alpha*update[switch0Idx + phaseIdx];
            }
        }

        // limit pressure reference change to 20% of the total value per iteration
        clampValue_(nextValue[pressure0Idx],
                    currentValue[pressure0Idx]*0.8,
                    currentValue[pressure0Idx]*1.2);
    }

    /*!
     * \copydoc NewtonMethod::endIteration_
     */
    void endIteration_(SolutionVector& uCurrentIter,
                       const SolutionVector& uLastIter)
    {
        ParentType::endIteration_(uCurrentIter, uLastIter);
        this->problem().model().switchPrimaryVars_();
    }

    void clampValue_(Scalar& val, Scalar minVal, Scalar maxVal) const
    { val = std::max(minVal, std::min(val, maxVal)); }
};
} // namespace Ewoms

#endif
