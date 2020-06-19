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
 * \copydoc Ewoms::FvBaseNewtonMethod
 */
#ifndef EWOMS_FV_BASE_NEWTON_METHOD_HH
#define EWOMS_FV_BASE_NEWTON_METHOD_HH

#include "fvbasenewtonconvergencewriter.hh"

#include <ewoms/numerics/nonlinear/newtonmethod.hh>
#include <ewoms/common/propertysystem.hh>

namespace Ewoms {

template <class TypeTag>
class FvBaseNewtonMethod;

template <class TypeTag>
class FvBaseNewtonConvergenceWriter;
} // namespace Ewoms

BEGIN_PROPERTIES

//! create a type tag for the Newton method of the finite-volume discretization
NEW_TYPE_TAG(FvBaseNewtonMethod, INHERITS_FROM(NewtonMethod));

//! The class dealing with the balance equations
NEW_PROP_TAG(Model);

//! The class storing primary variables plus pseudo primary variables
NEW_PROP_TAG(PrimaryVariables);

//! The class storing values of conservation equations (e.g., a "naked" primary varible
//! vector)
NEW_PROP_TAG(EqVector);

//! The number of balance equations.
NEW_PROP_TAG(NumEq);

//! The discretization specific part of he implementing the Newton algorithm
NEW_PROP_TAG(DiscNewtonMethod);

//! The class implementing the Newton algorithm
NEW_PROP_TAG(NewtonMethod);

// set default values
SET_TYPE_PROP(FvBaseNewtonMethod, DiscNewtonMethod,
              Ewoms::FvBaseNewtonMethod<TypeTag>);
SET_TYPE_PROP(FvBaseNewtonMethod, NewtonMethod,
              GET_PROP_TYPE(TypeTag, DiscNewtonMethod));
SET_TYPE_PROP(FvBaseNewtonMethod, NewtonConvergenceWriter,
              Ewoms::FvBaseNewtonConvergenceWriter<TypeTag>);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief A Newton method for models using a finite volume discretization.
 *
 * This class is sufficient for most models which use an Element or a
 * Vertex Centered Finite Volume discretization.
 */
template <class TypeTag>
class FvBaseNewtonMethod : public NewtonMethod<TypeTag>
{
    using ParentType = Ewoms::NewtonMethod<TypeTag>;
    using Implementation = GET_PROP_TYPE(TypeTag, NewtonMethod);

    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using Model = GET_PROP_TYPE(TypeTag, Model);
    using NewtonMethod = GET_PROP_TYPE(TypeTag, NewtonMethod);
    using GlobalEqVector = GET_PROP_TYPE(TypeTag, GlobalEqVector);
    using SolutionVector = GET_PROP_TYPE(TypeTag, SolutionVector);
    using PrimaryVariables = GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using EqVector = GET_PROP_TYPE(TypeTag, EqVector);

public:
    FvBaseNewtonMethod(Simulator& simulator)
        : ParentType(simulator)
    { }

protected:
    friend class Ewoms::NewtonMethod<TypeTag>;

    /*!
     * \brief Update the current solution with a delta vector.
     *
     * The error estimates required for the converged() and
     * proceed() methods should be updated inside this method.
     *
     * Different update strategies, such as line search and chopped
     * updates can be implemented. The default behavior is just to
     * subtract deltaU from uLastIter, i.e.
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param nextSolution The solution vector at the end of the current iteration
     * \param currentSolution The solution vector at the beginning of the current iteration
     * \param solutionUpdate The delta as calculated by solving the linear system of
     *                       equations. This parameter also stores the updated solution.
     * \param currentResidual The residual (i.e., right-hand-side) of the current solution.
     */
    void update_(SolutionVector& nextSolution,
                 const SolutionVector& currentSolution,
                 const GlobalEqVector& solutionUpdate,
                 const GlobalEqVector& currentResidual)
    {
        ParentType::update_(nextSolution, currentSolution, solutionUpdate, currentResidual);

        // make sure that the intensive quantities get recalculated at the next
        // linearization
        if (model_().storeIntensiveQuantities()) {
            for (unsigned dofIdx = 0; dofIdx < model_().numGridDof(); ++dofIdx)
                model_().setIntensiveQuantitiesCacheEntryValidity(dofIdx,
                                                                  /*timeIdx=*/0,
                                                                  /*valid=*/false);
        }
    }

    /*!
     * \brief Indicates the beginning of a Newton iteration.
     */
    void beginIteration_()
    {
        model_().syncOverlap();

        ParentType::beginIteration_();
    }

    /*!
     * \brief Returns a reference to the model.
     */
    Model& model_()
    { return ParentType::model(); }

    /*!
     * \brief Returns a reference to the model.
     */
    const Model& model_() const
    { return ParentType::model(); }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }
};
} // namespace Ewoms

#endif
