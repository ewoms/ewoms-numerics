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
 * \copydoc Ewoms::NcpModel
 */
#ifndef EWOMS_NCP_MODEL_HH
#define EWOMS_NCP_MODEL_HH

#include <ewoms/common/densead/math.hh>

#include "ncpproperties.hh"
#include "ncplocalresidual.hh"
#include "ncpextensivequantities.hh"
#include "ncpprimaryvariables.hh"
#include "ncpboundaryratevector.hh"
#include "ncpratevector.hh"
#include "ncpintensivequantities.hh"
#include "ncpnewtonmethod.hh"
#include "ncpindices.hh"

#include <ewoms/numerics/common/multiphasebasemodel.hh>
#include <ewoms/numerics/common/energymodule.hh>
#include <ewoms/numerics/common/diffusionmodule.hh>
#include <ewoms/numerics/io/vtkcompositionmodule.hh>
#include <ewoms/numerics/io/vtkenergymodule.hh>
#include <ewoms/numerics/io/vtkdiffusionmodule.hh>

#include <ewoms/common/valgrind.hh>
#include <ewoms/common/exceptions.hh>

#include <dune/common/fvector.hh>

#include <sstream>
#include <string>
#include <vector>
#include <array>

namespace Ewoms {
template <class TypeTag>
class NcpModel;
}

BEGIN_PROPERTIES

/*!
 * \brief Define the type tag for the compositional NCP model.
 */
NEW_TYPE_TAG(NcpModel, INHERITS_FROM(MultiPhaseBaseModel,
                                     VtkComposition,
                                     VtkEnergy,
                                     VtkDiffusion));

//! Use the Ncp local jacobian operator for the compositional NCP model
SET_TYPE_PROP(NcpModel,
              LocalResidual,
              Ewoms::NcpLocalResidual<TypeTag>);

//! Use the Ncp specific newton method for the compositional NCP model
SET_TYPE_PROP(NcpModel, NewtonMethod, Ewoms::NcpNewtonMethod<TypeTag>);

//! the Model property
SET_TYPE_PROP(NcpModel, Model, Ewoms::NcpModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(NcpModel, BaseProblem, Ewoms::MultiPhaseBaseProblem<TypeTag>);

//! Disable the energy equation by default
SET_BOOL_PROP(NcpModel, EnableEnergy, false);

//! disable diffusion by default
SET_BOOL_PROP(NcpModel, EnableDiffusion, false);

//! the RateVector property
SET_TYPE_PROP(NcpModel, RateVector, Ewoms::NcpRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(NcpModel, BoundaryRateVector, Ewoms::NcpBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(NcpModel, PrimaryVariables, Ewoms::NcpPrimaryVariables<TypeTag>);

//! the IntensiveQuantities property
SET_TYPE_PROP(NcpModel, IntensiveQuantities, Ewoms::NcpIntensiveQuantities<TypeTag>);

//! the ExtensiveQuantities property
SET_TYPE_PROP(NcpModel, ExtensiveQuantities, Ewoms::NcpExtensiveQuantities<TypeTag>);

//! The indices required by the compositional NCP model
SET_TYPE_PROP(NcpModel, Indices, Ewoms::NcpIndices<TypeTag, 0>);

//! The unmodified weight for the pressure primary variable
SET_SCALAR_PROP(NcpModel, NcpPressureBaseWeight, 1.0);
//! The weight for the saturation primary variables
SET_SCALAR_PROP(NcpModel, NcpSaturationsBaseWeight, 1.0);
//! The unmodified weight for the fugacity primary variables
SET_SCALAR_PROP(NcpModel, NcpFugacitiesBaseWeight, 1.0e-6);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup NcpModel
 *
 * \brief A compositional multi-phase model based on non-linear
 *        complementarity functions.
 *
 * This model implements a \f$M\f$-phase flow of a fluid mixture
 * composed of \f$N\f$ chemical species. The phases are denoted by
 * lower index \f$\alpha \in \{ 1, \dots, M \}\f$. All fluid phases
 * are mixtures of \f$N \geq M - 1\f$ chemical species which are
 * denoted by the upper index \f$\kappa \in \{ 1, \dots, N \} \f$.
 *
 *
 * By default, the standard multi-phase Darcy approach is used to determine
 * the velocity, i.e.
 * \f[
 *   \mathbf{v}_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 *   \left(\mathbf{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right) \;,
 * \f]
 * although the actual approach which is used can be specified via the
 * \c FluxModule property. For example, the velocity model can by
 * changed to the Forchheimer approach by
 * \code
 * SET_TYPE_PROP(MyProblemTypeTag, FluxModule, Ewoms::ForchheimerFluxModule<TypeTag>);
 * \endcode
 *
 * The core of the model is the conservation mass of each component by
 * means of the equation
 * \f[
 * \sum_\alpha \frac{\partial\;\phi c_\alpha^\kappa S_\alpha }{\partial t}
 * - \sum_\alpha \mathrm{div} \left\{ c_\alpha^\kappa \mathbf{v}_\alpha \right\}
 * - q^\kappa = 0 \;.
 * \f]
 *
 * For the missing \f$M\f$ model assumptions, the model uses
 * non-linear complementarity functions. These are based on the
 * observation that if a fluid phase is not present, the sum of the
 * mole fractions of this fluid phase is smaller than \f$1\f$, i.e.
 * \f[ \forall \alpha: S_\alpha = 0 \implies \sum_\kappa
 * x_\alpha^\kappa \leq 1 \f]
 *
 * Also, if a fluid phase may be present at a given spatial location
 * its saturation must be non-negative:
 * \f[ \forall \alpha: \sum_\kappa x_\alpha^\kappa = 1 \implies S_\alpha \geq 0
 *\f]
 *
 * Since at any given spatial location, a phase is always either
 * present or not present, one of the strict equalities on the
 * right hand side is always true, i.e.
 * \f[
 * \forall \alpha: S_\alpha \left( \sum_\kappa x_\alpha^\kappa - 1 \right) = 0
 * \f]
 * always holds.
 *
 * These three equations constitute a non-linear complementarity
 * problem, which can be solved using so-called non-linear
 * complementarity functions \f$\Phi(a, b)\f$. Such functions have the property
 * \f[\Phi(a,b) = 0 \iff a \geq0 \land b \geq0  \land a \cdot b = 0 \f]
 *
 * Several non-linear complementarity functions have been suggested,
 * e.g. the Fischer-Burmeister function
 * \f[ \Phi(a,b) = a + b - \sqrt{a^2 + b^2} \;. \f]
 * This model uses
 * \f[ \Phi(a,b) = \min \{a,  b \}\;, \f]
 * because of its piecewise linearity.
 *
 * The model assumes local thermodynamic equilibrium and uses the
 * following primary variables:
 * - The pressure of the first phase \f$p_1\f$
 * - The component fugacities \f$f^1, \dots, f^{N}\f$
 * - The saturations of the first \f$M-1\f$ phases \f$S_1, \dots, S_{M-1}\f$
 * - Temperature \f$T\f$ if the energy equation is enabled
 */
template <class TypeTag>
class NcpModel
    : public MultiPhaseBaseModel<TypeTag>
{
    using ParentType = MultiPhaseBaseModel<TypeTag>;

    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { fugacity0Idx = Indices::fugacity0Idx };
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { ncp0EqIdx = Indices::ncp0EqIdx };
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;

    using Toolbox = Ewoms::MathToolbox<Evaluation>;

    using EnergyModule = Ewoms::EnergyModule<TypeTag, enableEnergy>;
    using DiffusionModule = Ewoms::DiffusionModule<TypeTag, enableDiffusion>;

public:
    NcpModel(Simulator& simulator)
        : ParentType(simulator)
    {}

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        DiffusionModule::registerParameters();
        EnergyModule::registerParameters();

        // register runtime parameters of the VTK output modules
        Ewoms::VtkCompositionModule<TypeTag>::registerParameters();

        if (enableDiffusion)
            Ewoms::VtkDiffusionModule<TypeTag>::registerParameters();

        if (enableEnergy)
            Ewoms::VtkEnergyModule<TypeTag>::registerParameters();
    }

    /*!
     * \copydoc FvBaseDiscretization::finishInit()
     */
    void finishInit()
    {
        ParentType::finishInit();

        minActivityCoeff_.resize(this->numGridDof());
        std::fill(minActivityCoeff_.begin(), minActivityCoeff_.end(), 1.0);
    }

    void adaptGrid()
    {
        ParentType::adaptGrid();
        minActivityCoeff_.resize(this->numGridDof());
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    static std::string name()
    { return "ncp"; }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarName
     */
    std::string primaryVarName(unsigned pvIdx) const
    {
        std::string s;
        if (!(s = EnergyModule::primaryVarName(pvIdx)).empty())
            return s;

        std::ostringstream oss;
        if (pvIdx == pressure0Idx)
            oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/0);
        else if (saturation0Idx <= pvIdx && pvIdx < saturation0Idx + (numPhases - 1))
            oss << "saturation_" << FluidSystem::phaseName(/*phaseIdx=*/pvIdx - saturation0Idx);
        else if (fugacity0Idx <= pvIdx && pvIdx < fugacity0Idx + numComponents)
            oss << "fugacity^" << FluidSystem::componentName(pvIdx - fugacity0Idx);
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::eqName
     */
    std::string eqName(unsigned eqIdx) const
    {
        std::string s;
        if (!(s = EnergyModule::eqName(eqIdx)).empty())
            return s;

        std::ostringstream oss;
        if (conti0EqIdx <= eqIdx && eqIdx < conti0EqIdx + numComponents)
            oss << "continuity^" << FluidSystem::componentName(eqIdx - conti0EqIdx);
        else if (ncp0EqIdx <= eqIdx && eqIdx < ncp0EqIdx + numPhases)
            oss << "ncp_" << FluidSystem::phaseName(/*phaseIdx=*/eqIdx - ncp0EqIdx);
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::updateBegin
     */
    void updateBegin()
    {
        ParentType::updateBegin();

        // find the a reference pressure. The first degree of freedom
        // might correspond to non-interior entities which would lead
        // to an undefined value, so we have to iterate...
        for (unsigned dofIdx = 0; dofIdx < this->numGridDof(); ++ dofIdx) {
            if (this->isLocalDof(dofIdx)) {
                referencePressure_ =
                    this->solution(/*timeIdx=*/0)[dofIdx][/*pvIdx=*/Indices::pressure0Idx];
                break;
            }
        }
    }

    /*!
     * \copydoc FvBaseDiscretization::updatePVWeights
     */
    void updatePVWeights(const ElementContext& elemCtx) const
    {
        for (unsigned dofIdx = 0; dofIdx < elemCtx.numDof(/*timeIdx=*/0); ++dofIdx) {
            unsigned globalIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                minActivityCoeff_[globalIdx][compIdx] = 1e100;
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    const auto& fs = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0).fluidState();

                    minActivityCoeff_[globalIdx][compIdx] =
                        std::min(minActivityCoeff_[globalIdx][compIdx],
                                 Toolbox::value(fs.fugacityCoefficient(phaseIdx, compIdx))
                                 * Toolbox::value(fs.pressure(phaseIdx)));
                    Ewoms::Valgrind::CheckDefined(minActivityCoeff_[globalIdx][compIdx]);
                }
                if (minActivityCoeff_[globalIdx][compIdx] <= 0)
                    throw Ewoms::NumericalIssue("The minumum activity coefficient for component "+std::to_string(compIdx)
                                                +" on DOF "+std::to_string(globalIdx)+" is negative or zero!");
            }
        }
    }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarWeight
     */
    Scalar primaryVarWeight(unsigned globalDofIdx, unsigned pvIdx) const
    {
        Scalar tmp = EnergyModule::primaryVarWeight(*this, globalDofIdx, pvIdx);
        Scalar result;
        if (tmp > 0)
            // energy related quantity
            result = tmp;
        else if (fugacity0Idx <= pvIdx && pvIdx < fugacity0Idx + numComponents) {
            // component fugacity
            unsigned compIdx = pvIdx - fugacity0Idx;
            assert(0 <= compIdx && compIdx <= numComponents);

            Ewoms::Valgrind::CheckDefined(minActivityCoeff_[globalDofIdx][compIdx]);
            static const Scalar fugacityBaseWeight =
                GET_PROP_VALUE(TypeTag, NcpFugacitiesBaseWeight);
            result = fugacityBaseWeight / minActivityCoeff_[globalDofIdx][compIdx];
        }
        else if (Indices::pressure0Idx == pvIdx) {
            static const Scalar pressureBaseWeight = GET_PROP_VALUE(TypeTag, NcpPressureBaseWeight);
            result = pressureBaseWeight / referencePressure_;
        }
        else {
#ifndef NDEBUG
            unsigned phaseIdx = pvIdx - saturation0Idx;
            assert(0 <= phaseIdx && phaseIdx < numPhases - 1);
#endif

            // saturation
            static const Scalar saturationsBaseWeight =
                GET_PROP_VALUE(TypeTag, NcpSaturationsBaseWeight);
            result = saturationsBaseWeight;
        }

        assert(std::isfinite(result));
        assert(result > 0);

        return result;
    }

    /*!
     * \copydoc FvBaseDiscretization::eqWeight
     */
    Scalar eqWeight(unsigned globalDofIdx, unsigned eqIdx) const
    {
        Scalar tmp = EnergyModule::eqWeight(*this, globalDofIdx, eqIdx);
        if (tmp > 0)
            // an energy related equation
            return tmp;
        // an NCP
        else if (ncp0EqIdx <= eqIdx && eqIdx < Indices::ncp0EqIdx + numPhases)
            return 1.0;

        // a mass conservation equation
        unsigned compIdx = eqIdx - Indices::conti0EqIdx;
        assert(0 <= compIdx && compIdx <= numComponents);

        // make all kg equal
        return FluidSystem::molarMass(compIdx);
    }

    /*!
     * \brief Returns the smallest activity coefficient of a component for the
     *        most current solution at a vertex.
     *
     * \param globalDofIdx The global index of the vertex (i.e. finite volume) of interest.
     * \param compIdx The index of the component of interest.
     */
    Scalar minActivityCoeff(unsigned globalDofIdx, unsigned compIdx) const
    { return minActivityCoeff_[globalDofIdx][compIdx]; }

    /*!
     * \internal
     */
    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        this->addOutputModule(new Ewoms::VtkCompositionModule<TypeTag>(this->simulator_));
        if (enableDiffusion)
            this->addOutputModule(new Ewoms::VtkDiffusionModule<TypeTag>(this->simulator_));
        if (enableEnergy)
            this->addOutputModule(new Ewoms::VtkEnergyModule<TypeTag>(this->simulator_));
    }

    mutable Scalar referencePressure_;
    mutable std::vector<ComponentVector> minActivityCoeff_;
};

} // namespace Ewoms

#endif
