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
 * \copydoc Ewoms::ImmiscibleModel
 */
#ifndef EWOMS_IMMISCIBLE_MODEL_HH
#define EWOMS_IMMISCIBLE_MODEL_HH

#include <ewoms/common/densead/math.hh>
#include "immiscibleproperties.hh"
#include "immiscibleindices.hh"
#include "immiscibleextensivequantities.hh"
#include "immiscibleprimaryvariables.hh"
#include "immiscibleintensivequantities.hh"
#include "immiscibleratevector.hh"
#include "immiscibleboundaryratevector.hh"
#include "immisciblelocalresidual.hh"

#include <ewoms/numerics/common/multiphasebasemodel.hh>
#include <ewoms/numerics/common/energymodule.hh>
#include <ewoms/numerics/io/vtkenergymodule.hh>
#include <ewoms/material/components/nullcomponent.hh>
#include <ewoms/material/fluidsystems/gasphase.hh>
#include <ewoms/material/fluidsystems/liquidphase.hh>
#include <ewoms/material/fluidsystems/singlephasefluidsystem.hh>
#include <ewoms/material/fluidsystems/twophaseimmisciblefluidsystem.hh>

#include <sstream>
#include <string>

namespace Ewoms {
template <class TypeTag>
class ImmiscibleModel;
}

BEGIN_PROPERTIES

//! The generic type tag for problems using the immiscible multi-phase model
NEW_TYPE_TAG(ImmiscibleModel, INHERITS_FROM(MultiPhaseBaseModel, VtkEnergy));
//! The type tag for single-phase immiscible problems
NEW_TYPE_TAG(ImmiscibleSinglePhaseModel, INHERITS_FROM(ImmiscibleModel));
//! The type tag for two-phase immiscible problems
NEW_TYPE_TAG(ImmiscibleTwoPhaseModel, INHERITS_FROM(ImmiscibleModel));

//! Use the immiscible multi-phase local jacobian operator for the immiscible multi-phase model
SET_TYPE_PROP(ImmiscibleModel, LocalResidual,
              Ewoms::ImmiscibleLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(ImmiscibleModel, Model, Ewoms::ImmiscibleModel<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(ImmiscibleModel, RateVector, Ewoms::ImmiscibleRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(ImmiscibleModel, BoundaryRateVector, Ewoms::ImmiscibleBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(ImmiscibleModel, PrimaryVariables, Ewoms::ImmisciblePrimaryVariables<TypeTag>);

//! the IntensiveQuantities property
SET_TYPE_PROP(ImmiscibleModel, IntensiveQuantities, Ewoms::ImmiscibleIntensiveQuantities<TypeTag>);

//! the ExtensiveQuantities property
SET_TYPE_PROP(ImmiscibleModel, ExtensiveQuantities, Ewoms::ImmiscibleExtensiveQuantities<TypeTag>);

//! The indices required by the isothermal immiscible multi-phase model
SET_TYPE_PROP(ImmiscibleModel, Indices, Ewoms::ImmiscibleIndices<TypeTag, /*PVOffset=*/0>);

//! Disable the energy equation by default
SET_BOOL_PROP(ImmiscibleModel, EnableEnergy, false);

/////////////////////
// set slightly different properties for the single-phase case
/////////////////////

//! The fluid system to use by default
SET_PROP(ImmiscibleSinglePhaseModel, FluidSystem)
{ private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Fluid = GET_PROP_TYPE(TypeTag, Fluid);
public:
    using type = Ewoms::SinglePhaseFluidSystem<Scalar , Fluid>;
};

SET_PROP(ImmiscibleSinglePhaseModel, Fluid)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);

public:
    using type = Ewoms::LiquidPhase<Scalar, Ewoms::NullComponent<Scalar> >;
};

// disable output of a few quantities which make sense in a
// multi-phase but not in a single-phase context
SET_BOOL_PROP(ImmiscibleSinglePhaseModel, VtkWriteSaturations, false);
SET_BOOL_PROP(ImmiscibleSinglePhaseModel, VtkWriteMobilities, false);
SET_BOOL_PROP(ImmiscibleSinglePhaseModel, VtkWriteRelativePermeabilities, false);

/////////////////////
// set slightly different properties for the two-phase case
/////////////////////
SET_PROP(ImmiscibleTwoPhaseModel, WettingPhase)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);

public:
    using type = Ewoms::LiquidPhase<Scalar, Ewoms::NullComponent<Scalar> >;
};

SET_PROP(ImmiscibleTwoPhaseModel, NonwettingPhase)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);

public:
    using type = Ewoms::LiquidPhase<Scalar, Ewoms::NullComponent<Scalar> >;
};

SET_PROP(ImmiscibleTwoPhaseModel, FluidSystem)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using WettingPhase = GET_PROP_TYPE(TypeTag, WettingPhase);
    using NonwettingPhase = GET_PROP_TYPE(TypeTag, NonwettingPhase);

public:
    using type = Ewoms::TwoPhaseImmiscibleFluidSystem<Scalar, WettingPhase, NonwettingPhase>;
};

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup ImmiscibleModel
 * \brief A fully-implicit multi-phase flow model which assumes
 *        immiscibility of the phases.
 *
 * This model implements multi-phase flow of \f$M > 0\f$ immiscible
 * fluids \f$\alpha\f$. By default, the standard multi-phase Darcy
 * approach is used to determine the velocity, i.e.
 * \f[
 * \mathbf{v}_\alpha =
 * - \frac{k_{r\alpha}}{\mu_\alpha}
 * \mathbf{K}\left(\mathbf{grad}\, p_\alpha -
 *                 \varrho_{\alpha} \mathbf{g} \right) \;,
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
 * \frac{\partial\;\phi S_\alpha \rho_\alpha }{\partial t}
 * - \mathrm{div} \left\{ \rho_\alpha \mathbf{v}_\alpha  \right\}
 * - q_\alpha = 0 \;.
 * \f]
 *
 * The model uses the following primary variables:
 * - The pressure \f$p_0\f$ in Pascal of the phase with the lowest index
 * - The saturations \f$S_\alpha\f$ of the \f$M - 1\f$ phases that
 *   exhibit the lowest indices
 * - The absolute temperature \f$T\f$ in Kelvin if energy is conserved
 *   via the energy equation
 */
template <class TypeTag>
class ImmiscibleModel
    : public Ewoms::MultiPhaseBaseModel<TypeTag>
{
    using ParentType = Ewoms::MultiPhaseBaseModel<TypeTag>;
    using Implementation = GET_PROP_TYPE(TypeTag, Model);
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);

    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);

    enum { numComponents = FluidSystem::numComponents };

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    using EnergyModule = Ewoms::EnergyModule<TypeTag, enableEnergy>;

public:
    ImmiscibleModel(Simulator& simulator)
        : ParentType(simulator)
    {}

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        if (enableEnergy)
            Ewoms::VtkEnergyModule<TypeTag>::registerParameters();
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    static std::string name()
    { return "immiscible"; }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarName
     */
    std::string primaryVarName(unsigned pvIdx) const
    {
        std::string s;
        if (!(s = EnergyModule::primaryVarName(pvIdx)).empty())
            return s;

        std::ostringstream oss;

        if (pvIdx == Indices::pressure0Idx) {
            oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/0);
        }
        else if (Indices::saturation0Idx <= pvIdx
                 && pvIdx < Indices::saturation0Idx + numPhases - 1) {
            unsigned phaseIdx = pvIdx - Indices::saturation0Idx;
            oss << "saturation_" << FluidSystem::phaseName(phaseIdx);
        }
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

        if (Indices::conti0EqIdx <= eqIdx && eqIdx < Indices::conti0EqIdx + numComponents)
            oss << "conti_" << FluidSystem::phaseName(eqIdx - Indices::conti0EqIdx);
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
        size_t nDof = this->numTotalDof();
        for (unsigned dofIdx = 0; dofIdx < nDof; ++ dofIdx) {
            if (this->isLocalDof(dofIdx)) {
                referencePressure_ =
                    this->solution(/*timeIdx=*/0)[dofIdx][/*pvIdx=*/Indices::pressure0Idx];
                break;
            }
        }
    }

    /*!
     * \copydetails FvBaseDiscretization::primaryVarWeight
     */
    Scalar primaryVarWeight(unsigned globalDofIdx, unsigned pvIdx) const
    {
        assert(referencePressure_ > 0);

        Scalar tmp = EnergyModule::primaryVarWeight(asImp_(), globalDofIdx, pvIdx);
        if (tmp > 0)
            // energy related quantity
            return tmp;
        if (Indices::pressure0Idx == pvIdx) {
            return 10 / referencePressure_;
        }
        return 1.0;
    }

    /*!
     * \copydetails FvBaseDiscretization::eqWeight
     */
    Scalar eqWeight(unsigned globalDofIdx, unsigned eqIdx) const
    {
        Scalar tmp = EnergyModule::eqWeight(asImp_(), globalDofIdx, eqIdx);
        if (tmp > 0)
            // energy related equation
            return tmp;

#ifndef NDEBUG
        unsigned compIdx = eqIdx - Indices::conti0EqIdx;
        assert(0 <= compIdx && compIdx <= numPhases);
#endif

        // make all kg equal
        return 1.0;
    }

    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        if (enableEnergy)
            this->addOutputModule(new Ewoms::VtkEnergyModule<TypeTag>(this->simulator_));
    }

private:
    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }

    mutable Scalar referencePressure_;
};
} // namespace Ewoms

#endif
