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
 * \copydoc Ewoms::BlackOilModel
 */
#ifndef EWOMS_BLACK_OIL_MODEL_HH
#define EWOMS_BLACK_OIL_MODEL_HH

#include <ewoms/common/densead/math.hh>

#include "blackoilproblem.hh"
#include "blackoilindices.hh"
#include "blackoiltwophaseindices.hh"
#include "blackoilextensivequantities.hh"
#include "blackoilprimaryvariables.hh"
#include "blackoilintensivequantities.hh"
#include "blackoilratevector.hh"
#include "blackoilboundaryratevector.hh"
#include "blackoillocalresidual.hh"
#include "blackoilnewtonmethod.hh"
#include "blackoilproperties.hh"
#include "blackoilsolventmodules.hh"
#include "blackoilpolymermodules.hh"
#include "blackoilfoammodules.hh"
#include "blackoilbrinemodules.hh"
#include "blackoildarcyfluxmodule.hh"

#include <ewoms/numerics/common/multiphasebasemodel.hh>
#include <ewoms/numerics/io/vtkcompositionmodule.hh>
#include <ewoms/numerics/io/vtkblackoilmodule.hh>

#include <ewoms/material/fluidsystems/blackoilfluidsystem.hh>
#include <ewoms/common/unused.hh>
#include <ewoms/common/exceptions.hh>

#include <sstream>
#include <string>

namespace Ewoms {
template <class TypeTag>
class BlackOilModel;

template <class TypeTag>
class EclVanguard;
}

BEGIN_PROPERTIES

//! The type tag for the black-oil problems
NEW_TYPE_TAG(BlackOilModel, INHERITS_FROM(MultiPhaseBaseModel,
                                          VtkBlackOil,
                                          VtkBlackOilSolvent,
                                          VtkBlackOilPolymer,
                                          VtkBlackOilEnergy,
                                          VtkComposition));

//! Set the local residual function
SET_TYPE_PROP(BlackOilModel, LocalResidual,
              Ewoms::BlackOilLocalResidual<TypeTag>);

//! Use the black-oil specific newton method
SET_TYPE_PROP(BlackOilModel, NewtonMethod, Ewoms::BlackOilNewtonMethod<TypeTag>);

//! The Model property
SET_TYPE_PROP(BlackOilModel, Model, Ewoms::BlackOilModel<TypeTag>);

//! The Problem property
SET_TYPE_PROP(BlackOilModel, BaseProblem, Ewoms::BlackOilProblem<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(BlackOilModel, RateVector, Ewoms::BlackOilRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(BlackOilModel, BoundaryRateVector, Ewoms::BlackOilBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(BlackOilModel, PrimaryVariables, Ewoms::BlackOilPrimaryVariables<TypeTag>);

//! the IntensiveQuantities property
SET_TYPE_PROP(BlackOilModel, IntensiveQuantities, Ewoms::BlackOilIntensiveQuantities<TypeTag>);

//! the ExtensiveQuantities property
SET_TYPE_PROP(BlackOilModel, ExtensiveQuantities, Ewoms::BlackOilExtensiveQuantities<TypeTag>);

//! Use the the velocity module which is aware of the black-oil specific model extensions
//! (i.e., the polymer and solvent extensions)
SET_TYPE_PROP(BlackOilModel, FluxModule, Ewoms::BlackOilDarcyFluxModule<TypeTag>);

//! The indices required by the model
SET_TYPE_PROP(BlackOilModel, Indices,
              Ewoms::BlackOilIndices<GET_PROP_VALUE(TypeTag, EnableSolvent),
                                   GET_PROP_VALUE(TypeTag, EnablePolymer),
                                   GET_PROP_VALUE(TypeTag, EnableEnergy),
                                   GET_PROP_VALUE(TypeTag, EnableFoam),
                                   GET_PROP_VALUE(TypeTag, EnableBrine),
                                   /*PVOffset=*/0>);

//! Set the fluid system to the black-oil fluid system by default
SET_PROP(BlackOilModel, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Ewoms::BlackOilFluidSystem<Scalar> type;
};

// by default, all ECL extension modules are disabled
SET_BOOL_PROP(BlackOilModel, EnableSolvent, false);
SET_BOOL_PROP(BlackOilModel, EnablePolymer, false);
SET_BOOL_PROP(BlackOilModel, EnablePolymerMW, false);
SET_BOOL_PROP(BlackOilModel, EnableFoam, false);
SET_BOOL_PROP(BlackOilModel, EnableBrine, false);

//! By default, the blackoil model is isothermal and does not conserve energy
SET_BOOL_PROP(BlackOilModel, EnableTemperature, false);
SET_BOOL_PROP(BlackOilModel, EnableEnergy, false);

//! by default, scale the energy equation by the inverse of the energy required to heat
//! up one kg of water by 30 Kelvin. If we conserve surface volumes, this must be divided
//! by the weight of one cubic meter of water. This is required to make the "dumb" linear
//! solvers that do not weight the components of the solutions do the right thing.
//! by default, don't scale the energy equation, i.e. assume that a reasonable linear
//! solver is used. (Not scaling it makes debugging quite a bit easier.)
SET_PROP(BlackOilModel, BlackOilEnergyScalingFactor)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    static constexpr Scalar alpha = GET_PROP_VALUE(TypeTag, BlackoilConserveSurfaceVolume) ? 1000.0 : 1.0;

public:
    typedef Scalar type;
    static const Scalar value;
};

PROP_STATIC_CONST_MEMBER_DEFINITION_PREFIX_(BlackOilModel, BlackOilEnergyScalingFactor)
    ::value = 1.0/(30*4184.0*alpha);

// by default, eebos formulates the conservation equations in terms of mass not surface
// volumes
SET_BOOL_PROP(BlackOilModel, BlackoilConserveSurfaceVolume, false);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup BlackOilModel
 * \brief A fully-implicit black-oil flow model.
 *
 * The black-oil model is a three-phase, three-component model widely
 * used for oil reservoir simulation.  The phases are denoted by lower
 * index \f$\alpha \in \{ w, g, o \}\f$ ("water", "gas" and "oil") and
 * the components by upper index \f$\kappa \in \{ W, G, O \}\f$
 * ("Water", "Gas" and "Oil"). The model assumes partial miscibility:
 *
 * - Water and the gas phases are immisicible and are assumed to be
 *   only composed of the water and gas components respectively-
 * - The oil phase is assumed to be a mixture of the gas and the oil
 *  components.
 *
 * The densities of the phases are determined by so-called
 * <i>formation volume factors</i>:
 *
 * \f[
 * B_\alpha := \frac{\varrho_\alpha(1\,\text{bar})}{\varrho_\alpha(p_\alpha)}
 * \f]
 *
 * Since the gas and water phases are assumed to be immiscible, this
 * is sufficint to calculate their density. For the formation volume
 * factor of the the oil phase \f$B_o\f$ determines the density of
 * *saturated* oil, i.e. the density of the oil phase if some gas
 * phase is present.
 *
 * The composition of the oil phase is given by the <i>gas dissolution factor</i>
 * \f$R_s\f$, which defined as the volume of gas at atmospheric pressure that is
 * dissolved in a given amount of oil at reservoir pressure:
 *
 * \f[
 * R_s := \frac{\varrho_{o}^G}{\varrho_o^O}\;.
 * \f]
 *
 * This allows to calculate all quantities required for the
 * mass-conservation equations for each component, i.e.
 *
 * \f[
 * \sum_\alpha \frac{\partial\;\phi c_\alpha^\kappa S_\alpha }{\partial t}
 * - \sum_\alpha \mathrm{div} \left\{ c_\alpha^\kappa \mathbf{v}_\alpha \right\}
 * - q^\kappa = 0 \;,
 * \f]
 * where \f$\mathrm{v}_\alpha\f$ is the filter velocity of the phase
 * \f$\alpha\f$.
 *
 * By default \f$\mathrm{v}_\alpha\f$ is determined by using the
 * standard multi-phase Darcy approach, i.e.
 * \f[ \mathbf{v}_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 *\left(\mathbf{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right) \;, \f]
 * although the actual approach which is used can be specified via the
 * \c FluxModule property. For example, the velocity model can by
 * changed to the Forchheimer approach by
 * \code
 * SET_TYPE_PROP(MyProblemTypeTag, FluxModule, Ewoms::ForchheimerFluxModule<TypeTag>);
 * \endcode
 *
 * The primary variables used by this model are:
 * - The pressure of the phase with the lowest index
 * - The two saturations of the phases with the lowest indices
 */
template<class TypeTag >
class BlackOilModel
    : public MultiPhaseBaseModel<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef MultiPhaseBaseModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    enum { numComponents = FluidSystem::numComponents };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    static const bool compositionSwitchEnabled = Indices::gasEnabled;
    static const bool waterEnabled = Indices::waterEnabled;
    static const bool gasEnabled = Indices::gasEnabled;

    static const int waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static const int gasPhaseIdx = FluidSystem::gasPhaseIdx;

    static constexpr Scalar maxWaterPriority = 1.0;

    typedef BlackOilSolventModule<TypeTag> SolventModule;
    typedef BlackOilPolymerModule<TypeTag> PolymerModule;
    typedef BlackOilEnergyModule<TypeTag> EnergyModule;
public:
    BlackOilModel(Simulator& simulator)
        : ParentType(simulator)
    {
        waterPriority_ = maxWaterPriority;
    }

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        SolventModule::registerParameters();
        PolymerModule::registerParameters();
        EnergyModule::registerParameters();

        // register runtime parameters of the VTK output modules
        Ewoms::VtkBlackOilModule<TypeTag>::registerParameters();
        Ewoms::VtkCompositionModule<TypeTag>::registerParameters();
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    static std::string name()
    { return "blackoil"; }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarName
     */
    std::string primaryVarName(int pvIdx) const
    {
        std::ostringstream oss;

        if (pvIdx == Indices::waterSaturationIdx)
            oss << "saturation_" << FluidSystem::phaseName(FluidSystem::waterPhaseIdx);
        else if (pvIdx == Indices::pressureSwitchIdx)
            oss << "pressure_switching";
        else if (static_cast<int>(pvIdx) == Indices::compositionSwitchIdx)
            oss << "composition_switching";
        else if (SolventModule::primaryVarApplies(pvIdx))
            return SolventModule::primaryVarName(pvIdx);
        else if (PolymerModule::primaryVarApplies(pvIdx))
            return PolymerModule::primaryVarName(pvIdx);
        else if (EnergyModule::primaryVarApplies(pvIdx))
            return EnergyModule::primaryVarName(pvIdx);
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::eqName
     */
    std::string eqName(int eqIdx) const
    {
        std::ostringstream oss;

        if (Indices::conti0EqIdx <= eqIdx && eqIdx < Indices::conti0EqIdx + numComponents)
            oss << "conti_" << FluidSystem::phaseName(eqIdx - Indices::conti0EqIdx);
        else if (SolventModule::eqApplies(eqIdx))
            return SolventModule::eqName(eqIdx);
        else if (PolymerModule::eqApplies(eqIdx))
            return PolymerModule::eqName(eqIdx);
        else if (EnergyModule::eqApplies(eqIdx))
            return EnergyModule::eqName(eqIdx);
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarWeight
     */
    Scalar primaryVarWeight(unsigned globalDofIdx, unsigned pvIdx) const
    {
        // do not care about the auxiliary equations as they are supposed to scale
        // themselves
        if (globalDofIdx >= this->numGridDof())
            return 1.0;

        // saturations are always in the range [0, 1]!
        if (int(Indices::waterSaturationIdx) == int(pvIdx))
            return 1.0;

        // oil pressures usually are in the range of 100 to 500 bars for typical oil
        // reservoirs (which is the only relevant application for the black-oil model).
        else if (int(Indices::pressureSwitchIdx) == int(pvIdx))
            return 1.0/300e5;

        // deal with primary variables stemming from the solvent module
        else if (SolventModule::primaryVarApplies(pvIdx))
            return SolventModule::primaryVarWeight(pvIdx);

        // deal with primary variables stemming from the polymer module
        else if (PolymerModule::primaryVarApplies(pvIdx))
            return PolymerModule::primaryVarWeight(pvIdx);

        // deal with primary variables stemming from the energy module
        else if (EnergyModule::primaryVarApplies(pvIdx))
            return EnergyModule::primaryVarWeight(pvIdx);

        // if the primary variable is either the gas saturation, Rs or Rv
        assert(int(Indices::compositionSwitchIdx) == int(pvIdx));

        auto pvMeaning = this->solution(0)[globalDofIdx].primaryVarsMeaning();
        if (pvMeaning == PrimaryVariables::Sw_po_Sg)
            return 1.0; // gas saturation
        else if (pvMeaning == PrimaryVariables::Sw_po_Rs)
            return 1.0/250.; // gas dissolution factor
        else {
            assert(pvMeaning == PrimaryVariables::Sw_pg_Rv);
            return 1.0/0.025; // oil vaporization factor
        }

    }

    /*!
     * \copydoc FvBaseDiscretization::eqWeight
     */
    Scalar eqWeight(unsigned globalDofIdx, unsigned eqIdx) const
    {
        // do not care about the auxiliary equations as they are supposed to scale
        // themselves
        if (globalDofIdx >= this->numGridDof())
            return 1.0;

        if (GET_PROP_VALUE(TypeTag, BlackoilConserveSurfaceVolume)) {
            // Roughly convert the surface volume of the fluids from m^3 to kg. (in this
            // context, it does not really matter if the actual densities are off by a
            // factor of two or three.)
            switch (eqIdx) {
            case Indices::conti0EqIdx + FluidSystem::waterCompIdx:
                return 1000.0*waterPriority_;

            case Indices::conti0EqIdx + FluidSystem::gasCompIdx:
                return 1.0;

            case Indices::conti0EqIdx + FluidSystem::oilCompIdx:
                return 650.0;
            }
        }

        if (SolventModule::eqApplies(eqIdx))
            return SolventModule::eqWeight(eqIdx);

        else if (PolymerModule::eqApplies(eqIdx))
            return PolymerModule::eqWeight(eqIdx);

        else if (EnergyModule::eqApplies(eqIdx))
            return EnergyModule::eqWeight(eqIdx);

        // it is said that all kilograms are born equal (except water)!
        if (eqIdx == Indices::conti0EqIdx + FluidSystem::waterCompIdx)
            return waterPriority_;
        return 1.0;
    }

    /*!
     * \brief Write the current solution for a degree of freedom to a
     *        restart file.
     *
     * \param outstream The stream into which the vertex data should
     *                  be serialized to
     * \param dof The Dune entity which's data should be serialized
     */
    template <class DofEntity>
    void serializeEntity(std::ostream& outstream, const DofEntity& dof)
    {
        unsigned dofIdx = static_cast<unsigned>(asImp_().dofMapper().index(dof));

        // write phase state
        if (!outstream.good())
            throw std::runtime_error("Could not serialize degree of freedom "+std::to_string(dofIdx));

        // write the primary variables
        const auto& priVars = this->solution(/*timeIdx=*/0)[dofIdx];
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx)
            outstream << priVars[eqIdx] << " ";

        // write the pseudo primary variables
        outstream << priVars.primaryVarsMeaning() << " ";
        outstream << priVars.pvtRegionIndex() << " ";

        SolventModule::serializeEntity(*this, outstream, dof);
        PolymerModule::serializeEntity(*this, outstream, dof);
        EnergyModule::serializeEntity(*this, outstream, dof);
    }

    /*!
     * \brief Reads the current solution variables for a degree of
     *        freedom from a restart file.
     *
     * \param instream The stream from which the vertex data should
     *                  be deserialized from
     * \param dof The Dune entity which's data should be deserialized
     */
    template <class DofEntity>
    void deserializeEntity(std::istream& instream,
                           const DofEntity& dof)
    {
        unsigned dofIdx = static_cast<unsigned>(asImp_().dofMapper().index(dof));

        // read in the "real" primary variables of the DOF
        auto& priVars = this->solution(/*timeIdx=*/0)[dofIdx];
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (!instream.good())
                throw std::runtime_error("Could not deserialize degree of freedom "+std::to_string(dofIdx));
            instream >> priVars[eqIdx];
        }

        // read the pseudo primary variables
        unsigned primaryVarsMeaning;
        instream >> primaryVarsMeaning;

        unsigned pvtRegionIdx;
        instream >> pvtRegionIdx;

        if (!instream.good())
            throw std::runtime_error("Could not deserialize degree of freedom "+std::to_string(dofIdx));

        SolventModule::deserializeEntity(*this, instream, dof);
        PolymerModule::deserializeEntity(*this, instream, dof);
        EnergyModule::deserializeEntity(*this, instream, dof);

        typedef typename PrimaryVariables::PrimaryVarsMeaning PVM;
        priVars.setPrimaryVarsMeaning(static_cast<PVM>(primaryVarsMeaning));
        priVars.setPvtRegionIndex(pvtRegionIdx);
    }

    /*!
     * \brief Calculate the priority of the water component relative to the other
     *        components.
     *
     * This assumes that we're not that much interested in conservation of water. It can
     * considerably improve numerical performance, but it is strongly dependent on the
     * actual reservoir conditions if this is a good idea.
     *
     * Note that the only purpose of this method is to improve performance, i.e., calling
     * it is completely optional from a correctness point of view.
     */
    void updateWaterPriority()
    {
        if (!gasEnabled || !waterEnabled)
            // there is little point in determining a water priority if either water or
            // gas are not active...
            return;

        Scalar sumPv = 0.0;
        Scalar sumGasPv = 0.0;

        ElementContext elemCtx(this->simulator_);
        auto elemIt = this->gridView().template begin</*codim=*/0>();
        auto elemEndIt = this->gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++ elemIt) {
            const auto& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                // ignore non-interior elements
                continue;

            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            for (size_t dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                const auto& iq = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);

                Scalar dofPv = Ewoms::decay<Scalar>(iq.porosity())*elemCtx.dofVolume(dofIdx, /*timeIdx=*/0);
                sumPv += dofPv;
                sumGasPv += dofPv*Ewoms::decay<Scalar>(iq.fluidState().saturation(gasPhaseIdx));
            }
        }

        // take the other processes into account
        sumPv = this->gridView().comm().sum(sumPv);
        sumGasPv = this->gridView().comm().sum(sumGasPv);

        // linear interpolation with crash barriers
        Scalar x = sumGasPv/sumPv;
        Scalar b = maxWaterPriority;
        Scalar m = (0.001 - b)/0.01;

        waterPriority_ = m*x + b;
        waterPriority_ = std::max(0.001, std::min(1.0, waterPriority_));
    }

    /*!
     * \brief Deserializes the state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void deserialize(Restarter& res)
    {
        ParentType::deserialize(res);

        // set the PVT indices of the primary variables. This is also done by writing
        // them into the restart file and re-reading them, but it is better to calculate
        // them from scratch because the input could have been changed in this regard...
        ElementContext elemCtx(this->simulator_);
        auto elemIt = this->gridView().template begin</*codim=*/0>();
        auto elemEndIt = this->gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++ elemIt) {
            elemCtx.updateStencil(*elemIt);
            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timIdx=*/0); ++dofIdx) {
                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timIdx=*/0);
                updatePvtRegionIndex_(this->solution(/*timeIdx=*/0)[globalDofIdx],
                                      elemCtx,
                                      dofIdx,
                                      /*timeIdx=*/0);
            }
        }

        this->solution(/*timeIdx=*/1) = this->solution(/*timeIdx=*/0);
    }

/*
    // hack: this interferes with the static polymorphism trick
protected:
    friend ParentType;
    friend Discretization;
*/

    template <class Context>
    void supplementInitialSolution_(PrimaryVariables& priVars,
                                    const Context& context,
                                    unsigned dofIdx,
                                    unsigned timeIdx)
    { updatePvtRegionIndex_(priVars, context, dofIdx, timeIdx); }

    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        // add the VTK output modules which make sense for the blackoil model
        SolventModule::registerOutputModules(*this, this->simulator_);
        PolymerModule::registerOutputModules(*this, this->simulator_);
        EnergyModule::registerOutputModules(*this, this->simulator_);

        this->addOutputModule(new Ewoms::VtkBlackOilModule<TypeTag>(this->simulator_));
        this->addOutputModule(new Ewoms::VtkCompositionModule<TypeTag>(this->simulator_));
    }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    template <class Context>
    void updatePvtRegionIndex_(PrimaryVariables& priVars,
                               const Context& context,
                               unsigned dofIdx,
                               unsigned timeIdx)
    {
        unsigned regionIdx = context.problem().pvtRegionIndex(context, dofIdx, timeIdx);
        priVars.setPvtRegionIndex(regionIdx);
    }

    Scalar waterPriority_;
};
} // namespace Ewoms

#endif
