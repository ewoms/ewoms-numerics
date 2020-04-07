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
 * \copydoc Ewoms::ImmisciblePrimaryVariables
 */
#ifndef EWOMS_IMMISCIBLE_PRIMARY_VARIABLES_HH
#define EWOMS_IMMISCIBLE_PRIMARY_VARIABLES_HH

#include "immiscibleproperties.hh"

#include <ewoms/numerics/discretizations/common/fvbaseprimaryvariables.hh>
#include <ewoms/numerics/common/energymodule.hh>

#include <ewoms/material/constraintsolvers/immiscibleflash.hh>
#include <ewoms/material/fluidstates/immisciblefluidstate.hh>
#include <ewoms/common/valgrind.hh>

#include <dune/common/fvector.hh>

namespace Ewoms {

/*!
 * \ingroup ImmiscibleModel
 *
 * \brief Represents the primary variables used by the immiscible
 *        multi-phase, model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
template <class TypeTag>
class ImmisciblePrimaryVariables : public FvBasePrimaryVariables<TypeTag>
{
    typedef FvBasePrimaryVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    // primary variable indices
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;
    typedef Ewoms::ImmiscibleFlash<Scalar, FluidSystem> ImmiscibleFlash;
    typedef Ewoms::EnergyModule<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy)> EnergyModule;

public:
    /*!
     * \brief Default constructor
     */
    ImmisciblePrimaryVariables() : ParentType()
    { Ewoms::Valgrind::SetUndefined(*this); }

    /*!
     * \brief Constructor with assignment from scalar
     *
     * \param value The scalar value to which all entries of the vector will be set.
     */
    ImmisciblePrimaryVariables(Scalar value) : ParentType(value)
    {}

    /*!
     * \brief Copy constructor
     *
     * \param value The primary variables that will be duplicated.
     */
    ImmisciblePrimaryVariables(const ImmisciblePrimaryVariables& value) = default;

    /*!
     * \brief Assignment operator
     *
     * \param value The primary variables that will be duplicated.
     */
    ImmisciblePrimaryVariables& operator=(const ImmisciblePrimaryVariables& value) = default;

    /*!
     * \brief Set the primary variables from an arbitrary fluid state
     *        in a mass conservative way.
     *
     * If an energy equation is included, the fluid temperatures are
     * the same as the one given in the fluid state, *not* the
     * enthalpy.
     *
     * \param fluidState The fluid state which should be represented
     *                   by the primary variables. The temperatures,
     *                   pressures, compositions and densities of all
     *                   phases must be defined.
     * \param matParams The capillary pressure law parameters
     * \param isInEquilibrium If true, the fluid state expresses
     *                        thermodynamic equilibrium assuming the
     *                        relations expressed by the fluid
     *                        system. This implies that in addition to
     *                        the quantities mentioned above, the
     *                        fugacities are also defined.
     */
    template <class FluidState>
    void assignMassConservative(const FluidState& fluidState,
                                const MaterialLawParams& matParams,
                                bool isInEquilibrium = false)
    {
        #ifndef NDEBUG
        // make sure the temperature is the same in all fluid phases
        for (unsigned phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx) {
            assert(std::abs(fluidState.temperature(0) - fluidState.temperature(phaseIdx)) < 1e-30);
        }
#endif // NDEBUG

        // for the equilibrium case, we don't need complicated
        // computations.
        if (isInEquilibrium) {
            assignNaive(fluidState);
            return;
        }

        // use a flash calculation to calculate a fluid state in
        // thermodynamic equilibrium
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        Ewoms::ImmiscibleFluidState<Scalar, FluidSystem> fsFlash;

        // use the externally given fluid state as initial value for
        // the flash calculation
        fsFlash.assign(fluidState);

        // calculate the phase densities
        paramCache.updateAll(fsFlash);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar rho = FluidSystem::density(fsFlash, paramCache, phaseIdx);
            fsFlash.setDensity(phaseIdx, rho);
        }

        // calculate the "global molarities"
        ComponentVector globalMolarities(0.0);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                globalMolarities[compIdx] +=
                    fsFlash.saturation(phaseIdx) * fsFlash.molarity(phaseIdx, compIdx);
            }
        }

        // run the flash calculation
        ImmiscibleFlash::template solve<MaterialLaw>(fsFlash, matParams, paramCache, globalMolarities);

        // use the result to assign the primary variables
        assignNaive(fsFlash);
    }

    /*!
     * \brief Directly retrieve the primary variables from an
     *        arbitrary fluid state.
     *
     * This method retrieves all primary variables from an abitrary
     * fluid state without careing whether the state which is
     * represented by the resulting primary variables features the
     * equivalent mass as the given fluid state. This method is
     * massively cheaper and simpler than assignMassConservative() but
     * it should be used with care!
     *
     * \param fluidState The fluid state which should be represented
     *                   by the primary variables. The temperatures,
     *                   pressures, compositions and densities of all
     *                   phases must be defined.
     */
    template <class FluidState>
    void assignNaive(const FluidState& fluidState)
    {
        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(asImp_(), fluidState);

        (*this)[pressure0Idx] = fluidState.pressure(/*phaseIdx=*/0);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            (*this)[saturation0Idx + phaseIdx] = fluidState.saturation(phaseIdx);
    }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }
};

} // namespace Ewoms

#endif
