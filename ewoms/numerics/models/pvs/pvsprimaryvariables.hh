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
 * \copydoc Ewoms::PvsPrimaryVariables
 */
#ifndef EWOMS_PVS_PRIMARY_VARIABLES_HH
#define EWOMS_PVS_PRIMARY_VARIABLES_HH

#include "pvsindices.hh"
#include "pvsproperties.hh"

#include <ewoms/numerics/discretizations/common/fvbaseprimaryvariables.hh>
#include <ewoms/numerics/common/energymodule.hh>

#include <ewoms/material/constraintsolvers/ncpflash.hh>
#include <ewoms/material/fluidstates/compositionalfluidstate.hh>
#include <ewoms/common/valgrind.hh>
#include <ewoms/common/exceptions.hh>

#include <dune/common/fvector.hh>

#include <iostream>

namespace Ewoms {

/*!
 * \ingroup PvsModel
 *
 * \brief Represents the primary variables used in the primary
 *        variable switching compositional model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
template <class TypeTag>
class PvsPrimaryVariables : public FvBasePrimaryVariables<TypeTag>
{
    using ParentType = FvBasePrimaryVariables<TypeTag>;
    using ThisType = PvsPrimaryVariables<TypeTag>;
    using Implementation = GET_PROP_TYPE(TypeTag, PrimaryVariables);

    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
    using MaterialLaw = GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = GET_PROP_TYPE(TypeTag, MaterialLawParams);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);

    // primary variable indices
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { switch0Idx = Indices::switch0Idx };

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    using Toolbox = typename Ewoms::MathToolbox<Evaluation>;
    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;
    using EnergyModule = Ewoms::EnergyModule<TypeTag, enableEnergy>;
    using NcpFlash = Ewoms::NcpFlash<Scalar, FluidSystem>;

public:
    PvsPrimaryVariables() : ParentType()
    { Ewoms::Valgrind::SetDefined(*this); }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
     */
    explicit PvsPrimaryVariables(Scalar value) : ParentType(value)
    {
        Ewoms::Valgrind::CheckDefined(value);
        Ewoms::Valgrind::SetDefined(*this);

        phasePresence_ = 0;
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const
     * ImmisciblePrimaryVariables& )
     */
    PvsPrimaryVariables(const PvsPrimaryVariables& value) : ParentType(value)
    {
        Ewoms::Valgrind::SetDefined(*this);

        phasePresence_ = value.phasePresence_;
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignMassConservative
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
        Ewoms::CompositionalFluidState<Scalar, FluidSystem> fsFlash;

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
        NcpFlash::template solve<MaterialLaw>(fsFlash, matParams, paramCache, globalMolarities);

        // use the result to assign the primary variables
        assignNaive(fsFlash);
    }

    /*!
     * \brief Return the fluid phases which are present in a given
     *        control volume.
     */
    short phasePresence() const
    { return phasePresence_; }

    /*!
     * \brief Set which fluid phases are present in a given control volume.
     *
     * \param value The new phase presence. The phase with index i is
     *              present if the i-th bit of \c value is 1.
     */
    void setPhasePresence(short value)
    { phasePresence_ = value; }

    /*!
     * \brief Set whether a given indivividual phase should be present
     *        or not.
     *
     * \param phaseIdx The index of the phase which's presence ought to be set or reset.
     * \param yesno If true, the presence of the phase is set, else it is reset
     */
    void setPhasePresent(unsigned phaseIdx, bool yesno = true)
    {
        if (yesno)
            setPhasePresence(phasePresence_ | (1 << phaseIdx));
        else
            setPhasePresence(phasePresence_&  ~(1 << phaseIdx));
    }

    /*!
     * \brief Returns the index of the phase with's its saturation is
     *        determined by the closure condition of saturation.
     */
    unsigned implicitSaturationIdx() const
    { return lowestPresentPhaseIdx(); }

    /*!
     * \brief Returns true iff a phase is present for a given phase
     *        presence.
     *
     * \param phaseIdx The index of the phase which's presence is
     *                 queried.
     * \param phasePresence The bit-map of present phases.
     */
    static bool phaseIsPresent(unsigned phaseIdx, short phasePresence)
    { return phasePresence&  (1 << phaseIdx); }

    /*!
     * \brief Returns true iff a phase is present for the current
     *        phase presence.
     *
     * \copydoc Doxygen::phaseIdxParam
     */
    bool phaseIsPresent(unsigned phaseIdx) const
    { return phasePresence_&  (1 << phaseIdx); }

    /*!
     * \brief Returns the phase with the lowest index that is present.
     */
    unsigned lowestPresentPhaseIdx() const
    { return static_cast<unsigned>(ffs(phasePresence_) - 1); }

    /*!
     * \brief Assignment operator from an other primary variables object
     */
    ThisType& operator=(const Implementation& value)
    {
        ParentType::operator=(value);
        phasePresence_ = value.phasePresence_;

        return *this;
    }

    /*!
     * \brief Assignment operator from a scalar value
     */
    ThisType& operator=(Scalar value)
    {
        ParentType::operator=(value);

        phasePresence_ = 0;
        return *this;
    }

    /*!
     * \brief Returns an explcitly stored saturation for a given phase.
     *
     * (or 0 if the saturation is not explicitly stored.)
     *
     * \copydoc Doxygen::phaseIdxParam
     */
    Evaluation explicitSaturationValue(unsigned phaseIdx, unsigned timeIdx) const
    {
        if (!phaseIsPresent(phaseIdx) || phaseIdx == lowestPresentPhaseIdx())
            // non-present phases have saturation 0
            return 0.0;

        unsigned varIdx = switch0Idx + phaseIdx - 1;
        if (std::is_same<Evaluation, Scalar>::value)
            return (*this)[varIdx]; // finite differences
        else {
            // automatic differentiation
            if (timeIdx != 0)
                Toolbox::createConstant((*this)[varIdx]);
            return Toolbox::createVariable((*this)[varIdx], varIdx);
        }
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignNaive
     */
    template <class FluidState>
    void assignNaive(const FluidState& fluidState)
    {
        using FsToolbox = Ewoms::MathToolbox<typename FluidState::Scalar>;

        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(*this, fluidState);

        // set the pressure of the first phase
        (*this)[pressure0Idx] = FsToolbox::value(fluidState.pressure(/*phaseIdx=*/0));
        Ewoms::Valgrind::CheckDefined((*this)[pressure0Idx]);

        // determine the phase presence.
        phasePresence_ = 0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // use a NCP condition to determine if the phase is
            // present or not
            Scalar a = 1;
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                a -= FsToolbox::value(fluidState.moleFraction(phaseIdx, compIdx));
            }
            Scalar b = FsToolbox::value(fluidState.saturation(phaseIdx));

            if (b > a)
                phasePresence_ |= (1 << phaseIdx);
        }

        // some phase must be present
        if (phasePresence_ == 0)
            throw Ewoms::NumericalIssue("Phase state was 0, i.e., no fluid is present");

        // set the primary variables which correspond to mole
        // fractions of the present phase which has the lowest index.
        unsigned lowestPhaseIdx = lowestPresentPhaseIdx();
        for (unsigned switchIdx = 0; switchIdx < numPhases - 1; ++switchIdx) {
            unsigned phaseIdx = switchIdx;
            unsigned compIdx = switchIdx + 1;
            if (switchIdx >= lowestPhaseIdx)
                ++phaseIdx;

            if (phaseIsPresent(phaseIdx)) {
                (*this)[switch0Idx + switchIdx] = FsToolbox::value(fluidState.saturation(phaseIdx));
                Ewoms::Valgrind::CheckDefined((*this)[switch0Idx + switchIdx]);
            }
            else {
                (*this)[switch0Idx + switchIdx] =
                    FsToolbox::value(fluidState.moleFraction(lowestPhaseIdx, compIdx));
                Ewoms::Valgrind::CheckDefined((*this)[switch0Idx + switchIdx]);
            }
        }

        // set the mole fractions in of the remaining components in
        // the phase with the lowest index
        for (unsigned compIdx = numPhases - 1; compIdx < numComponents - 1; ++compIdx) {
            (*this)[switch0Idx + compIdx] =
                FsToolbox::value(fluidState.moleFraction(lowestPhaseIdx, compIdx + 1));
            Ewoms::Valgrind::CheckDefined((*this)[switch0Idx + compIdx]);
        }
    }

    /*!
     * \copydoc FlashPrimaryVariables::print
     */
    void print(std::ostream& os = std::cout) const
    {
        os << "(p_" << FluidSystem::phaseName(0) << " = "
           << this->operator[](pressure0Idx);
        unsigned lowestPhaseIdx = lowestPresentPhaseIdx();
        for (unsigned switchIdx = 0; switchIdx < numPhases - 1; ++switchIdx) {
            unsigned phaseIdx = switchIdx;
            unsigned compIdx = switchIdx + 1;
            if (phaseIdx >= lowestPhaseIdx)
                ++phaseIdx; // skip the saturation of the present
                            // phase with the lowest index

            if (phaseIsPresent(phaseIdx)) {
                os << ", S_" << FluidSystem::phaseName(phaseIdx) << " = "
                   << (*this)[switch0Idx + switchIdx];
            }
            else {
                os << ", x_" << FluidSystem::phaseName(lowestPhaseIdx) << "^"
                   << FluidSystem::componentName(compIdx) << " = "
                   << (*this)[switch0Idx + switchIdx];
            }
        }
        for (unsigned compIdx = numPhases - 1; compIdx < numComponents - 1;
             ++compIdx) {
            os << ", x_" << FluidSystem::phaseName(lowestPhaseIdx) << "^"
               << FluidSystem::componentName(compIdx + 1) << " = "
               << (*this)[switch0Idx + compIdx];
        }
        os << ")";
        os << ", phase presence: " << static_cast<int>(phasePresence_) << std::flush;
    }

private:
    short phasePresence_;
};

} // namespace Ewoms

#endif
