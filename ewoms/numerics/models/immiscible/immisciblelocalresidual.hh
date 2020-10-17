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
 * \copydoc Ewoms::ImmiscibleLocalResidual
 */
#ifndef EWOMS_IMMISCIBLE_LOCAL_RESIDUAL_BASE_HH
#define EWOMS_IMMISCIBLE_LOCAL_RESIDUAL_BASE_HH

#include "immiscibleproperties.hh"

#include <ewoms/numerics/common/energymodule.hh>

#include <ewoms/common/valgrind.hh>

namespace Ewoms {
/*!
 * \ingroup ImmiscibleModel
 *
 * \brief Calculates the local residual of the immiscible multi-phase
 *        model.
 */
template <class TypeTag>
class ImmiscibleLocalResidual : public GET_PROP_TYPE(TypeTag, DiscLocalResidual)
{
    using Implementation = GET_PROP_TYPE(TypeTag, LocalResidual);

    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);
    using IntensiveQuantities = GET_PROP_TYPE(TypeTag, IntensiveQuantities);
    using ExtensiveQuantities = GET_PROP_TYPE(TypeTag, ExtensiveQuantities);
    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);
    using RateVector = GET_PROP_TYPE(TypeTag, RateVector);

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    using EnergyModule = Ewoms::EnergyModule<TypeTag, enableEnergy>;
    using Toolbox = Ewoms::MathToolbox<Evaluation>;

public:
    /*!
     * \brief Adds the amount all conservation quantities (e.g. phase
     *        mass) within a single fluid phase
     *
     * \copydetails Doxygen::storageParam
     * \copydetails Doxygen::dofCtxParams
     * \copydetails Doxygen::phaseIdxParam
     */
    template <class LhsEval>
    void addPhaseStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                         const ElementContext& elemCtx,
                         unsigned dofIdx,
                         unsigned timeIdx,
                         unsigned phaseIdx) const
    {
        // retrieve the intensive quantities for the SCV at the specified
        // point in time
        const IntensiveQuantities& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto& fs = intQuants.fluidState();

        storage[conti0EqIdx + phaseIdx] =
            Toolbox::template decay<LhsEval>(intQuants.porosity())
            * Toolbox::template decay<LhsEval>(fs.saturation(phaseIdx))
            * Toolbox::template decay<LhsEval>(fs.density(phaseIdx));

        EnergyModule::addPhaseStorage(storage, intQuants, phaseIdx);
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeStorage
     */
    template <class LhsEval>
    void computeStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                        const ElementContext& elemCtx,
                        unsigned dofIdx,
                        unsigned timeIdx) const
    {
        storage = 0.0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            asImp_().addPhaseStorage(storage, elemCtx, dofIdx, timeIdx, phaseIdx);

        EnergyModule::addSolidEnergyStorage(storage, elemCtx.intensiveQuantities(dofIdx, timeIdx));
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeFlux
     */
    void computeFlux(RateVector& flux,
                     const ElementContext& elemCtx,
                     unsigned scvfIdx,
                     unsigned timeIdx) const
    {
        flux = 0.0;
        asImp_().addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        asImp_().addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \brief Add the advective mass flux at a given flux integration point
     *
     * \copydetails computeFlux
     */
    void addAdvectiveFlux(RateVector& flux,
                          const ElementContext& elemCtx,
                          unsigned scvfIdx,
                          unsigned timeIdx) const
    {
        const ExtensiveQuantities& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        ////////
        // advective fluxes of all components in all phases
        ////////
        unsigned focusDofIdx = elemCtx.focusDofIndex();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // data attached to upstream DOF of the current phase.
            unsigned upIdx = static_cast<unsigned>(extQuants.upstreamIndex(phaseIdx));

            const IntensiveQuantities& up = elemCtx.intensiveQuantities(upIdx, /*timeIdx=*/0);

            // add advective flux of current component in current phase.
            const Evaluation& rho = up.fluidState().density(phaseIdx);
            if (focusDofIdx == upIdx)
                flux[conti0EqIdx + phaseIdx] += extQuants.volumeFlux(phaseIdx)*rho;
            else
                flux[conti0EqIdx + phaseIdx] += extQuants.volumeFlux(phaseIdx)*Toolbox::value(rho);
        }

        EnergyModule::addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \brief Adds the diffusive flux at a given flux integration point.
     *
     * For the immiscible model, this is a no-op for mass fluxes. For energy it adds the
     * contribution of thermal conduction to the enthalpy flux.
     *
     * \copydetails computeFlux
     */
    void addDiffusiveFlux(RateVector& flux,
                          const ElementContext& elemCtx,
                          unsigned scvfIdx,
                          unsigned timeIdx) const
    {
        // no diffusive mass fluxes for the immiscible model

        // thermal conduction
        EnergyModule::addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeSource
     *
     * By default, this method only asks the problem to specify a
     * source term.
     */
    void computeSource(RateVector& source,
                       const ElementContext& elemCtx,
                       unsigned dofIdx,
                       unsigned timeIdx) const
    {
        Ewoms::Valgrind::SetUndefined(source);
        elemCtx.problem().source(source, elemCtx, dofIdx, timeIdx);
        Ewoms::Valgrind::CheckDefined(source);
    }

private:
    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // namespace Ewoms

#endif
