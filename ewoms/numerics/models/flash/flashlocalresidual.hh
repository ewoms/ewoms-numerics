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
 * \copydoc Ewoms::FlashLocalResidual
 */
#ifndef EWOMS_FLASH_LOCAL_RESIDUAL_HH
#define EWOMS_FLASH_LOCAL_RESIDUAL_HH

#include "flashproperties.hh"

#include <ewoms/numerics/common/diffusionmodule.hh>
#include <ewoms/numerics/common/energymodule.hh>

#include <ewoms/common/valgrind.hh>

namespace Ewoms {
/*!
 * \ingroup FlashModel
 *
 * \brief Calculates the local residual of the compositional multi-phase
 *        model based on flash calculations.
 */
template <class TypeTag>
class FlashLocalResidual: public GET_PROP_TYPE(TypeTag, DiscLocalResidual)
{
    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);
    using RateVector = GET_PROP_TYPE(TypeTag, RateVector);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);
    using IntensiveQuantities = GET_PROP_TYPE(TypeTag, IntensiveQuantities);
    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { conti0EqIdx = Indices::conti0EqIdx };

    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    using DiffusionModule = Ewoms::DiffusionModule<TypeTag, enableDiffusion>;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    using EnergyModule = Ewoms::EnergyModule<TypeTag, enableEnergy>;

    using Toolbox = Ewoms::MathToolbox<Evaluation>;

public:
    /*!
     * \copydoc ImmiscibleLocalResidual::addPhaseStorage
     */
    template <class LhsEval>
    void addPhaseStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                         const ElementContext& elemCtx,
                         unsigned dofIdx,
                         unsigned timeIdx,
                         unsigned phaseIdx) const
    {
        const IntensiveQuantities& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto& fs = intQuants.fluidState();

        // compute storage term of all components within all phases
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            unsigned eqIdx = conti0EqIdx + compIdx;
            storage[eqIdx] +=
                Toolbox::template decay<LhsEval>(fs.molarity(phaseIdx, compIdx))
                * Toolbox::template decay<LhsEval>(fs.saturation(phaseIdx))
                * Toolbox::template decay<LhsEval>(intQuants.porosity());
        }

        EnergyModule::addPhaseStorage(storage, elemCtx.intensiveQuantities(dofIdx, timeIdx), phaseIdx);
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
        storage = 0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            addPhaseStorage(storage, elemCtx, dofIdx, timeIdx, phaseIdx);

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
        addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Ewoms::Valgrind::CheckDefined(flux);

        addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Ewoms::Valgrind::CheckDefined(flux);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::addAdvectiveFlux
     */
    void addAdvectiveFlux(RateVector& flux,
                          const ElementContext& elemCtx,
                          unsigned scvfIdx,
                          unsigned timeIdx) const
    {
        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        unsigned focusDofIdx = elemCtx.focusDofIndex();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // data attached to upstream and the finite volume of the current phase
            unsigned upIdx = static_cast<unsigned>(extQuants.upstreamIndex(phaseIdx));
            const IntensiveQuantities& up = elemCtx.intensiveQuantities(upIdx, timeIdx);

            // this is a bit hacky because it is specific to the element-centered
            // finite volume scheme. (N.B. that if finite differences are used to
            // linearize the system of equations, it does not matter.)
            if (upIdx == focusDofIdx) {
                Evaluation tmp =
                    up.fluidState().molarDensity(phaseIdx)
                    * extQuants.volumeFlux(phaseIdx);

                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    flux[conti0EqIdx + compIdx] +=
                        tmp*up.fluidState().moleFraction(phaseIdx, compIdx);
                }
            }
            else {
                Evaluation tmp =
                    Toolbox::value(up.fluidState().molarDensity(phaseIdx))
                    * extQuants.volumeFlux(phaseIdx);

                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    flux[conti0EqIdx + compIdx] +=
                        tmp*Toolbox::value(up.fluidState().moleFraction(phaseIdx, compIdx));
                }
            }
        }

        EnergyModule::addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::addDiffusiveFlux
     */
    void addDiffusiveFlux(RateVector& flux,
                          const ElementContext& elemCtx,
                          unsigned scvfIdx,
                          unsigned timeIdx) const
    {
        DiffusionModule::addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        EnergyModule::addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::computeSource
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
};

} // namespace Ewoms

#endif
