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
 * \copydoc Ewoms::Co2InjectionFlash
 */
#ifndef EWOMS_CO2_INJECTION_FLASH_HH
#define EWOMS_CO2_INJECTION_FLASH_HH

#include <ewoms/material/constraintsolvers/ncpflash.hh>

namespace Ewoms {
/*!
 * \brief Flash solver used by the CO2 injection problem.
 *
 * This class is just the NCP flash solver with the guessInitial()
 * method that is adapted to the pressure regime of the CO2 injection
 * problem.
 */
template <class Scalar, class FluidSystem>
class Co2InjectionFlash : public Ewoms::NcpFlash<Scalar, FluidSystem>
{
    typedef Ewoms::NcpFlash<Scalar, FluidSystem> ParentType;

    enum { numPhases = FluidSystem::numPhases };

public:
    /*!
     * \brief Guess initial values for all quantities.
     */
    template <class FluidState, class ComponentVector>
    static void guessInitial(FluidState& fluidState, const ComponentVector& globalMolarities)
    {
        ParentType::guessInitial(fluidState, globalMolarities);

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // pressure. use something close to the reservoir pressure as initial guess
            fluidState.setPressure(phaseIdx, 100e5);
        }
    }
};

} // namespace Ewoms

#endif // EWOMS_CO2_INJECTION_FLASH_HH
