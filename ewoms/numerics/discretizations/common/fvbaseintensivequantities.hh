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
 * \copydoc Ewoms::FvBaseIntensiveQuantities
 */
#ifndef EWOMS_FV_BASE_INTENSIVE_QUANTITIES_HH
#define EWOMS_FV_BASE_INTENSIVE_QUANTITIES_HH

#include "fvbaseproperties.hh"

#include <ewoms/common/valgrind.hh>
#include <ewoms/common/unused.hh>

namespace Ewoms {

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief Base class for the model specific class which provides access to all intensive
 *        (i.e., volume averaged) quantities.
 */
template <class TypeTag>
class FvBaseIntensiveQuantities
{
    using Implementation = GET_PROP_TYPE(TypeTag, IntensiveQuantities);
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);

public:
    // default constructor
    FvBaseIntensiveQuantities()
    { }

    // copy constructor
    FvBaseIntensiveQuantities(const FvBaseIntensiveQuantities& v) = default;

    /*!
     * \brief Register all run-time parameters for the intensive quantities.
     */
    static void registerParameters()
    { }

    /*!
     * \brief Update all quantities for a given control volume.
     */
    void update(const ElementContext& elemCtx,
                unsigned dofIdx,
                unsigned timeIdx)
    { extrusionFactor_ = elemCtx.problem().extrusionFactor(elemCtx, dofIdx, timeIdx); }

    /*!
     * \brief Return how much a given sub-control volume is extruded.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    /*!
     * \brief If running in valgrind this makes sure that all
     *        quantities in the intensive quantities are defined.
     */
    void checkDefined() const
    { }

private:
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Scalar extrusionFactor_;
};

} // namespace Ewoms

#endif
