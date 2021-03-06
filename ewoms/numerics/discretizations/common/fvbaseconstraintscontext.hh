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
 * \copydoc Ewoms::FvBaseConstraintsContext
 */
#ifndef EWOMS_FV_BASE_CONSTRAINTS_CONTEXT_HH
#define EWOMS_FV_BASE_CONSTRAINTS_CONTEXT_HH

#include "fvbaseproperties.hh"

#include <dune/common/fvector.hh>

namespace Ewoms {

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief Represents all quantities which available for calculating constraints
 */
template<class TypeTag>
class FvBaseConstraintsContext
{
    using Problem = GET_PROP_TYPE(TypeTag, Problem);
    using Model = GET_PROP_TYPE(TypeTag, Model);
    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);
    using GridView = GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    enum { dimWorld = GridView::dimensionworld };

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

public:
    /*!
     * \brief The constructor.
     */
    explicit FvBaseConstraintsContext(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    { }

    /*!
     * \copydoc Ewoms::ElementContext::problem()
     */
    const Problem& problem() const
    { return elemCtx_.problem(); }

    /*!
     * \copydoc Ewoms::ElementContext::model()
     */
    const Model& model() const
    { return elemCtx_.model(); }

    /*!
     * \copydoc Ewoms::ElementContext::gridView()
     */
    const GridView& gridView() const
    { return elemCtx_.gridView(); }

    /*!
     * \copydoc Ewoms::ElementContext::element()
     */
    const Element& element() const
    { return elemCtx_.element(); }

    /*!
     * \copydoc Ewoms::ElementContext::numDof()
     */
    int numDof(int timeIdx) const
    { return elemCtx_.numDof(timeIdx); }

    /*!
     * \copydoc Ewoms::ElementContext::numInteriorFaces()
     */
    int numInteriorFaces(int timeIdx) const
    { return elemCtx_.numInteriorFaces(timeIdx); }

    /*!
     * \copydoc Ewoms::ElementContext::globalSpaceIndex
     */
    int globalSpaceIndex(int dofIdx, int timeIdx) const
    { return elemCtx_.globalSpaceIndex(dofIdx, timeIdx); }

    /*!
     * \copydoc Ewoms::ElementContext::pos
     */
    GlobalPosition pos(int dofIdx, int timeIdx) const
    { return elemCtx_.pos(dofIdx, timeIdx); }

protected:
    const ElementContext& elemCtx_;
};

} // namespace Ewoms

#endif
