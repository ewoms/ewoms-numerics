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
 * \copydoc Ewoms::FvBaseConstraints
 */
#ifndef EWOMS_FV_BASE_CONSTRAINTS_HH
#define EWOMS_FV_BASE_CONSTRAINTS_HH

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/valgrind.hh>

BEGIN_PROPERTIES

NEW_PROP_TAG(PrimaryVariables);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief Class to specify constraints for a finite volume spatial discretization.
 *
 * Note that eWoms does not support "partial" constraints. (I.e., either all equations
 * must be constraint or none.)
 */
template <class TypeTag>
class FvBaseConstraints : public GET_PROP_TYPE(TypeTag, PrimaryVariables)
{

public:
    FvBaseConstraints()
    { setActive(false); }

    FvBaseConstraints(const FvBaseConstraints&) = default;

//! \cond SKIP
    /*!
     * \brief Use the default assignment operator
     */
    FvBaseConstraints &operator=(const FvBaseConstraints&) = default;
//! \endcond

    /*!
     * \brief Returns true if the constraints are enforced.
     */
    bool isActive() const
    { return isActive_; }

    /*!
     * \brief Specify whether the constraints should be enforced or not
     */
    void setActive(bool yesno)
    { isActive_ = yesno; }

private:
    bool isActive_;
};

} // namespace Ewoms

#endif
