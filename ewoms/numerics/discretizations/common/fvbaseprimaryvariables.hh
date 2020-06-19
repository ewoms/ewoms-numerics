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
 * \copydoc Ewoms::FvBasePrimaryVariables
 */
#ifndef EWOMS_FV_BASE_PRIMARY_VARIABLES_HH
#define EWOMS_FV_BASE_PRIMARY_VARIABLES_HH

#include <type_traits>

#include "fvbaseproperties.hh"

#include <ewoms/common/valgrind.hh>
#include <ewoms/common/unused.hh>
#include <ewoms/common/exceptions.hh>

#include <dune/common/fvector.hh>

namespace Ewoms {

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief Represents the primary variables used by the a model.
 */
template <class TypeTag>
class FvBasePrimaryVariables
    : public Dune::FieldVector<GET_PROP_TYPE(TypeTag, Scalar),
                               GET_PROP_VALUE(TypeTag, NumEq)>
{
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    using Toolbox = Ewoms::MathToolbox<Evaluation>;
    using ParentType = Dune::FieldVector<Scalar, numEq>;

public:
    FvBasePrimaryVariables()
        : ParentType()
    { Ewoms::Valgrind::SetUndefined(*this); }

    /*!
     * \brief Construction from a scalar value
     */
    FvBasePrimaryVariables(Scalar value)
        : ParentType(value)
    { }

    /*!
     * \brief Assignment from another primary variables object
     */
    FvBasePrimaryVariables(const FvBasePrimaryVariables& value) = default;

    /*!
     * \brief Assignment from another primary variables object
     */
    FvBasePrimaryVariables& operator=(const FvBasePrimaryVariables& value) = default;

    /*!
     * \brief Return a primary variable intensive evaluation.
     *
     * i.e., the result represents the function f = x_i if the time index is zero, else
     * it represents the a constant f = x_i. (the difference is that in the first case,
     * the derivative w.r.t. x_i is 1, while it is 0 in the second case.
     */
    Evaluation makeEvaluation(unsigned varIdx, unsigned timeIdx) const
    {
        if (std::is_same<Evaluation, Scalar>::value)
            return (*this)[varIdx]; // finite differences
        else {
            // automatic differentiation
            if (timeIdx == 0)
                return Toolbox::createVariable((*this)[varIdx], varIdx);
            else
                return Toolbox::createConstant((*this)[varIdx]);
        }
    }

    /*!
     * \brief Assign the primary variables "somehow" from a fluid state
     *
     * That is without considering any consistency issues which the
     * fluid state might have. This method is guaranteed to produce
     * consistent results if the fluid state is consistent to the
     * properties at a given spatial location. (Where "consistent
     * results" means that the same fluid state can be reconstructed
     * from the primary variables.)
     */
    template <class FluidState>
    void assignNaive(const FluidState& fluidState EWOMS_UNUSED)
    {
        throw std::runtime_error("The PrimaryVariables class does not define "
                                 "an assignNaive() method");
    }

    /*!
     * \brief Instruct valgrind to check the definedness of all attributes of this class.
     */
    void checkDefined() const
    {
        Ewoms::Valgrind::CheckDefined(*static_cast<const ParentType*>(this));
    }
};

} // namespace Ewoms

namespace Dune {

  /** Compatibility traits class for DenseVector and DenseMatrix.
   */
  template<class TypeTag, bool>
  struct FieldTraitsImpl;

  /** FieldTraitsImpl for classes derived from
   * Ewoms::FvBasePrimaryVariables: use FieldVector's FieldTraits implementation) */
  template<class TypeTag>
  struct FieldTraitsImpl< TypeTag, true >
      : public FieldTraits<FieldVector<GET_PROP_TYPE(TypeTag, Scalar),
                                       GET_PROP_VALUE(TypeTag, NumEq)> >
  {
  };

  /** FieldTraitsImpl for classes not derived from
   * Ewoms::FvBasePrimaryVariables, fall bakc to existing implementation */
  template<class T>
  struct FieldTraitsImpl< T, false >
    : public FieldTraits< T >
  {
  };

  /** Specialization of FieldTraits for all PrimaryVariables derived from Ewoms::FvBasePrimaryVariables */
  template<class TypeTag, template <class> class EwomsPrimaryVariable>
  struct FieldTraits< EwomsPrimaryVariable< TypeTag > >
    : public FieldTraitsImpl< TypeTag,
                              std::is_base_of< Ewoms::FvBasePrimaryVariables< TypeTag >,
                                               EwomsPrimaryVariable< TypeTag > > :: value >
  {
  };
}

#endif
