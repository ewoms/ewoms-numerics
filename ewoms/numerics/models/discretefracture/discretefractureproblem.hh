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
 * \copydoc Ewoms::DiscreteFractureProblem
 */
#ifndef EWOMS_DISCRETE_FRACTURE_PROBLEM_HH
#define EWOMS_DISCRETE_FRACTURE_PROBLEM_HH

#include "discretefractureproperties.hh"

#include <ewoms/numerics/common/multiphasebaseproblem.hh>

#include <ewoms/common/means.hh>
#include <ewoms/common/unused.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

BEGIN_PROPERTIES

NEW_PROP_TAG(ThermalConductionLawParams);
NEW_PROP_TAG(EnableGravity);
NEW_PROP_TAG(FluxModule);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup DiscreteFractureModel
 * \brief The base class for the problems of ECFV discretizations which deal
 *        with a multi-phase flow through a porous medium.
 */
template<class TypeTag>
class DiscreteFractureProblem
    : public MultiPhaseBaseProblem<TypeTag>
{
    using ParentType = Ewoms::MultiPhaseBaseProblem<TypeTag>;

    using Implementation = GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = GET_PROP_TYPE(TypeTag, GridView);
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);

    enum { dimWorld = GridView::dimensionworld };
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    /*!
     * \copydoc Problem::FvBaseProblem(Simulator& )
     */
    DiscreteFractureProblem(Simulator& simulator)
        : ParentType(simulator)
    {}

    /*!
     * \brief Returns the intrinsic permeability of a face due to a fracture.
     *
     * This method is specific to the finite volume discretizations. If left unspecified,
     * it calls the intrinsicPermeability() methods for the face's interior and exterior
     * finite volume and averages them harmonically. Note that if this function is
     * defined, the intrinsicPermeability() method does not need to be defined by the
     * problem (if a finite-volume discretization is used).
     */
    template <class Context>
    void fractureFaceIntrinsicPermeability(DimMatrix& result,
                                           const Context& context,
                                           unsigned localFaceIdx,
                                           unsigned timeIdx) const
    {
        const auto& scvf = context.stencil(timeIdx).interiorFace(localFaceIdx);
        unsigned interiorElemIdx = scvf.interiorIndex();
        unsigned exteriorElemIdx = scvf.exteriorIndex();
        const DimMatrix& K1 = asImp_().fractureIntrinsicPermeability(context, interiorElemIdx, timeIdx);
        const DimMatrix& K2 = asImp_().fractureIntrinsicPermeability(context, exteriorElemIdx, timeIdx);

        // entry-wise harmonic mean. this is almost certainly wrong if
        // you have off-main diagonal entries in your permeabilities!
        for (unsigned i = 0; i < dimWorld; ++i)
            for (unsigned j = 0; j < dimWorld; ++j)
                result[i][j] = Ewoms::harmonicMean(K1[i][j], K2[i][j]);
    }
    /*!
     * \brief Returns the intrinsic permeability tensor \f$[m^2]\f$ at a given position due to a fracture
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    const DimMatrix& fractureIntrinsicPermeability(const Context& context EWOMS_UNUSED,
                                                   unsigned spaceIdx EWOMS_UNUSED,
                                                   unsigned timeIdx EWOMS_UNUSED) const
    {
        throw std::logic_error("Not implemented: Problem::fractureIntrinsicPermeability()");
    }

    /*!
     * \brief Returns the porosity [] inside fractures for a given control volume.
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    Scalar fracturePorosity(const Context& context EWOMS_UNUSED,
                            unsigned spaceIdx EWOMS_UNUSED,
                            unsigned timeIdx EWOMS_UNUSED) const
    {
        throw std::logic_error("Not implemented: Problem::fracturePorosity()");
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }
    //! \copydoc asImp_()
    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // namespace Ewoms

#endif
