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
 * \copydoc Ewoms::VcfvDiscretization
 */
#ifndef EWOMS_VCFV_DISCRETIZATION_HH
#define EWOMS_VCFV_DISCRETIZATION_HH

#include <ewoms/common/densead/math.hh>

#include "vcfvproperties.hh"
#include "vcfvstencil.hh"
#include "p1fegradientcalculator.hh"
#include "vcfvgridcommhandlefactory.hh"
#include "vcfvbaseoutputmodule.hh"

#include <ewoms/numerics/linear/vertexborderlistfromgrid.hh>
#include <ewoms/numerics/discretizations/common/fvbasediscretization.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>
#endif

namespace Ewoms {
template <class TypeTag>
class VcfvDiscretization;

} // namespace Ewoms

BEGIN_PROPERTIES

//! Set the stencil
SET_PROP(VcfvDiscretization, Stencil)
{
private:
    using GridView = GET_PROP_TYPE(TypeTag, GridView);
    using CoordScalar = typename GridView::ctype;

public:
    using type = Ewoms::VcfvStencil<CoordScalar, GridView>;
};

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(VcfvDiscretization, DofMapper, GET_PROP_TYPE(TypeTag, VertexMapper));

//! The concrete class which manages the spatial discretization
SET_TYPE_PROP(VcfvDiscretization, Discretization, Ewoms::VcfvDiscretization<TypeTag>);

//! The base class for the output modules (decides whether to write
//! element or vertex based fields)
SET_TYPE_PROP(VcfvDiscretization, DiscBaseOutputModule,
              Ewoms::VcfvBaseOutputModule<TypeTag>);

//! Calculates the gradient of any quantity given the index of a flux approximation point
SET_TYPE_PROP(VcfvDiscretization, GradientCalculator,
              Ewoms::P1FeGradientCalculator<TypeTag>);

//! The class to create grid communication handles
SET_TYPE_PROP(VcfvDiscretization, GridCommHandleFactory,
              Ewoms::VcfvGridCommHandleFactory<TypeTag>);

//! Use two-point gradients by default for the vertex centered finite volume scheme.
SET_BOOL_PROP(VcfvDiscretization, UseP1FiniteElementGradients, false);

#if HAVE_DUNE_FEM
//! Set the DiscreteFunctionSpace
SET_PROP(VcfvDiscretization, DiscreteFunctionSpace)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar)  ;
    using GridPart = GET_PROP_TYPE(TypeTag, GridPart);
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::Fem::FunctionSpace<typename GridPart::GridType::ctype,
                                     Scalar,
                                     GridPart::GridType::dimensionworld,
                                     numEq> FunctionSpace;
public:
    // Lagrange discrete function space with unknowns at the cell vertices
    using type = Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, 1 >;
};
#endif

//! Set the border list creator for vertices
SET_PROP(VcfvDiscretization, BorderListCreator)
{ private:
    using VertexMapper = GET_PROP_TYPE(TypeTag, VertexMapper);
    using GridView = GET_PROP_TYPE(TypeTag, GridView);
public:
    using type = Ewoms::Linear::VertexBorderListFromGrid<GridView, VertexMapper>;
};

//! For the vertex centered finite volume method, ghost and overlap elements must _not_
//! be assembled to avoid accounting twice for the fluxes over the process boundary faces
//! of the local process' grid partition
SET_BOOL_PROP(VcfvDiscretization, LinearizeNonLocalElements, false);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup VcfvDiscretization
 *
 * \brief The base class for the vertex centered finite volume discretization scheme.
 */
template<class TypeTag>
class VcfvDiscretization : public FvBaseDiscretization<TypeTag>
{
    using ParentType = FvBaseDiscretization<TypeTag>;
    using Implementation = GET_PROP_TYPE(TypeTag, Model);
    using DofMapper = GET_PROP_TYPE(TypeTag, DofMapper);
    using GridView = GET_PROP_TYPE(TypeTag, GridView);
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);

    enum { dim = GridView::dimension };

public:
    VcfvDiscretization(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Returns a string of discretization's human-readable name
     */
    static std::string discretizationName()
    { return "vcfv"; }

    /*!
     * \brief Returns the number of global degrees of freedom (DOFs) due to the grid
     */
    size_t numGridDof() const
    { return static_cast<size_t>(this->gridView_.size(/*codim=*/dim)); }

    /*!
     * \brief Mapper to convert the Dune entities of the
     *        discretization's degrees of freedoms are to indices.
     */
    const DofMapper& dofMapper() const
    { return this->vertexMapper(); }

    /*!
     * \brief Serializes the current state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void serialize(Restarter& res)
    { res.template serializeEntities</*codim=*/dim>(asImp_(), this->gridView_); }

    /*!
     * \brief Deserializes the state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void deserialize(Restarter& res)
    {
        res.template deserializeEntities</*codim=*/dim>(asImp_(), this->gridView_);
        this->solution(/*timeIdx=*/1) = this->solution(/*timeIdx=*/0);
    }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }
};
} // namespace Ewoms

#endif
