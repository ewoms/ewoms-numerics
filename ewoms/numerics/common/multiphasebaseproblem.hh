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
 * \copydoc Ewoms::MultiPhaseBaseProblem
 */
#ifndef EWOMS_MULTI_PHASE_BASE_PROBLEM_HH
#define EWOMS_MULTI_PHASE_BASE_PROBLEM_HH

#include "multiphasebaseproperties.hh"

#include <ewoms/numerics/discretizations/common/fvbaseproblem.hh>
#include <ewoms/numerics/discretizations/common/fvbaseproperties.hh>

#include <ewoms/material/fluidmatrixinteractions/nullmaterial.hh>
#include <ewoms/common/means.hh>
#include <ewoms/common/unused.hh>
#include <ewoms/common/exceptions.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

BEGIN_PROPERTIES

NEW_PROP_TAG(SolidEnergyLawParams);
NEW_PROP_TAG(ThermalConductionLawParams);
NEW_PROP_TAG(EnableGravity);
NEW_PROP_TAG(FluxModule);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup Discretization
 *
 * \brief The base class for the problems of ECFV discretizations which deal
 *        with a multi-phase flow through a porous medium.
 */
template<class TypeTag>
class MultiPhaseBaseProblem
    : public FvBaseProblem<TypeTag>
    , public GET_PROP_TYPE(TypeTag, FluxModule)::FluxBaseProblem
{
//! \cond SKIP_THIS
    typedef Ewoms::FvBaseProblem<TypeTag> ParentType;

    typedef GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef GET_PROP_TYPE(TypeTag, SolidEnergyLawParams) SolidEnergyLawParams;
    typedef GET_PROP_TYPE(TypeTag, ThermalConductionLawParams) ThermalConductionLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params MaterialLawParams;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
//! \endcond

public:
    /*!
     * \copydoc Problem::FvBaseProblem(Simulator& )
     */
    MultiPhaseBaseProblem(Simulator& simulator)
        : ParentType(simulator)
    { init_(); }

    /*!
     * \brief Register all run-time parameters for the problem and the model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableGravity,
                             "Use the gravity correction for the pressure gradients.");
    }

    /*!
     * \brief Returns the intrinsic permeability of an intersection.
     *
     * This method is specific to the finite volume discretizations. If left unspecified,
     * it calls the intrinsicPermeability() method for the intersection's interior and
     * exterior finite volumes and averages them harmonically. Note that if this function
     * is defined, the intrinsicPermeability() method does not need to be defined by the
     * problem (if a finite-volume discretization is used).
     */
    template <class Context>
    void intersectionIntrinsicPermeability(DimMatrix& result,
                                           const Context& context,
                                           unsigned intersectionIdx,
                                           unsigned timeIdx) const
    {
        const auto& scvf = context.stencil(timeIdx).interiorFace(intersectionIdx);

        const DimMatrix& K1 = asImp_().intrinsicPermeability(context, scvf.interiorIndex(), timeIdx);
        const DimMatrix& K2 = asImp_().intrinsicPermeability(context, scvf.exteriorIndex(), timeIdx);

        // entry-wise harmonic mean. this is almost certainly wrong if
        // you have off-main diagonal entries in your permeabilities!
        for (unsigned i = 0; i < dimWorld; ++i)
            for (unsigned j = 0; j < dimWorld; ++j)
                result[i][j] = Ewoms::harmonicMean(K1[i][j], K2[i][j]);
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the intrinsic permeability tensor \f$[m^2]\f$ at a given position
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context EWOMS_UNUSED,
                                           unsigned spaceIdx EWOMS_UNUSED,
                                           unsigned timeIdx EWOMS_UNUSED) const
    {
        throw std::logic_error("Not implemented: Problem::intrinsicPermeability()");
    }

    /*!
     * \brief Returns the porosity [] of the porous medium for a given
     *        control volume.
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    Scalar porosity(const Context& context EWOMS_UNUSED,
                    unsigned spaceIdx EWOMS_UNUSED,
                    unsigned timeIdx EWOMS_UNUSED) const
    {
        throw std::logic_error("Not implemented: Problem::porosity()");
    }

    /*!
     * \brief Returns the parameter object for the energy storage law of the solid in a
     *        sub-control volume.
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    const SolidEnergyLawParams&
    solidEnergyParams(const Context& context EWOMS_UNUSED,
                      unsigned spaceIdx EWOMS_UNUSED,
                      unsigned timeIdx EWOMS_UNUSED) const
    {
        throw std::logic_error("Not implemented: Problem::solidEnergyParams()");
    }

    /*!
     * \brief Returns the parameter object for the thermal conductivity law in a
     *        sub-control volume.
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    const ThermalConductionLawParams&
    thermalConductionParams(const Context& context EWOMS_UNUSED,
                         unsigned spaceIdx EWOMS_UNUSED,
                         unsigned timeIdx EWOMS_UNUSED) const
    {
        throw std::logic_error("Not implemented: Problem::thermalConductionParams()");
    }

    /*!
     * \brief Define the tortuosity.
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    Scalar tortuosity(const Context& context EWOMS_UNUSED,
                      unsigned spaceIdx EWOMS_UNUSED,
                      unsigned timeIdx EWOMS_UNUSED) const
    {
        throw std::logic_error("Not implemented: Problem::tortuosity()");
    }

    /*!
     * \brief Define the dispersivity.
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    Scalar dispersivity(const Context& context EWOMS_UNUSED,
                        unsigned spaceIdx EWOMS_UNUSED,
                        unsigned timeIdx EWOMS_UNUSED) const
    {
        throw std::logic_error("Not implemented: Problem::dispersivity()");
    }

    /*!
     * \brief Returns the material law parameters \f$\mathrm{[K]}\f$ within a control volume.
     *
     * If you get a compiler error at this method, you set the
     * MaterialLaw property to something different than
     * Ewoms::NullMaterialLaw. In this case, you have to overload the
     * matererialLaw() method in the derived class!
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    const MaterialLawParams &
    materialLawParams(const Context& context EWOMS_UNUSED,
                      unsigned spaceIdx EWOMS_UNUSED,
                      unsigned timeIdx EWOMS_UNUSED) const
    {
        static MaterialLawParams dummy;
        return dummy;
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ within a control volume.
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    Scalar temperature(const Context& context EWOMS_UNUSED,
                       unsigned spaceIdx EWOMS_UNUSED,
                       unsigned timeIdx EWOMS_UNUSED) const
    { return asImp_().temperature(); }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is to be used.
     */
    Scalar temperature() const
    { throw std::logic_error("Not implemented:temperature() method not implemented by the actual problem"); }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    const DimVector& gravity(const Context& context EWOMS_UNUSED,
                             unsigned spaceIdx EWOMS_UNUSED,
                             unsigned timeIdx EWOMS_UNUSED) const
    { return asImp_().gravity(); }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * This method is used for problems where the gravitational
     * acceleration does not depend on the spatial position. The
     * default behaviour is that if the <tt>EnableGravity</tt>
     * property is true, \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$ holds,
     * else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$.
     */
    const DimVector& gravity() const
    { return gravity_; }

    /*!
     * \brief Mark grid cells for refinement or coarsening
     *
     * \return The number of elements marked for refinement or coarsening.
     */
    unsigned markForGridAdaptation()
    {
        typedef Ewoms::MathToolbox<Evaluation> Toolbox;

        unsigned numMarked = 0;
        ElementContext elemCtx( this->simulator() );
        auto gridView = this->simulator().vanguard().gridView();
        auto& grid = this->simulator().vanguard().grid();
        auto elemIt = gridView.template begin</*codim=*/0, Dune::Interior_Partition>();
        auto elemEndIt = gridView.template end</*codim=*/0, Dune::Interior_Partition>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            const auto& element = *elemIt ;
            elemCtx.updateAll( element );

            // HACK: this should better be part of an AdaptionCriterion class
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar minSat = 1e100 ;
                Scalar maxSat = -1e100;
                size_t nDofs = elemCtx.numDof(/*timeIdx=*/0);
                for (unsigned dofIdx = 0; dofIdx < nDofs; ++dofIdx)
                {
                    const auto& intQuant = elemCtx.intensiveQuantities( dofIdx, /*timeIdx=*/0 );
                    minSat = std::min(minSat,
                                      Toolbox::value(intQuant.fluidState().saturation(phaseIdx)));
                    maxSat = std::max(maxSat,
                                      Toolbox::value(intQuant.fluidState().saturation(phaseIdx)));
                }

                const Scalar indicator =
                    (maxSat - minSat)/(std::max<Scalar>(0.01, maxSat+minSat)/2);
                if( indicator > 0.2 && element.level() < 2 ) {
                    grid.mark( 1, element );
                    ++ numMarked;
                }
                else if ( indicator < 0.025 ) {
                    grid.mark( -1, element );
                    ++ numMarked;
                }
                else
                {
                    grid.mark( 0, element );
                }
            }
        }

        // get global sum so that every proc is on the same page
        numMarked = this->simulator().vanguard().grid().comm().sum( numMarked );

        return numMarked;
    }

    // \}

protected:
    /*!
     * \brief Converts a Scalar value to an isotropic Tensor
     *
     * This is convenient e.g. for specifying intrinsic permebilities:
     * \code{.cpp}
     * auto permTensor = this->toDimMatrix_(1e-12);
     * \endcode
     *
     * \param val The scalar value which should be expressed as a tensor
     */
    DimMatrix toDimMatrix_(Scalar val) const
    {
        DimMatrix ret(0.0);
        for (unsigned i = 0; i < DimMatrix::rows; ++i)
            ret[i][i] = val;
        return ret;
    }

    DimVector gravity_;

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }
    //! \copydoc asImp_()
    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }

    void init_()
    {
        gravity_ = 0.0;
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dimWorld-1]  = -9.81;
    }
};

} // namespace Ewoms

#endif
