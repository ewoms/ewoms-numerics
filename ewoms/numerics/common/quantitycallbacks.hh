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
 * \brief This method contains all callback classes for quantities
 *        that are required by some extensive quantities
 */
#ifndef EWOMS_QUANTITY_CALLBACKS_HH
#define EWOMS_QUANTITY_CALLBACKS_HH

#include <ewoms/numerics/discretizations/common/fvbaseproperties.hh>

#include <ewoms/common/mathtoolbox.hh>
#include <ewoms/common/valgrind.hh>

#include <type_traits>
#include <utility>

namespace Ewoms {
/*!
 * \ingroup Discretization
 *
 * \brief Callback class for temperature.
 */
template <class TypeTag>
class TemperatureCallback
{
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;

    typedef decltype(std::declval<IntensiveQuantities>().fluidState()) IQFluidState;
    typedef decltype(std::declval<IQFluidState>().temperature(0)) ResultRawType;

public:
    typedef typename std::remove_const<typename std::remove_reference<ResultRawType>::type>::type ResultType;
    typedef typename Ewoms::MathToolbox<ResultType>::ValueType ResultValueType;

    TemperatureCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    {}

    /*!
     * \brief Return the temperature given the index of a degree of freedom within an
     *        element context.
     *
     * In this context, we assume that thermal equilibrium applies, i.e. that the
     * temperature of all phases is equal.
     */
    ResultType operator()(unsigned dofIdx) const
    { return elemCtx_.intensiveQuantities(dofIdx, /*timeIdx=*/0).fluidState().temperature(/*phaseIdx=*/0); }

private:
    const ElementContext& elemCtx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for a phase pressure.
 */
template <class TypeTag>
class PressureCallback
{
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;

    typedef decltype(std::declval<IntensiveQuantities>().fluidState()) IQFluidState;
    typedef decltype(std::declval<IQFluidState>().pressure(0)) ResultRawType;

public:
    typedef typename std::remove_const<typename std::remove_reference<ResultRawType>::type>::type ResultType;
    typedef typename Ewoms::MathToolbox<ResultType>::ValueType ResultValueType;

    PressureCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    { Ewoms::Valgrind::SetUndefined(phaseIdx_); }

    PressureCallback(const ElementContext& elemCtx, unsigned phaseIdx)
        : elemCtx_(elemCtx)
        , phaseIdx_(static_cast<unsigned short>(phaseIdx))
    {}

    /*!
     * \brief Set the index of the fluid phase for which the pressure
     *        should be returned.
     */
    void setPhaseIndex(unsigned phaseIdx)
    { phaseIdx_ = static_cast<unsigned short>(phaseIdx); }

    /*!
     * \brief Return the pressure of the specified phase given the index of a degree of
     *        freedom within an element context.
     */
    ResultType operator()(unsigned dofIdx) const
    {
        Ewoms::Valgrind::CheckDefined(phaseIdx_);
        return elemCtx_.intensiveQuantities(dofIdx, /*timeIdx=*/0).fluidState().pressure(phaseIdx_);
    }

private:
    const ElementContext& elemCtx_;
    unsigned short phaseIdx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for a phase pressure.
 */
template <class TypeTag, class FluidState>
class BoundaryPressureCallback
{
    typedef GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;

    typedef decltype(std::declval<IntensiveQuantities>().fluidState()) IQRawFluidState;
    typedef typename std::remove_const<typename std::remove_reference<IQRawFluidState>::type>::type IQFluidState;
    typedef typename IQFluidState::Scalar IQScalar;

public:
    typedef IQScalar ResultType;

    BoundaryPressureCallback(const ElementContext& elemCtx, const FluidState& boundaryFs)
        : elemCtx_(elemCtx)
        , boundaryFs_(boundaryFs)
    { Ewoms::Valgrind::SetUndefined(phaseIdx_); }

    BoundaryPressureCallback(const ElementContext& elemCtx,
                             const FluidState& boundaryFs,
                             unsigned phaseIdx)
        : elemCtx_(elemCtx)
        , boundaryFs_(boundaryFs)
        , phaseIdx_(static_cast<unsigned short>(phaseIdx))
    {}

    /*!
     * \brief Set the index of the fluid phase for which the pressure
     *        should be returned.
     */
    void setPhaseIndex(unsigned phaseIdx)
    { phaseIdx_ = static_cast<unsigned short>(phaseIdx); }

    /*!
     * \brief Return the pressure of a phase given the index of a
     *        degree of freedom within an element context.
     */
    ResultType operator()(unsigned dofIdx) const
    {
        Ewoms::Valgrind::CheckDefined(phaseIdx_);
        return elemCtx_.intensiveQuantities(dofIdx, /*timeIdx=*/0).fluidState().pressure(phaseIdx_);
    }

    IQScalar boundaryValue() const
    {
        Ewoms::Valgrind::CheckDefined(phaseIdx_);
        return boundaryFs_.pressure(phaseIdx_);
    }

private:
    const ElementContext& elemCtx_;
    const FluidState& boundaryFs_;
    unsigned short phaseIdx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for the density of a phase.
 */
template <class TypeTag>
class DensityCallback
{
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;

    typedef decltype(std::declval<IntensiveQuantities>().fluidState()) IQFluidState;
    typedef decltype(std::declval<IQFluidState>().density(0)) ResultRawType;

public:
    typedef typename std::remove_const<typename std::remove_reference<ResultRawType>::type>::type ResultType;
    typedef typename Ewoms::MathToolbox<ResultType>::ValueType ResultValueType;

    DensityCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    { Ewoms::Valgrind::SetUndefined(phaseIdx_); }

    DensityCallback(const ElementContext& elemCtx, unsigned phaseIdx)
        : elemCtx_(elemCtx)
        , phaseIdx_(static_cast<unsigned short>(phaseIdx))
    {}

    /*!
     * \brief Set the index of the fluid phase for which the density
     *        should be returned.
     */
    void setPhaseIndex(unsigned phaseIdx)
    { phaseIdx_ = static_cast<unsigned short>(phaseIdx); }

    /*!
     * \brief Return the density of a phase given the index of a
     *        degree of freedom within an element context.
     */
    ResultType operator()(unsigned dofIdx) const
    {
        Ewoms::Valgrind::CheckDefined(phaseIdx_);
        return elemCtx_.intensiveQuantities(dofIdx, /*timeIdx=*/0).fluidState().density(phaseIdx_);
    }

private:
    const ElementContext& elemCtx_;
    unsigned short phaseIdx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for the molar density of a phase.
 */
template <class TypeTag>
class MolarDensityCallback
{
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;

    typedef decltype(std::declval<IntensiveQuantities>().fluidState()) IQFluidState;

public:
    typedef decltype(std::declval<IQFluidState>().molarDensity(0)) ResultType;
    typedef typename Ewoms::MathToolbox<ResultType>::ValueType ResultValueType;

    MolarDensityCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    { Ewoms::Valgrind::SetUndefined(phaseIdx_); }

    MolarDensityCallback(const ElementContext& elemCtx, unsigned phaseIdx)
        : elemCtx_(elemCtx)
        , phaseIdx_(static_cast<unsigned short>(phaseIdx))
    {}

    /*!
     * \brief Set the index of the fluid phase for which the molar
     *        density should be returned.
     */
    void setPhaseIndex(unsigned phaseIdx)
    { phaseIdx_ = static_cast<unsigned short>(phaseIdx); }

    /*!
     * \brief Return the molar density of a phase given the index of a
     *        degree of freedom within an element context.
     */
    ResultType operator()(unsigned dofIdx) const
    {
        Ewoms::Valgrind::CheckDefined(phaseIdx_);
        return elemCtx_.intensiveQuantities(dofIdx, /*timeIdx=*/0).fluidState().molarDensity(phaseIdx_);
    }

private:
    const ElementContext& elemCtx_;
    unsigned short phaseIdx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for the viscosity of a phase.
 */
template <class TypeTag>
class ViscosityCallback
{
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;

    typedef decltype(std::declval<IntensiveQuantities>().fluidState()) IQFluidState;
    typedef decltype(std::declval<IQFluidState>().viscosity(0)) ResultRawType;

public:
    typedef typename std::remove_const<typename std::remove_reference<ResultRawType>::type>::type ResultType;
    typedef typename Ewoms::MathToolbox<ResultType>::ValueType ResultValueType;

    ViscosityCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    { Ewoms::Valgrind::SetUndefined(phaseIdx_); }

    ViscosityCallback(const ElementContext& elemCtx, unsigned phaseIdx)
        : elemCtx_(elemCtx)
        , phaseIdx_(static_cast<unsigned short>(phaseIdx))
    {}

    /*!
     * \brief Set the index of the fluid phase for which the viscosity
     *        should be returned.
     */
    void setPhaseIndex(unsigned phaseIdx)
    { phaseIdx_ = static_cast<unsigned short>(phaseIdx); }

    /*!
     * \brief Return the viscosity of a phase given the index of a
     *        degree of freedom within an element context.
     */
    ResultType operator()(unsigned dofIdx) const
    {
        Ewoms::Valgrind::CheckDefined(phaseIdx_);
        return elemCtx_.intensiveQuantities(dofIdx, /*timeIdx=*/0).fluidState().viscosity(phaseIdx_);
    }

private:
    const ElementContext& elemCtx_;
    unsigned short phaseIdx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for the velocity of a phase at the center of a DOF.
 */
template <class TypeTag>
class VelocityCallback
{
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;

    typedef decltype(IntensiveQuantities().velocityCenter()) ResultRawType;


public:
    typedef typename std::remove_const<typename std::remove_reference<ResultRawType>::type>::type ResultType;
    typedef typename ResultType::field_type ResultFieldType;
    typedef typename Ewoms::MathToolbox<ResultFieldType>::ValueType ResultFieldValueType;

    VelocityCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    {}

    /*!
     * \brief Return the velocity of a phase given the index of a
     *        degree of freedom within an element context.
     */
    ResultType operator()(unsigned dofIdx) const
    { return elemCtx_.intensiveQuantities(dofIdx, /*timeIdx=*/0).velocityCenter(); }

private:
    const ElementContext& elemCtx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for the velocity of a phase at the center of a DOF.
 */
template <class TypeTag>
class VelocityComponentCallback
{
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;

    typedef decltype(IntensiveQuantities().velocityCenter()[0]) ResultRawType;

public:
    typedef typename std::remove_const<typename std::remove_reference<ResultRawType>::type>::type ResultType;
    typedef typename Ewoms::MathToolbox<ResultType>::ValueType ResultValueType;

    VelocityComponentCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    { Ewoms::Valgrind::SetUndefined(dimIdx_); }

    VelocityComponentCallback(const ElementContext& elemCtx, unsigned dimIdx)
        : elemCtx_(elemCtx)
        , dimIdx_(dimIdx)
    {}

    /*!
     * \brief Set the index of the component of the velocity
     *        which should be returned.
     */
    void setDimIndex(unsigned dimIdx)
    { dimIdx_ = dimIdx; }

    /*!
     * \brief Return the velocity of a phase given the index of a
     *        degree of freedom within an element context.
     */
    ResultType operator()(unsigned dofIdx) const
    {
        Ewoms::Valgrind::CheckDefined(dimIdx_);
        return elemCtx_.intensiveQuantities(dofIdx, /*timeIdx=*/0).velocityCenter()[dimIdx_];
    }

private:
    const ElementContext& elemCtx_;
    unsigned dimIdx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for a mole fraction of a component in a phase.
 */
template <class TypeTag>
class MoleFractionCallback
{
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;

    typedef decltype(std::declval<IntensiveQuantities>().fluidState()) IQFluidState;
    typedef decltype(std::declval<IQFluidState>().moleFraction(0, 0)) ResultRawType;

public:
    typedef typename std::remove_const<typename std::remove_reference<ResultRawType>::type>::type ResultType;
    typedef typename Ewoms::MathToolbox<ResultType>::ValueType ResultValueType;

    MoleFractionCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    {
        Ewoms::Valgrind::SetUndefined(phaseIdx_);
        Ewoms::Valgrind::SetUndefined(compIdx_);
    }

    MoleFractionCallback(const ElementContext& elemCtx, unsigned phaseIdx, unsigned compIdx)
        : elemCtx_(elemCtx)
        , phaseIdx_(static_cast<unsigned short>(phaseIdx))
        , compIdx_(static_cast<unsigned short>(compIdx))
    {}

    /*!
     * \brief Set the index of the fluid phase for which a mole fraction should be
     *        returned.
     */
    void setPhaseIndex(unsigned phaseIdx)
    { phaseIdx_ = static_cast<unsigned short>(phaseIdx); }

    /*!
     * \brief Set the index of the component for which the mole fraction should be
     *        returned.
     */
    void setComponentIndex(unsigned compIdx)
    { compIdx_ = static_cast<unsigned short>(compIdx); }

    /*!
     * \brief Return the mole fraction of a component in a phase given the index of a
     *        degree of freedom within an element context.
     */
    ResultType operator()(unsigned dofIdx) const
    {
        Ewoms::Valgrind::CheckDefined(phaseIdx_);
        Ewoms::Valgrind::CheckDefined(compIdx_);
        return elemCtx_.intensiveQuantities(dofIdx, /*timeIdx=*/0).fluidState().moleFraction(phaseIdx_, compIdx_);
    }

private:
    const ElementContext& elemCtx_;
    unsigned short phaseIdx_;
    unsigned short compIdx_;
};

} // namespace Ewoms

#endif
