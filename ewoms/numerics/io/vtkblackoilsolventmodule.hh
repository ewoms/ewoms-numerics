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
 * \copydoc Ewoms::VtkBlackOilSolventModule
 */
#ifndef EWOMS_VTK_BLACK_OIL_SOLVENT_MODULE_HH
#define EWOMS_VTK_BLACK_OIL_SOLVENT_MODULE_HH

#include <ewoms/common/densead/math.hh>

#include "vtkmultiwriter.hh"
#include "baseoutputmodule.hh"

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>
#include <ewoms/numerics/models/blackoil/blackoilproperties.hh>

#include <dune/common/fvector.hh>

#include <cstdio>

BEGIN_PROPERTIES

// create new type tag for the VTK multi-phase output
NEW_TYPE_TAG(VtkBlackOilSolvent);

// create the property tags needed for the solvent output module
NEW_PROP_TAG(EnableSolvent);
NEW_PROP_TAG(EnableVtkOutput);
NEW_PROP_TAG(VtkWriteSolventSaturation);
NEW_PROP_TAG(VtkWriteSolventDensity);
NEW_PROP_TAG(VtkWriteSolventViscosity);
NEW_PROP_TAG(VtkWriteSolventMobility);

// set default values for what quantities to output
SET_BOOL_PROP(VtkBlackOilSolvent, VtkWriteSolventSaturation, true);
SET_BOOL_PROP(VtkBlackOilSolvent, VtkWriteSolventDensity, true);
SET_BOOL_PROP(VtkBlackOilSolvent, VtkWriteSolventViscosity, true);
SET_BOOL_PROP(VtkBlackOilSolvent, VtkWriteSolventMobility, true);

END_PROPERTIES

namespace Ewoms {
/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the black oil model's solvent related quantities.
 */
template <class TypeTag>
class VtkBlackOilSolventModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using GridView = GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);
    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);

    static const int vtkFormat = GET_PROP_VALUE(TypeTag, VtkOutputFormat);
    using VtkMultiWriter = Ewoms::VtkMultiWriter<GridView, vtkFormat>;

    enum { enableSolvent = GET_PROP_VALUE(TypeTag, EnableSolvent) };

    using ScalarBuffer = typename ParentType::ScalarBuffer;

public:
    VtkBlackOilSolventModule(const Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output
     * module.
     */
    static void registerParameters()
    {
        if (!enableSolvent)
            return;

        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteSolventSaturation,
                             "Include the \"saturation\" of the solvent component "
                             "in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteSolventDensity,
                             "Include the \"density\" of the solvent component "
                             "in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteSolventViscosity,
                             "Include the \"viscosity\" of the solvent component "
                             "in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteSolventMobility,
                             "Include the \"mobility\" of the solvent component "
                             "in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
            return;

        if (!enableSolvent)
            return;

        if (solventSaturationOutput_())
            this->resizeScalarBuffer_(solventSaturation_);
        if (solventDensityOutput_())
            this->resizeScalarBuffer_(solventDensity_);
        if (solventViscosityOutput_())
            this->resizeScalarBuffer_(solventViscosity_);
        if (solventMobilityOutput_())
            this->resizeScalarBuffer_(solventMobility_);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
            return;

        if (!enableSolvent)
            return;

        using Toolbox = Ewoms::MathToolbox<Evaluation>;
        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
            unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

            if (solventSaturationOutput_())
                solventSaturation_[globalDofIdx] =
                    Toolbox::scalarValue(intQuants.solventSaturation());

            if (solventDensityOutput_())
                solventDensity_[globalDofIdx] =
                    Toolbox::scalarValue(intQuants.solventDensity());

            if (solventViscosityOutput_())
                solventViscosity_[globalDofIdx] =
                    Toolbox::scalarValue(intQuants.solventViscosity());

            if (solventMobilityOutput_())
                solventMobility_[globalDofIdx] =
                    Toolbox::scalarValue(intQuants.solventMobility());
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter& baseWriter)
    {
        VtkMultiWriter *vtkWriter = dynamic_cast<VtkMultiWriter*>(&baseWriter);
        if (!vtkWriter)
            return;

        if (!enableSolvent)
            return;

        if (solventSaturationOutput_())
            this->commitScalarBuffer_(baseWriter, "saturation_solvent", solventSaturation_);

        if (solventDensityOutput_())
            this->commitScalarBuffer_(baseWriter, "density_solvent", solventDensity_);

        if (solventViscosityOutput_())
            this->commitScalarBuffer_(baseWriter, "viscosity_solvent", solventViscosity_);

        if (solventMobilityOutput_())
            this->commitScalarBuffer_(baseWriter, "mobility_solvent", solventMobility_);
    }

private:
    static bool solventSaturationOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteSolventSaturation);
        return val;
    }

    static bool solventDensityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteSolventDensity);
        return val;
    }

    static bool solventViscosityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteSolventViscosity);
        return val;
    }

    static bool solventMobilityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteSolventMobility);
        return val;
    }

    ScalarBuffer solventSaturation_;
    ScalarBuffer solventDensity_;
    ScalarBuffer solventViscosity_;
    ScalarBuffer solventMobility_;
};
} // namespace Ewoms

#endif
