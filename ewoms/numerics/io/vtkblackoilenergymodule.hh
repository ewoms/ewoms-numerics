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
 * \copydoc Ewoms::VtkBlackOilEnergyModule
 */
#ifndef EWOMS_VTK_BLACK_OIL_ENERGY_MODULE_HH
#define EWOMS_VTK_BLACK_OIL_ENERGY_MODULE_HH

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
NEW_TYPE_TAG(VtkBlackOilEnergy);

// create the property tags needed for the energy module
NEW_PROP_TAG(VtkWriteRockInternalEnergy);
NEW_PROP_TAG(VtkWriteTotalThermalConductivity);
NEW_PROP_TAG(VtkWriteFluidInternalEnergies);
NEW_PROP_TAG(VtkWriteFluidEnthalpies);
NEW_PROP_TAG(VtkOutputFormat);
NEW_PROP_TAG(EnableVtkOutput);

// set default values for what quantities to output
SET_BOOL_PROP(VtkBlackOilEnergy, VtkWriteRockInternalEnergy, true);
SET_BOOL_PROP(VtkBlackOilEnergy, VtkWriteTotalThermalConductivity, true);
SET_BOOL_PROP(VtkBlackOilEnergy, VtkWriteFluidInternalEnergies, true);
SET_BOOL_PROP(VtkBlackOilEnergy, VtkWriteFluidEnthalpies, true);

END_PROPERTIES

namespace Ewoms {
/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the black oil model's energy related quantities.
 */
template <class TypeTag>
class VtkBlackOilEnergyModule : public BaseOutputModule<TypeTag>
{
    typedef BaseOutputModule<TypeTag> ParentType;

    typedef GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    static const int vtkFormat = GET_PROP_VALUE(TypeTag, VtkOutputFormat);
    typedef Ewoms::VtkMultiWriter<GridView, vtkFormat> VtkMultiWriter;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef typename ParentType::ScalarBuffer ScalarBuffer;
    typedef typename ParentType::PhaseBuffer PhaseBuffer;

public:
    VtkBlackOilEnergyModule(const Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output
     * module.
     */
    static void registerParameters()
    {
        if (!enableEnergy)
            return;

        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteRockInternalEnergy,
                             "Include the volumetric internal energy of rock "
                             "in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteTotalThermalConductivity,
                             "Include the total thermal conductivity of the medium and the fluids "
                             "in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFluidInternalEnergies,
                             "Include the internal energies of the fluids "
                             "in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFluidEnthalpies,
                             "Include the enthalpies of the fluids "
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

        if (!enableEnergy)
            return;

        if (rockInternalEnergyOutput_())
            this->resizeScalarBuffer_(rockInternalEnergy_);
        if (totalThermalConductivityOutput_())
            this->resizeScalarBuffer_(totalThermalConductivity_);
        if (fluidInternalEnergiesOutput_())
            this->resizePhaseBuffer_(fluidInternalEnergies_);
        if (fluidEnthalpiesOutput_())
            this->resizePhaseBuffer_(fluidEnthalpies_);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
            return;

        if (!enableEnergy)
            return;

        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
            unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

            if (rockInternalEnergyOutput_())
                rockInternalEnergy_[globalDofIdx] =
                    Ewoms::scalarValue(intQuants.rockInternalEnergy());

            if (totalThermalConductivityOutput_())
                totalThermalConductivity_[globalDofIdx] =
                    Ewoms::scalarValue(intQuants.totalThermalConductivity());

            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                if (fluidInternalEnergiesOutput_())
                    fluidInternalEnergies_[phaseIdx][globalDofIdx] =
                        Ewoms::scalarValue(intQuants.fluidState().internalEnergy(phaseIdx));

                if (fluidEnthalpiesOutput_())
                    fluidEnthalpies_[phaseIdx][globalDofIdx] =
                        Ewoms::scalarValue(intQuants.fluidState().enthalpy(phaseIdx));
            }
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

        if (!enableEnergy)
            return;

        if (rockInternalEnergyOutput_())
            this->commitScalarBuffer_(baseWriter, "volumetric internal energy rock", rockInternalEnergy_);

        if (totalThermalConductivityOutput_())
            this->commitScalarBuffer_(baseWriter, "total thermal conductivity", totalThermalConductivity_);

        if (fluidInternalEnergiesOutput_())
            this->commitPhaseBuffer_(baseWriter, "internal energy_%s", fluidInternalEnergies_);

        if (fluidEnthalpiesOutput_())
            this->commitPhaseBuffer_(baseWriter, "enthalpy_%s", fluidEnthalpies_);
    }

private:
    static bool rockInternalEnergyOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteRockInternalEnergy);
        return val;
    }

    static bool totalThermalConductivityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteTotalThermalConductivity);
        return val;
    }

    static bool fluidInternalEnergiesOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFluidInternalEnergies);
        return val;
    }

    static bool fluidEnthalpiesOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFluidEnthalpies);
        return val;
    }

    ScalarBuffer rockInternalEnergy_;
    ScalarBuffer totalThermalConductivity_;
    PhaseBuffer fluidInternalEnergies_;
    PhaseBuffer fluidEnthalpies_;
};
} // namespace Ewoms

#endif
