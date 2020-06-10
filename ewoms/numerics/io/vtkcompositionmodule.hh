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
 * \copydoc Ewoms::VtkCompositionModule
 */
#ifndef EWOMS_VTK_COMPOSITION_MODULE_HH
#define EWOMS_VTK_COMPOSITION_MODULE_HH

#include "vtkmultiwriter.hh"
#include "baseoutputmodule.hh"

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <ewoms/common/mathtoolbox.hh>

BEGIN_PROPERTIES

// create new type tag for the VTK composition output
NEW_TYPE_TAG(VtkComposition);

// create the property tags needed for the composition module
NEW_PROP_TAG(VtkWriteMassFractions);
NEW_PROP_TAG(VtkWriteMoleFractions);
NEW_PROP_TAG(VtkWriteTotalMassFractions);
NEW_PROP_TAG(VtkWriteTotalMoleFractions);
NEW_PROP_TAG(VtkWriteMolarities);
NEW_PROP_TAG(VtkWriteFugacities);
NEW_PROP_TAG(VtkWriteFugacityCoeffs);
NEW_PROP_TAG(VtkOutputFormat);
NEW_PROP_TAG(EnableVtkOutput);

// set default values for what quantities to output
SET_BOOL_PROP(VtkComposition, VtkWriteMassFractions, false);
SET_BOOL_PROP(VtkComposition, VtkWriteMoleFractions, true);
SET_BOOL_PROP(VtkComposition, VtkWriteTotalMassFractions, false);
SET_BOOL_PROP(VtkComposition, VtkWriteTotalMoleFractions, false);
SET_BOOL_PROP(VtkComposition, VtkWriteMolarities, false);
SET_BOOL_PROP(VtkComposition, VtkWriteFugacities, false);
SET_BOOL_PROP(VtkComposition, VtkWriteFugacityCoeffs, false);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the fluid composition
 *
 * This module deals with the following quantities:
 * - Mole fraction of a component in a fluid phase
 * - Mass fraction of a component in a fluid phase
 * - Molarity (i.e. molar concentration) of a component in a fluid phase
 * - Fugacity of all components
 * - FugacityCoefficient of all components in all phases
 */
template <class TypeTag>
class VtkCompositionModule : public BaseOutputModule<TypeTag>
{
    typedef BaseOutputModule<TypeTag> ParentType;

    typedef GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    static const int vtkFormat = GET_PROP_VALUE(TypeTag, VtkOutputFormat);
    typedef Ewoms::VtkMultiWriter<GridView, vtkFormat> VtkMultiWriter;

    typedef typename ParentType::ComponentBuffer ComponentBuffer;
    typedef typename ParentType::PhaseComponentBuffer PhaseComponentBuffer;

public:
    VtkCompositionModule(const Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteMassFractions,
                             "Include mass fractions in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteMoleFractions,
                             "Include mole fractions in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteTotalMassFractions,
                             "Include total mass fractions in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteTotalMoleFractions,
                             "Include total mole fractions in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteMolarities,
                             "Include component molarities in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFugacities,
                             "Include component fugacities in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFugacityCoeffs,
                             "Include component fugacity coefficients in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (moleFracOutput_())
            this->resizePhaseComponentBuffer_(moleFrac_);
        if (massFracOutput_())
            this->resizePhaseComponentBuffer_(massFrac_);
        if (totalMassFracOutput_())
            this->resizeComponentBuffer_(totalMassFrac_);
        if (totalMoleFracOutput_())
            this->resizeComponentBuffer_(totalMoleFrac_);
        if (molarityOutput_())
            this->resizePhaseComponentBuffer_(molarity_);

        if (fugacityOutput_())
            this->resizeComponentBuffer_(fugacity_);
        if (fugacityCoeffOutput_())
            this->resizePhaseComponentBuffer_(fugacityCoeff_);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant
     *        for an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        typedef Ewoms::MathToolbox<Evaluation> Toolbox;

        if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
            return;

        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(i, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    if (moleFracOutput_())
                        moleFrac_[phaseIdx][compIdx][I] = Toolbox::value(fs.moleFraction(phaseIdx, compIdx));
                    if (massFracOutput_())
                        massFrac_[phaseIdx][compIdx][I] = Toolbox::value(fs.massFraction(phaseIdx, compIdx));
                    if (molarityOutput_())
                        molarity_[phaseIdx][compIdx][I] = Toolbox::value(fs.molarity(phaseIdx, compIdx));

                    if (fugacityCoeffOutput_())
                        fugacityCoeff_[phaseIdx][compIdx][I] =
                            Toolbox::value(fs.fugacityCoefficient(phaseIdx, compIdx));
                }
            }

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                if (totalMassFracOutput_()) {
                    Scalar compMass = 0;
                    Scalar totalMass = 0;
                    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                        totalMass += Toolbox::value(fs.density(phaseIdx)) * Toolbox::value(fs.saturation(phaseIdx));
                        compMass +=
                            Toolbox::value(fs.density(phaseIdx))
                            *Toolbox::value(fs.saturation(phaseIdx))
                            *Toolbox::value(fs.massFraction(phaseIdx, compIdx));
                    }
                    totalMassFrac_[compIdx][I] = compMass / totalMass;
                }
                if (totalMoleFracOutput_()) {
                    Scalar compMoles = 0;
                    Scalar totalMoles = 0;
                    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                        totalMoles +=
                            Toolbox::value(fs.molarDensity(phaseIdx))
                            *Toolbox::value(fs.saturation(phaseIdx));
                        compMoles +=
                            Toolbox::value(fs.molarDensity(phaseIdx))
                            *Toolbox::value(fs.saturation(phaseIdx))
                            *Toolbox::value(fs.moleFraction(phaseIdx, compIdx));
                    }
                    totalMoleFrac_[compIdx][I] = compMoles / totalMoles;
                }
                if (fugacityOutput_())
                    fugacity_[compIdx][I] = Toolbox::value(intQuants.fluidState().fugacity(/*phaseIdx=*/0, compIdx));
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter& baseWriter)
    {
        VtkMultiWriter *vtkWriter = dynamic_cast<VtkMultiWriter*>(&baseWriter);
        if (!vtkWriter) {
            return;
        }

        if (moleFracOutput_())
            this->commitPhaseComponentBuffer_(baseWriter, "moleFrac_%s^%s", moleFrac_);
        if (massFracOutput_())
            this->commitPhaseComponentBuffer_(baseWriter, "massFrac_%s^%s", massFrac_);
        if (molarityOutput_())
            this->commitPhaseComponentBuffer_(baseWriter, "molarity_%s^%s", molarity_);
        if (totalMassFracOutput_())
            this->commitComponentBuffer_(baseWriter, "totalMassFrac^%s", totalMassFrac_);
        if (totalMoleFracOutput_())
            this->commitComponentBuffer_(baseWriter, "totalMoleFrac^%s", totalMoleFrac_);

        if (fugacityOutput_())
            this->commitComponentBuffer_(baseWriter, "fugacity^%s", fugacity_);
        if (fugacityCoeffOutput_())
            this->commitPhaseComponentBuffer_(baseWriter, "fugacityCoeff_%s^%s", fugacityCoeff_);
    }

private:
    static bool massFracOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteMassFractions);
        return val;
    }

    static bool moleFracOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteMoleFractions);
        return val;
    }

    static bool totalMassFracOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteTotalMassFractions);
        return val;
    }

    static bool totalMoleFracOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteTotalMoleFractions);
        return val;
    }

    static bool molarityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteMolarities);
        return val;
    }

    static bool fugacityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFugacities);
        return val;
    }

    static bool fugacityCoeffOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFugacityCoeffs);
        return val;
    }

    PhaseComponentBuffer moleFrac_;
    PhaseComponentBuffer massFrac_;
    PhaseComponentBuffer molarity_;
    ComponentBuffer totalMassFrac_;
    ComponentBuffer totalMoleFrac_;

    ComponentBuffer fugacity_;
    PhaseComponentBuffer fugacityCoeff_;
};

} // namespace Ewoms

#endif
