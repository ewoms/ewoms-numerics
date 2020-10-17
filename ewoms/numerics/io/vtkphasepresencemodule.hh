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
 * \copydoc Ewoms::VtkPhasePresenceModule
 */
#ifndef EWOMS_VTK_PHASE_PRESENCE_MODULE_HH
#define EWOMS_VTK_PHASE_PRESENCE_MODULE_HH

#include "vtkmultiwriter.hh"
#include "baseoutputmodule.hh"

#include <ewoms/common/parametersystem.hh>
#include <ewoms/common/propertysystem.hh>

BEGIN_PROPERTIES

// create new type tag for the VTK primary variables output
NEW_TYPE_TAG(VtkPhasePresence);

// create the property tags needed for the primary variables module
NEW_PROP_TAG(VtkWritePhasePresence);
NEW_PROP_TAG(VtkOutputFormat);
NEW_PROP_TAG(EnableVtkOutput);

SET_BOOL_PROP(VtkPhasePresence, VtkWritePhasePresence, false);

END_PROPERTIES

namespace Ewoms {
/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the fluid composition
 */
template<class TypeTag>
class VtkPhasePresenceModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);
    using GridView = GET_PROP_TYPE(TypeTag, GridView);

    static const int vtkFormat = GET_PROP_VALUE(TypeTag, VtkOutputFormat);
    using VtkMultiWriter = Ewoms::VtkMultiWriter<GridView, vtkFormat>;

    using ScalarBuffer = typename ParentType::ScalarBuffer;

public:
    VtkPhasePresenceModule(const Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePhasePresence,
                             "Include the phase presence pseudo primary "
                             "variable in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (phasePresenceOutput_()) this->resizeScalarBuffer_(phasePresence_);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quanties relevant
     *        for an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
            return;

        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            // calculate the phase presence
            int phasePresence = elemCtx.primaryVars(i, /*timeIdx=*/0).phasePresence();
            unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);

            if (phasePresenceOutput_())
                phasePresence_[I] = phasePresence;
        }
    }

    /*!
     * \brief Add all buffers to the output writer.
     */
    void commitBuffers(BaseOutputWriter& baseWriter)
    {
        VtkMultiWriter *vtkWriter = dynamic_cast<VtkMultiWriter*>(&baseWriter);
        if (!vtkWriter) {
            return;
        }

        if (phasePresenceOutput_())
            this->commitScalarBuffer_(baseWriter, "phase presence", phasePresence_);
    }

private:
    static bool phasePresenceOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWritePhasePresence);
        return val;
    }

    ScalarBuffer phasePresence_;
};

} // namespace Ewoms

#endif
