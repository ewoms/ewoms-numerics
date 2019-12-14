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
 * \copydoc Ewoms::DiscreteFractureModel
 */
#ifndef EWOMS_DISCRETE_FRACTURE_MODEL_HH
#define EWOMS_DISCRETE_FRACTURE_MODEL_HH

#include <ewoms/common/densead/math.hh>

#include "discretefractureproperties.hh"
#include "discretefractureprimaryvariables.hh"
#include "discretefractureintensivequantities.hh"
#include "discretefractureextensivequantities.hh"
#include "discretefracturelocalresidual.hh"
#include "discretefractureproblem.hh"

#include <ewoms/numerics/models/immiscible/immisciblemodel.hh>
#include <ewoms/numerics/io/vtkdiscretefracturemodule.hh>

#include <ewoms/common/exceptions.hh>

#include <string>

namespace Ewoms {
template <class TypeTag>
class DiscreteFractureModel;
}

BEGIN_PROPERTIES

//! The generic type tag for problems using the immiscible multi-phase model
NEW_TYPE_TAG(DiscreteFractureModel, INHERITS_FROM(ImmiscibleTwoPhaseModel, VtkDiscreteFracture));

//! The class for the model
SET_TYPE_PROP(DiscreteFractureModel, Model, Ewoms::DiscreteFractureModel<TypeTag>);

//! The class for the model
SET_TYPE_PROP(DiscreteFractureModel, BaseProblem, Ewoms::DiscreteFractureProblem<TypeTag>);

//! Use the immiscible multi-phase local jacobian operator for the immiscible multi-phase model
SET_TYPE_PROP(DiscreteFractureModel, LocalResidual, Ewoms::DiscreteFractureLocalResidual<TypeTag>);

// The type of the base base class for actual problems.
// TODO!?
//SET_TYPE_PROP(DiscreteFractureModel BaseProblem, DiscreteFractureBaseProblem<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(DiscreteFractureModel, PrimaryVariables,
              Ewoms::DiscreteFracturePrimaryVariables<TypeTag>);

//! the IntensiveQuantities property
SET_TYPE_PROP(DiscreteFractureModel, IntensiveQuantities,
              Ewoms::DiscreteFractureIntensiveQuantities<TypeTag>);

//! the ExtensiveQuantities property
SET_TYPE_PROP(DiscreteFractureModel, ExtensiveQuantities,
              Ewoms::DiscreteFractureExtensiveQuantities<TypeTag>);

//! For the discrete fracture model, we need to use two-point flux approximation or it
//! will converge very poorly
SET_BOOL_PROP(DiscreteFractureModel, UseTwoPointGradients, true);

// The intensive quantity cache cannot be used by the discrete fracture model, because
// the intensive quantities of a control degree of freedom are not identical to the
// intensive quantities of the other intensive quantities of the same of the same degree
// of freedom. This is because the fracture properties (volume, permeability, etc) are
// specific for each...
SET_BOOL_PROP(DiscreteFractureModel, EnableIntensiveQuantityCache, false);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup DiscreteFractureModel
 * \brief A fully-implicit multi-phase flow model which assumes
 *        immiscibility of the phases and is able to include fractures
 *        in the domain.
 *
 * This model implements multi-phase flow of \f$M > 0\f$ immiscible
 * fluids \f$\alpha\f$. It also can consider edges of the
 * computational grid as fractures i.e. as a porous medium with
 * different higher permeability than the rest of the domain.
 *
 * \todo So far, the discrete fracture model only works for 2D grids
 *       and without energy. Also only the Darcy velocity model is
 *       supported for the fractures.
 *
 * \sa ImmiscibleModel
 */
template <class TypeTag>
class DiscreteFractureModel : public ImmiscibleModel<TypeTag>
{
    typedef ImmiscibleModel<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

public:
    DiscreteFractureModel(Simulator& simulator)
        : ParentType(simulator)
    {
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableIntensiveQuantityCache)) {
            throw std::runtime_error("The discrete fracture model does not work in conjunction "
                                     "with intensive quantities caching");
        }
    }

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        // register runtime parameters of the VTK output modules
        Ewoms::VtkDiscreteFractureModule<TypeTag>::registerParameters();
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    static std::string name()
    { return "discretefracture"; }

    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        this->addOutputModule(new Ewoms::VtkDiscreteFractureModule<TypeTag>(this->simulator_));
    }
};
} // namespace Ewoms

#endif
