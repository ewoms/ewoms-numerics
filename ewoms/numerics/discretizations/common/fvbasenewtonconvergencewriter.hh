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
 * \copydoc Ewoms::FvBaseNewtonConvergenceWriter
 */
#ifndef EWOMS_FV_BASE_NEWTON_CONVERGENCE_WRITER_HH
#define EWOMS_FV_BASE_NEWTON_CONVERGENCE_WRITER_HH

#include <ewoms/numerics/io/vtkmultiwriter.hh>
#include <ewoms/common/propertysystem.hh>

#include <iostream>

//! \cond SKIP_THIS
BEGIN_PROPERTIES

// forward declaration of the required property tags
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(NewtonMethod);
NEW_PROP_TAG(SolutionVector);
NEW_PROP_TAG(GlobalEqVector);
NEW_PROP_TAG(VtkOutputFormat);

END_PROPERTIES
//! \endcond

namespace Ewoms {
/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief Writes the intermediate solutions during the Newton scheme
 *        for models using a finite volume discretization
 */
template <class TypeTag>
class FvBaseNewtonConvergenceWriter
{
    using GridView = GET_PROP_TYPE(TypeTag, GridView);

    using SolutionVector = GET_PROP_TYPE(TypeTag, SolutionVector);
    using GlobalEqVector = GET_PROP_TYPE(TypeTag, GlobalEqVector);
    using NewtonMethod = GET_PROP_TYPE(TypeTag, NewtonMethod);

    static const int vtkFormat = GET_PROP_VALUE(TypeTag, VtkOutputFormat);
    using VtkMultiWriter = Ewoms::VtkMultiWriter<GridView, vtkFormat>;

public:
    FvBaseNewtonConvergenceWriter(NewtonMethod& nm)
        : newtonMethod_(nm)
    {
        timeStepIdx_ = 0;
        iteration_ = 0;
        vtkMultiWriter_ = 0;
    }

    ~FvBaseNewtonConvergenceWriter()
    { delete vtkMultiWriter_; }

    /*!
     * \brief Called by the Newton method before the actual algorithm
     *        is started for any given timestep.
     */
    void beginTimeStep()
    {
        ++timeStepIdx_;
        iteration_ = 0;
    }

    /*!
     * \brief Called by the Newton method before an iteration of the
     *        Newton algorithm is started.
     */
    void beginIteration()
    {
        ++ iteration_;
        if (!vtkMultiWriter_)
            vtkMultiWriter_ =
                new VtkMultiWriter(/*async=*/false,
                                   newtonMethod_.problem().gridView(),
                                   newtonMethod_.problem().outputDir(),
                                   "convergence");
        vtkMultiWriter_->beginWrite(timeStepIdx_ + iteration_ / 100.0);
    }

    /*!
     * \brief Write the Newton update to disk.
     *
     * Called after the linear solution is found for an iteration.
     *
     * \param uLastIter The solution vector of the previous iteration.
     * \param deltaU The negative difference between the solution
     *        vectors of the previous and the current iteration.
     */
    void writeFields(const SolutionVector& uLastIter,
                     const GlobalEqVector& deltaU)
    {
        try {
            newtonMethod_.problem().model().addConvergenceVtkFields(*vtkMultiWriter_,
                                                                    uLastIter,
                                                                    deltaU);
        }
        catch (...) {
            std::cout << "Oops: exception thrown on rank "
                      << newtonMethod_.problem().gridView().comm().rank()
                      << " while writing the convergence\n"  << std::flush;
        }
    }

    /*!
     * \brief Called by the Newton method after an iteration of the
     *        Newton algorithm has been completed.
     */
    void endIteration()
    { vtkMultiWriter_->endWrite(); }

    /*!
     * \brief Called by the Newton method after Newton algorithm
     *        has been completed for any given timestep.
     *
     * This method is called regardless of whether the Newton method
     * converged or not.
     */
    void endTimeStep()
    { iteration_ = 0; }

private:
    int timeStepIdx_;
    int iteration_;
    VtkMultiWriter *vtkMultiWriter_;
    NewtonMethod& newtonMethod_;
};

} // namespace Ewoms

#endif
