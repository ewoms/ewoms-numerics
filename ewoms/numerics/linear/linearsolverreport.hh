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
 * \copydoc Ewoms::Linear::SolverReport
 */
#ifndef EWOMS_LINEAR_SOLVER_REPORT_HH
#define EWOMS_LINEAR_SOLVER_REPORT_HH

#include "convergencecriterion.hh"

#include <ewoms/common/timer.hh>
#include <ewoms/common/timerguard.hh>

namespace Ewoms {
namespace Linear {

/*!
 * \brief Collects summary information about the execution of the linear solver.
 */
class SolverReport
{
public:
    SolverReport()
    { reset(); }

    void reset()
    {
        timer_.halt();
        iterations_ = 0;
        converged_ = 0;
    }

    const Ewoms::Timer& timer() const
    { return timer_; }

    Ewoms::Timer& timer()
    { return timer_; }

    unsigned iterations() const
    { return iterations_; }

    void increment()
    { ++iterations_; }

    SolverReport& operator++()
    { ++iterations_; return *this; }

    bool converged() const
    { return converged_; }

    void setConverged(bool value)
    { converged_ = value; }

private:
    Ewoms::Timer timer_;
    unsigned iterations_;
    bool converged_;
};

}} // end namespace Linear, Ewoms

#endif
