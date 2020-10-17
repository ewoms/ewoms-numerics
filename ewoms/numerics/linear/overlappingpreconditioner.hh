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
 * \copydoc Ewoms::Linear::OverlappingPreconditioner
 */
#ifndef EWOMS_OVERLAPPING_PRECONDITIONER_HH
#define EWOMS_OVERLAPPING_PRECONDITIONER_HH

#include "overlappingscalarproduct.hh"

#include <ewoms/common/exceptions.hh>

#include <dune/istl/preconditioners.hh>

#include <dune/common/version.hh>

namespace Ewoms {
namespace Linear {

/*!
 * \brief An overlap aware preconditioner for any ISTL linear solver.
 */
template <class SeqPreCond, class Overlap>
class OverlappingPreconditioner
    : public Dune::Preconditioner<typename SeqPreCond::domain_type,
                                  typename SeqPreCond::range_type>
{
public:
    typedef typename SeqPreCond::domain_type domain_type;
    typedef typename SeqPreCond::range_type range_type;

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2,6)
    //! the kind of computations supported by the operator. Either overlapping or non-overlapping
    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::overlapping; }
#else
    // redefine the category
    enum { category = Dune::SolverCategory::overlapping };
#endif

    OverlappingPreconditioner(SeqPreCond& seqPreCond, const Overlap& overlap)
        : seqPreCond_(seqPreCond), overlap_(&overlap)
    {}

    void pre(domain_type& x, range_type& y) override
    {
#if HAVE_MPI
        short success;
        try
        {
            seqPreCond_.pre(x, y);
            short localSuccess = 1;
            MPI_Allreduce(&localSuccess,   // source buffer
                          &success,        // destination buffer
                          1,               // number of objects in buffers
                          MPI_SHORT,       // data type
                          MPI_MIN,         // operation
                          MPI_COMM_WORLD); // communicator
        }
        catch (...)
        {
            short localSuccess = 0;
            MPI_Allreduce(&localSuccess,   // source buffer
                          &success,        // destination buffer
                          1,               // number of objects in buffers
                          MPI_SHORT,       // data type
                          MPI_MIN,         // operation
                          MPI_COMM_WORLD); // communicator
        }

        if (success) {
            x.sync();
        }
        else
            throw Ewoms::NumericalIssue("Preconditioner threw an exception in pre() method on some process.");
#else
        seqPreCond_.pre(x, y);
#endif

        // communicate the results on the overlap
        x.sync();
        y.sync();
    }

    void apply(domain_type& x, const range_type& d) override
    {
#if HAVE_MPI
        if (overlap_->peerSet().size() > 0) {
            // make sure that all processes react the same if the
            // sequential preconditioner on one process throws an
            // exception
            short success;
            try
            {
                // execute the sequential preconditioner
                seqPreCond_.apply(x, d);
                short localSuccess = 1;
                MPI_Allreduce(&localSuccess,   // source buffer
                              &success,        // destination buffer
                              1,               // number of objects in buffers
                              MPI_SHORT,       // data type
                              MPI_MIN,         // operation
                              MPI_COMM_WORLD); // communicator
            }
            catch (...)
            {
                short localSuccess = 0;
                MPI_Allreduce(&localSuccess,   // source buffer
                              &success,        // destination buffer
                              1,               // number of objects in buffers
                              MPI_SHORT,       // data type
                              MPI_MIN,         // operation
                              MPI_COMM_WORLD); // communicator
            }

            if (success) {
                x.sync();
            }
            else
                throw Ewoms::NumericalIssue("Preconditioner threw an exception on some process.");
        }
        else
#endif // HAVE_MPI
            seqPreCond_.apply(x, d);
    }

    void post(domain_type& x) override
    {
#if HAVE_MPI
        short success;
        try
        {
            seqPreCond_.post(x);
            short localSuccess = 1;
            MPI_Allreduce(&localSuccess,   // source buffer
                          &success,        // destination buffer
                          1,               // number of objects in buffers
                          MPI_SHORT,       // data type
                          MPI_MIN,         // operation
                          MPI_COMM_WORLD); // communicator
        }
        catch (...)
        {
            short localSuccess = 0;
            MPI_Allreduce(&localSuccess,   // source buffer
                          &success,        // destination buffer
                          1,               // number of objects in buffers
                          MPI_SHORT,       // data type
                          MPI_MIN,         // operation
                          MPI_COMM_WORLD); // communicator
        }

        if (success) {
            x.sync();
        }
        else
            throw Ewoms::NumericalIssue("Preconditioner threw an exception in post() method on "
                                        "some process.");
#else
        seqPreCond_.post(x);
#endif
    }

private:
    SeqPreCond& seqPreCond_;
    const Overlap *overlap_;
};

} // namespace Linear
} // namespace Ewoms

#endif
