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
 * \copydoc Ewoms::Linear::SuperLUBackend
 */
#ifndef EWOMS_SUPER_LU_BACKEND_HH
#define EWOMS_SUPER_LU_BACKEND_HH

#if HAVE_SUPERLU

#include <ewoms/numerics/linear/istlsparsematrixbackend.hh>
#include <ewoms/common/parametersystem.hh>

#include <ewoms/common/unused.hh>

#include <dune/istl/superlu.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>

BEGIN_PROPERTIES

// forward declaration of the required property tags
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(NumEq);
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(SparseMatrixAdapter);
NEW_PROP_TAG(GlobalEqVector);
NEW_PROP_TAG(LinearSolverVerbosity);
NEW_PROP_TAG(LinearSolverBackend);
NEW_TYPE_TAG(SuperLULinearSolver);

END_PROPERTIES

namespace Ewoms {
namespace Linear {
template <class Scalar, class TypeTag, class Matrix, class Vector>
class SuperLUSolve_;

/*!
 * \ingroup Linear
 * \brief A linear solver backend for the SuperLU sparse matrix library.
 */
template <class TypeTag>
class SuperLUBackend
{
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using SparseMatrixAdapter = GET_PROP_TYPE(TypeTag, SparseMatrixAdapter);
    using Matrix = typename SparseMatrixAdapter::block_type;

    static_assert(std::is_same<SparseMatrixAdapter, IstlSparseMatrixAdapter<MatrixBlock>::value,
                  "The SuperLU linear solver backend requires the IstlSparseMatrixAdapter");

public:
    SuperLUBackend(Simulator& simulator EWOMS_UNUSED)
    {}

    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverVerbosity,
                             "The verbosity level of the linear solver");
    }

    /*!
     * \brief Causes the solve() method to discared the structure of the linear system of
     *        equations the next time it is called.
     *
     * Since the SuperLU backend does not create any internal matrices, this is a no-op.
     */
    void eraseMatrix()
    { }

    void prepare(const SparseMatrixAdapter& M, const Vector& b)
    { }

    void setResidual(const Vector& b)
    { b_ = &b; }

    void getResidual(Vector& b) const
    { b = *b_; }

    void setMatrix(const SparseMatrixAdapter& M)
    { M_ = &M; }

    bool solve(Vector& x)
    { return SuperLUSolve_<Scalar, TypeTag, Matrix, Vector>::solve_(*M_, x, *b_); }

private:
    const Matrix* M_;
    Vector* b_;
};

template <class Scalar, class TypeTag, class Matrix, class Vector>
class SuperLUSolve_
{
public:
    static bool solve_(const Matrix& A, Vector& x, const Vector& b)
    {
        Vector bTmp(b);

        int verbosity = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);
        Dune::InverseOperatorResult result;
        Dune::SuperLU<Matrix> solver(A, verbosity > 0);
        solver.apply(x, bTmp, result);

        if (result.converged) {
            // make sure that the result only contains finite values.
            Scalar tmp = 0;
            for (unsigned i = 0; i < x.size(); ++i) {
                const auto& xi = x[i];
                for (unsigned j = 0; j < Vector::block_type::dimension; ++j)
                    tmp += xi[j];
            }
            result.converged = std::isfinite(tmp);
        }

        return result.converged;
    }
};

// the following is required to make the SuperLU adapter of dune-istl happy with
// quadruple precision math on Dune 2.4. this is because the most which SuperLU can
// handle is double precision (i.e., the linear systems of equations are always solved
// with at most double precision if chosing SuperLU as the linear solver...)
#if HAVE_QUAD
template <class TypeTag, class Matrix, class Vector>
class SuperLUSolve_<__float128, TypeTag, Matrix, Vector>
{
public:
    static bool solve_(const Matrix& A,
                       Vector& x,
                       const Vector& b)
    {
        static const int numEq = GET_PROP_VALUE(TypeTag, NumEq);
        using DoubleEqVector = Dune::FieldVector<double, numEq>;
        using DoubleEqMatrix = Dune::FieldMatrix<double, numEq, numEq>;
        using DoubleVector = Dune::BlockVector<DoubleEqVector>;
        using DoubleMatrix = Dune::BCRSMatrix<DoubleEqMatrix>;

        // copy the inputs into the double precision data structures
        DoubleVector bDouble(b);
        DoubleVector xDouble(x);
        DoubleMatrix ADouble(A);

        bool res =
            SuperLUSolve_<double, TypeTag, Matrix, Vector>::solve_(ADouble,
                                                                   xDouble,
                                                                   bDouble);

        // copy the result back into the quadruple precision vector.
        x = xDouble;

        return res;
    }
};
#endif

} // namespace Linear
} // namespace Ewoms

BEGIN_PROPERTIES

SET_INT_PROP(SuperLULinearSolver, LinearSolverVerbosity, 0);
SET_TYPE_PROP(SuperLULinearSolver, LinearSolverBackend,
              Ewoms::Linear::SuperLUBackend<TypeTag>);

END_PROPERTIES

#endif // HAVE_SUPERLU

#endif
