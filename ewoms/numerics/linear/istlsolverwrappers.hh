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
 * \brief Provides wrapper classes for the iterative linear solvers available in
 *        dune-istl.
 *
 * In conjunction with a suitable solver backend, solver wrappers work by specifying the
 * "SolverWrapper" property:
 * \code
 * SET_TYPE_PROP(YourTypeTag, LinearSolverWrapper,
 *               Ewoms::Linear::SolverWrapper$SOLVER<TypeTag>);
 * \endcode
 *
 * The possible choices for '\c $SOLVER' are:
 * - \c Richardson: A fixpoint solver using the Richardson iteration
 * - \c SteepestDescent: The steepest descent solver
 * - \c ConjugatedGradients: A conjugated gradients solver
 * - \c BiCGStab: A stabilized bi-conjugated gradients solver
 * - \c MinRes: A solver based on the  minimized residual algorithm
 * - \c RestartedGMRes: A restarted GMRES solver
 */
#ifndef EWOMS_ISTL_SOLVER_WRAPPERS_HH
#define EWOMS_ISTL_SOLVER_WRAPPERS_HH

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <dune/istl/solvers.hh>

BEGIN_PROPERTIES

NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(SparseMatrixAdapter);
NEW_PROP_TAG(OverlappingMatrix);
NEW_PROP_TAG(OverlappingVector);
NEW_PROP_TAG(GMResRestart);
NEW_PROP_TAG(LinearSolverTolerance);
NEW_PROP_TAG(LinearSolverMaxIterations);
NEW_PROP_TAG(LinearSolverVerbosity);

END_PROPERTIES

namespace Ewoms {

namespace Linear {

/*!
 * \brief Macro to create a wrapper around an ISTL solver
 */
#define EWOMS_WRAP_ISTL_SOLVER(SOLVER_NAME, ISTL_SOLVER_NAME)                      \
    template <class TypeTag>                                                       \
    class SolverWrapper##SOLVER_NAME                                               \
    {                                                                              \
        using Scalar = GET_PROP_TYPE(TypeTag, Scalar);                    \
        using OverlappingMatrix = GET_PROP_TYPE(TypeTag, OverlappingMatrix);       \
        using OverlappingVector = GET_PROP_TYPE(TypeTag, OverlappingVector);       \
                                                                                   \
    public:                                                                        \
        using RawSolver = ISTL_SOLVER_NAME<OverlappingVector>;                     \
                                                                                   \
        SolverWrapper##SOLVER_NAME()                                               \
        {}                                                                         \
                                                                                   \
        static void registerParameters()                                           \
        {}                                                                         \
                                                                                   \
        template <class LinearOperator, class ScalarProduct, class Preconditioner> \
        std::shared_ptr<RawSolver> get(LinearOperator& parOperator,                \
                                       ScalarProduct& parScalarProduct,            \
                                       Preconditioner& parPreCond)                 \
        {                                                                          \
            Scalar tolerance = EWOMS_GET_PARAM(TypeTag, Scalar,                    \
                                               LinearSolverTolerance);             \
            int maxIter = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIterations);\
                                                                                   \
            int verbosity = 0;                                                     \
            if (parOperator.overlap().myRank() == 0)                               \
                verbosity = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);  \
            solver_ = std::make_shared<RawSolver>(parOperator, parScalarProduct,   \
                                                  parPreCond, tolerance, maxIter,  \
                                                  verbosity);                      \
                                                                                   \
            return solver_;                                                        \
        }                                                                          \
                                                                                   \
        void cleanup()                                                             \
        { solver_.reset(); }                                                       \
                                                                                   \
    private:                                                                       \
        std::shared_ptr<RawSolver> solver_;                                        \
    };

EWOMS_WRAP_ISTL_SOLVER(Richardson, Dune::LoopSolver)
EWOMS_WRAP_ISTL_SOLVER(SteepestDescent, Dune::GradientSolver)
EWOMS_WRAP_ISTL_SOLVER(ConjugatedGradients, Dune::CGSolver)
EWOMS_WRAP_ISTL_SOLVER(BiCGStab, Dune::BiCGSTABSolver)
EWOMS_WRAP_ISTL_SOLVER(MinRes, Dune::MINRESSolver)

/*!
 * \brief Solver wrapper for the restarted GMRES solver of dune-istl.
 *
 * dune-istl uses a slightly different API for this solver than for the others...
 */
template <class TypeTag>
class SolverWrapperRestartedGMRes
{
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using OverlappingMatrix = GET_PROP_TYPE(TypeTag, OverlappingMatrix);
    using OverlappingVector = GET_PROP_TYPE(TypeTag, OverlappingVector);

public:
    using RawSolver = Dune::RestartedGMResSolver<OverlappingVector>;

    SolverWrapperRestartedGMRes()
    {}

    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, int, GMResRestart,
                             "Number of iterations after which the GMRES linear solver is restarted");
    }

    template <class LinearOperator, class ScalarProduct, class Preconditioner>
    std::shared_ptr<RawSolver> get(LinearOperator& parOperator,
                                   ScalarProduct& parScalarProduct,
                                   Preconditioner& parPreCond)
    {
        Scalar tolerance = EWOMS_GET_PARAM(TypeTag, Scalar, LinearSolverTolerance);
        int maxIter = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIterations);

        int verbosity = 0;
        if (parOperator.overlap().myRank() == 0)
            verbosity = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);
        int restartAfter = EWOMS_GET_PARAM(TypeTag, int, GMResRestart);
        solver_ = std::make_shared<RawSolver>(parOperator,
                                              parScalarProduct,
                                              parPreCond,
                                              tolerance,
                                              restartAfter,
                                              maxIter,
                                              verbosity);

        return solver_;
    }

    void cleanup()
    { solver_.reset(); }

private:
    std::shared_ptr<RawSolver> solver_;
};

#undef EWOMS_WRAP_ISTL_SOLVER

}} // namespace Linear, Ewoms

#endif
