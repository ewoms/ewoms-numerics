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
#include <config.h>

#include <ewoms/numerics/parallel/mpiutil.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <cassert>

#if HAVE_MPI
struct MPIError
{
    MPIError(std::string s, int e) : errorstring(std::move(s)), errorcode(e){}
    std::string errorstring;
    int errorcode;
};

void mpiErrorHandler(MPI_Comm*, int* errCode, ...)
{
    std::vector<char> errString(MPI_MAX_ERROR_STRING);
    int errLength;
    MPI_Error_string(*errCode, errString.data(), &errLength);
    std::string s(errString.data(), errLength);
    std::cerr << "An MPI Error ocurred:" << std::endl << s << std::endl;
    throw MPIError(s, *errCode);
}
#endif

bool noStrings(int, int)
{
    std::string empty;
    auto res = Ewoms::gatherStrings(empty);
    assert(res.empty());
    return true;
}

bool oddRankStrings(int size, int rank)
{
    std::string what = (rank % 2 == 1) ? "An error on rank " + std::to_string(rank) : std::string();
    auto res = Ewoms::gatherStrings(what);
    assert(int(res.size()) == size/2);
    for (int i = 0; i < size/2; ++i) {
        assert(res[i] == "An error on rank " + std::to_string(2*i + 1));
    }
    return true;
}

bool allRankStrings(int size, int rank)
{
    std::string what = "An error on rank " + std::to_string(rank);
    auto res = Ewoms::gatherStrings(what);
    assert(int(res.size()) == size);
    for (int i = 0; i < size; ++i) {
        assert(res[i] == "An error on rank " + std::to_string(i));
    }
    return true;
}

int testMain(int size, int rank)
{
    bool ok = noStrings(size, rank);
    ok = ok && oddRankStrings(size, rank);
    ok = ok && allRankStrings(size, rank);
    if (ok) {
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }
}

int main(int argc, char** argv)
{
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);
    int mpiSize = mpiHelper.size();
    int mpiRank = mpiHelper.rank();
#if HAVE_MPI
    // register a throwing error handler to allow for
    // debugging with "catch throw" in gdb
    MPI_Errhandler handler;
    MPI_Comm_create_errhandler(mpiErrorHandler, &handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
#endif
    return testMain(mpiSize, mpiRank);
}
