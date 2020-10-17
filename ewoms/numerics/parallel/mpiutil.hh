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
 * \copydoc Ewoms::MpiBuffer
 */
#ifndef EWOMS_MATERIAL_MPIUTIL_HH
#define EWOMS_MATERIAL_MPIUTIL_HH

#include <dune/common/parallel/mpitraits.hh>

#include <cassert>
#include <numeric>
#include <string>
#include <vector>

#if HAVE_MPI

#include <mpi.h>

namespace Ewoms
{

    template <typename T>
    int packSize()
    {
        int packSize;
        MPI_Pack_size(1, Dune::MPITraits<T>::getType(), MPI_COMM_WORLD, &packSize);
        return packSize;
    }

    // --------  Packer --------
    template <typename T>
    struct Packer
    {
        static int size(const T&)
        {
            return packSize<T>();
        }

        static void pack(const T& content, std::vector<char>& buf, int& offset)
        {
            MPI_Pack(&content, 1, Dune::MPITraits<T>::getType(), buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
        }

        static T unpack(const std::vector<char>& recvBuffer, int& offset)
        {
            T content;
            auto* data = const_cast<char*>(recvBuffer.data());
            MPI_Unpack(data, recvBuffer.size(), &offset, &content, 1, Dune::MPITraits<T>::getType(), MPI_COMM_WORLD);
            return content;
        }
    };

    // --------  Packer, string specialization --------
    template <>
    struct Packer<std::string>
    {
        static int size(const std::string& content)
        {
            return packSize<unsigned int>() + content.size()*packSize<char>();
        }

        static void pack(const std::string& content, std::vector<char>& buf, int& offset)
        {
            unsigned int size = content.size();
            Packer<unsigned int>::pack(size, buf, offset);
            if (size > 0) {
                MPI_Pack(const_cast<char*>(content.c_str()), size, MPI_CHAR, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
            }
        }

        static std::string unpack(const std::vector<char>& recvBuffer, int& offset)
        {
            unsigned int size = Packer<unsigned int>::unpack(recvBuffer, offset);
            std::string text;
            if (size > 0) {
                auto* data = const_cast<char*>(recvBuffer.data());
                std::vector<char> chars(size);
                MPI_Unpack(data, recvBuffer.size(), &offset, chars.data(), size, MPI_CHAR, MPI_COMM_WORLD);
                text = std::string(chars.data(), size);
            }
            return text;
        }
    };

    // --------  Packer, vector partial specialization --------
    template <typename T>
    struct Packer<std::vector<T>>
    {
        static int size(const std::string& content)
        {
            int sz = 0;
            sz += packSize<unsigned int>();
            for (const T& elem : content) {
                sz += Packer<T>::size(elem);
            }
            return sz;
        }

        static void pack(const std::vector<T>& content, std::vector<char>& buf, int& offset)
        {
            unsigned int size = content.size();
            Packer<unsigned int>::pack(size, buf, offset);
            for (const T& elem : content) {
                Packer<T>::pack(elem);
            }
        }

        static std::vector<T> unpack(const std::vector<char>& recvBuffer, int& offset)
        {
            unsigned int size = Packer<T>::unpack(recvBuffer, offset);
            std::vector<T> content;
            content.reserve(size);
            for (unsigned int i = 0; i < size; ++i) {
                content.push_back(Packer<T>::unpack(recvBuffer, offset));
            }
            return content;
        }
    };

} // anonymous namespace

namespace Ewoms
{

    /// From each rank, gather its string (if not empty) into a vector.
    inline std::vector<std::string> gatherStrings(const std::string& localString)
    {
        using StringPacker = Ewoms::Packer<std::string>;

        // Pack local messages.
        const int messageSize = StringPacker::size(localString);
        std::vector<char> buffer(messageSize);
        int offset = 0;
        StringPacker::pack(localString, buffer, offset);
        assert(offset == messageSize);

        // Get message sizes and create offset/displacement array for gathering.
        int numProcesses = -1;
        MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
        std::vector<int> messageSizes(numProcesses);
        MPI_Allgather(&messageSize, 1, MPI_INT, messageSizes.data(), 1, MPI_INT, MPI_COMM_WORLD);
        std::vector<int> displ(numProcesses + 1, 0);
        std::partial_sum(messageSizes.begin(), messageSizes.end(), displ.begin() + 1);

        // Gather.
        std::vector<char> recvBuffer(displ.back());
        MPI_Allgatherv(buffer.data(), buffer.size(), MPI_PACKED,
                       const_cast<char*>(recvBuffer.data()), messageSizes.data(),
                       displ.data(), MPI_PACKED,
                       MPI_COMM_WORLD);

        // Unpack and return.
        std::vector<std::string> ret;
        for (int process = 0; process < numProcesses; ++process) {
            offset = displ[process];
            std::string s = StringPacker::unpack(recvBuffer, offset);
            if (!s.empty()) {
                ret.push_back(s);
            }
            assert(offset == displ[process + 1]);
        }
        return ret;
    }

} // namespace Ewoms

#else // HAVE_MPI

namespace Ewoms
{
    inline std::vector<std::string> gatherStrings(const std::string& localString)
    {
        if (localString.empty()) {
            return {};
        } else {
            return { localString };
        }
    }
} // namespace Ewoms

#endif // HAVE_MPI

#endif // EWOMS_MATERIAL_MPIUTIL_HH

