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
 *
 * \copydoc Ewoms::BaseOutputWriter
 */
#ifndef EWOMS_BASE_OUTPUT_WRITER_HH
#define EWOMS_BASE_OUTPUT_WRITER_HH

#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>

#include <vector>

namespace Ewoms {
/*!
 * \brief The base class for all output writers.
 *
 * The sole purpose of this class is to enable RTTI
 * (i.e. dynamic_cast) on writer objects.
 */
class BaseOutputWriter
{
public:
    using Scalar = double;
    using Vector = Dune::DynamicVector<double>;
    using Tensor = Dune::DynamicMatrix<double>;
    using ScalarBuffer = std::vector<Scalar>;
    using VectorBuffer = std::vector<Vector>;
    using TensorBuffer = std::vector<Tensor>;

    BaseOutputWriter()
    {}

    virtual ~BaseOutputWriter()
    {}

    /*!
     * \brief Called when ever a new time step or a new grid
     *        must be written.
     */
    virtual void beginWrite(double t) = 0;

    /*!
     * \brief Add a scalar vertex centered vector field to the output.
     */
    virtual void attachScalarVertexData(ScalarBuffer& buf, std::string name) = 0;

    /*!
     * \brief Add a scalar element centered quantity to the output.
     */
    virtual void attachScalarElementData(ScalarBuffer& buf, std::string name) = 0;

    /*!
     * \brief Add a vectorial vertex centered vector field to the output.
     */
    virtual void attachVectorVertexData(VectorBuffer& buf, std::string name) = 0;

    /*!
     * \brief Add a vectorial element centered quantity to the output.
     */
    virtual void attachVectorElementData(VectorBuffer& buf, std::string name) = 0;

    /*!
     * \brief Add a tensorial vertex centered tensor field to the output.
     */
    virtual void attachTensorVertexData(TensorBuffer& buf, std::string name) = 0;

    /*!
     * \brief Add a tensorial element centered quantity to the output.
     */
    virtual void attachTensorElementData(TensorBuffer& buf, std::string name) = 0;

    /*!
     * \brief Finalizes the current writer.
     *
     * This means that everything will be written to disk, except if
     * the onlyDiscard argument is true. In this case only all managed
     * buffers are deleted, but no output is written.
     */
    virtual void endWrite(bool onlyDiscard = false) = 0;
};
} // namespace Ewoms

#endif
