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
 * \copydoc Ewoms::EcfvBaseOutputModule
 */
#ifndef EWOMS_ECFV_VTK_BASE_OUTPUT_MODULE_HH
#define EWOMS_ECFV_VTK_BASE_OUTPUT_MODULE_HH

#include "ecfvproperties.hh"

#include <ewoms/numerics/io/vtkmultiwriter.hh>

namespace Ewoms {
/*!
 * \ingroup EcfvDiscretization
 *
 * \brief Implements the discretization specific parts of writing files.
 */
template<class TypeTag>
class EcfvBaseOutputModule
{
public:
    typedef BaseOutputWriter::Scalar Scalar;
    typedef BaseOutputWriter::Vector Vector;
    typedef BaseOutputWriter::ScalarBuffer ScalarBuffer;
    typedef BaseOutputWriter::VectorBuffer VectorBuffer;
    typedef BaseOutputWriter::TensorBuffer TensorBuffer;

    /*!
     * \brief Add a buffer where the data is associated with the
     *        degrees of freedom to the current VTK output file.
     */
    static void attachScalarDofData_(BaseOutputWriter& baseWriter,
                                     ScalarBuffer& buffer,
                                     const std::string& name)
    { baseWriter.attachScalarElementData(buffer, name.c_str()); }

    /*!
     * \brief Add a buffer where the data is associated with the
     *        degrees of freedom to the current VTK output file.
     */
    static void attachVectorDofData_(BaseOutputWriter& baseWriter,
                                     VectorBuffer& buffer,
                                     const std::string& name)
    { baseWriter.attachVectorElementData(buffer, name.c_str()); }

    /*!
     * \brief Add a buffer where the data is associated with the
     *        degrees of freedom to the current VTK output file.
     */
    static void attachTensorDofData_(BaseOutputWriter& baseWriter,
                                     TensorBuffer& buffer,
                                     const std::string& name)
    { baseWriter.attachTensorElementData(buffer, name.c_str()); }
};

} // namespace Ewoms

#endif
