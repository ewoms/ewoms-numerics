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
 * \brief Main file of the tutorial problem using the model which assumes
 *        immisciblility.
 */
#include "config.h"              /*@\label{tutorial1:include-begin}@*/
#include <ewoms/numerics/utils/start.hh> /*@\label{tutorial1:include-end}@*/
#include "tutorial1problem.hh" /*@\label{tutorial1:include-problem-header}@*/

int main(int argc, char **argv)
{
    using TypeTag = TTAG(Tutorial1Problem); /*@\label{tutorial1:set-type-tag}@*/
    return Ewoms::start<TypeTag>(argc, argv); /*@\label{tutorial1:call-start}@*/
}
