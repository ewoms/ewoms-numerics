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
 * \brief This test is identical to the simulation of the lens problem that uses the
 *        element centered finite volume discretization in conjunction with automatic
 *        differentiation (lens_immiscible_ecfv_ad).
 *
 * The only difference is that it uses multiple compile units in order to ensure that
 * eWoms code can be used within libraries that use the same type tag within multiple
 * compile units. This file calls contains main() and just calls the main entry point
 * defined by the first compile unit.
 */
#include "config.h"

int mainCU1(int argc, char **argv);
int mainCU2(int argc, char **argv);

int main(int argc, char **argv)
{
    return mainCU1(argc, argv);
}
