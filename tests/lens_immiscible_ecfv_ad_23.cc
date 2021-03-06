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
 * \brief Two-phase test for the immiscible model which uses the element-centered finite
 *        volume discretization in conjunction with automatic differentiation
 */
#include "config.h"

#include "lens_immiscible_ecfv_ad.hh"

#include <dune/grid/geometrygrid.hh>
#include <dune/grid/io/file/dgfparser/dgfgeogrid.hh>

BEGIN_PROPERTIES

// Use Dune-grid's GeometryGrid< YaspGrid >
SET_PROP(LensProblemEcfvAd, Grid )
{
  template< class ctype, unsigned int dim, unsigned int dimworld >
  class IdentityCoordFct
    : public Dune::AnalyticalCoordFunction
      < ctype, dim, dimworld, IdentityCoordFct< ctype, dim, dimworld > >
  {
    using This = IdentityCoordFct< ctype, dim, dimworld >;
    using Base = Dune::AnalyticalCoordFunction< ctype, dim, dimworld, This >;

  public:
    using DomainVector = typename Base :: DomainVector;
    using RangeVector = typename Base :: RangeVector ;

    template< typename... Args >
    IdentityCoordFct( Args&... )
    {}

    RangeVector operator()(const DomainVector& x) const
    {
      RangeVector y;
      evaluate( x, y );
      return y;
    }

    void evaluate( const DomainVector &x, RangeVector &y ) const
    {
      y = 0;
      for( unsigned int i = 0; i<dim; ++i )
        y[ i ] = x[ i ];
    }

  };

  using MyYaspGrid = Dune::YaspGrid< 2 >;

public:
  //using type = MyYaspGrid;
  typedef Dune::GeometryGrid< MyYaspGrid,
                              IdentityCoordFct< typename MyYaspGrid::ctype,
                                                MyYaspGrid::dimension,
                                                MyYaspGrid::dimensionworld+1> >  type;
};

END_PROPERTIES

#include <ewoms/numerics/utils/start.hh>

int main(int argc, char **argv)
{
    using ProblemTypeTag = TTAG(LensProblemEcfvAd);
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}
