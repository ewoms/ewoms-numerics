// $Id$

#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "test_fractionalflow_soilproperties.hh"
#include <dumux/material/fluids/water.hh>
#include <dumux/material/fluids/oil.hh>
#include <dumux/material/twophaserelations.hh>
#include "test_fractionalflow_testproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/diffusion/fv/fvtotalvelocity2p.hh"
#include "dumux/transport/fv/fvsaturationwetting2p.hh"
#include "dumux/transport/fv/capillarydiffusion.hh"
#include "dumux/fractionalflow/impes/impes.hh"
#include <dune/disc/operators/boundaryconditions.hh>
#include "dumux/timedisc/expliciteulerstep.hh"
#include "dumux/fractionalflow/variableclass2p.hh"


int main(int argc, char** argv)
{
    try{
        // define the problem dimensions
        const int dim=2;

        // create a grid object
        typedef double Scalar;
        typedef Dune::SGrid<dim,dim> Grid;
        typedef Grid::LevelGridView GridView;
        typedef Dune::FieldVector<Grid::ctype,dim> FieldVector;

        Dune::FieldVector<int,dim> N(4); N[0] = 30;
        FieldVector L(0);
        FieldVector H(60); H[0] = 300;
        Grid grid(N,L,H);

        grid.globalRefine(0);
        GridView gridView(grid.levelView(0));

        Dune::Water wetmat;
        Dune::Oil nonwetmat;

        Dune::FractionalFlowTestSoil<Grid, Scalar> soil;

        Dune::TwoPhaseRelations<Grid, Scalar> materialLaw(soil, wetmat, nonwetmat);

        typedef Dune::VariableClass<GridView, Scalar> VariableType;

        VariableType variables(gridView);

        typedef Dune::FractionalFlowTestProblem<GridView, Scalar, VariableType> Problem;
        Problem problem(variables, wetmat, nonwetmat, soil, materialLaw,L, H);

        typedef Dune::FVTotalVelocity2P<GridView, Scalar, VariableType, Problem> DiffusionType;
        DiffusionType diffusion(gridView, problem, "pglobal");

        Dune::CapillaryDiffusion<GridView, Scalar, VariableType, Problem> capillaryDiffusion(problem, soil);

        typedef Dune::FVSaturationWetting2P<GridView, Scalar, VariableType, Problem> TransportType;
        TransportType transport(gridView, problem, "vt", capillaryDiffusion);

        int iterFlag = 2;
        int nIter = 30;
        double maxDefect = 1e-5;
        typedef Dune::IMPES<GridView, DiffusionType, TransportType, VariableType> IMPESType;
        IMPESType impes(diffusion, transport, iterFlag, nIter, maxDefect);

        double tStart = 0;
        double tEnd = 5e5;
        const char* fileName = "test_fractionalflow";
        int modulo = 20;
        double cFLFactor = 0.3;
        Dune::TimeLoop<Grid, IMPESType > timeloop(tStart, tEnd, fileName, modulo, cFLFactor);

        timeloop.execute(impes);

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
        return 1;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
        return 1;
    }
}
