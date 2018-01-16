// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<fstream>

#include<vector>
#include<map>
#include<string>
#include<sstream>
#include <mpi.h>

#include <config.h>


#include <dune/common/parametertree.hh>
Dune::ParameterTree config;

#include <dune/common/bitsetvector.hh>

#include <dune/grid/yaspgrid.hh> // Checked Inclusion
#include <dune/grid/common/gridview.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/gridinfo.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/io.hh>
#include <dune/istl/matrixmarket.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/collectivecommunication.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>




#include <dune/pdelab/newton/newton.hh>
#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/finiteelementmap/rannacherturekfem.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/common/instationaryfilenamehelper.hh>
#include <dune/pdelab/instationary/implicitonestep.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
//#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/instationary/onestepparameter.hh>

#include <dune/pdelab/adaptivity/adaptivity.hh>

// -- Problem Specific Includes

#include "user_inputs/Random_Coeff/random_coeff.hh"
#include "user_inputs/user_model_cosserat.hh"

#include "user_inputs/onesolution.hh"

int main(int argc, char** argv)
{
  
      MPI_Init(&argc,&argv);

      //Read ini file
      Dune::ParameterTreeParser parser;
      parser.readINITree(argv[1],config);

      const int dim = 2;

      std::vector<double> model_param(15);

      const int Coarse_Elements = config.get<int>("MLMC.CoarseMeshElements",8);

      Dune::FieldVector<double,dim> L(config.get<double>("Domain.L",1.0)); // Dimension of Grid (Square)
      Dune::array<int,dim> N(Dune::fill_array<int,2>(Coarse_Elements));
      std::bitset<dim> periodic(false);
      int overlap=0;

      typedef Dune::YaspGrid<dim> GRID;

      GRID grid(L,N,periodic,overlap);

      int maxLevel = config.get<int>("maxLevel",5);

      grid.globalRefine(maxLevel);

      typedef MODEL<dim,GRID> MY_MODEL;

      MY_MODEL model(grid,L);

      testSamples<MY_MODEL> myTest(model);

      myTest.apply(5);



}
