/*
  ====== QuantityOfInterest.hh
  This header file defines the Quantity of Interest for the UQ problem
  and if required the right hand side of the dual problem for this quantity of interest.
  ------
  Example:

  Q is the integral of the pressure over a sub domain of Omega
  sub domain is defined by a ball by center (y) and radius r

  Q = int_Omega k(x)p(x) dV

  where k(x) = 1.0 if dist(x,y) < r else k(x) = 0.0

  ------
  Dr Tim Dodwell - University of Exeter - t.dodwell@exeter.ac.uk
  #1 - last updated 9th May 2016.
*/

#include "QoI.hh" // Include LocalOperator for Quantity of Interest
#include "Cosserat/localStress.hh"
#include "evalTensor.hh"

template <class PARAM,class GFS0, class GFS1, class MBE, class U,class GV,int nedof>
double inline QuantityofInterest(PARAM& param,GFS0& gfs0, GFS1& gfs1, MBE& mbe, U& x,GV& gv, double PatchSize,double Ratio, double tau_y){

  using Dune::PDELab::Backend::native;

  Dune::FieldVector<double,2> y(PatchSize); double radius = 0.5 * PatchSize;

  //  Construct Linear Operator on FEM Space
  typedef Dune::PDELab::FailureCriterion<PARAM,nedof,2> QLOP;
  QLOP Qlop(param,y,radius,Ratio,tau_y);

  // Construct Linear Operator for Stress
  typedef  Dune::PDELab::LocalStress<PARAM,nedof,2> SLOP;
  SLOP Slop(param,y,radius);

  // Construct linear Operator to Evaluate Tensor
  typedef Dune::PDELab::evalTensor<PARAM,2> PLOP;
  PLOP plop(param);

  typedef typename GFS1::template ConstraintsContainer<double>::Type CC;
  CC cc;

  typedef typename GFS0::template ConstraintsContainer<double>::Type C;
  C c;

  typedef Dune::PDELab::GridOperator<GFS1,GFS0,QLOP,MBE,double,double,double,CC,C> QGO;
    QGO qgo(gfs1,cc,gfs0,c,Qlop,mbe);

  typedef Dune::PDELab::GridOperator<GFS1,GFS0,SLOP,MBE,double,double,double,CC,C> SGO;
    SGO sgo(gfs1,cc,gfs0,c,Slop,mbe);

  typedef Dune::PDELab::GridOperator<GFS1,GFS0,PLOP,MBE,double,double,double,CC,C> PGO;
    PGO pgo(gfs1,cc,gfs0,c,plop,mbe);

  typedef typename QGO::Traits::Range V;
    V fc(gfs0,0.0);

  typedef typename SGO::Traits::Range VV;
    VV ss(gfs0,0.0);

  typedef typename PGO::Traits::Range VVV;
    VVV k(gfs0,0.0);


  qgo.residual(x,fc);
  sgo.residual(x,ss);
  pgo.residual(x,k);

  /*Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs0,fc);
  vtkwriter.write("failure_criterion",Dune::VTK::appendedraw);

  Dune::VTKWriter<GV> vtkwriter2(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter2,gfs0,ss);
  vtkwriter2.write("sn",Dune::VTK::appendedraw);

  */

  /*Dune::VTKWriter<GV> vtkwriter3(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter3,gfs0,k);
  vtkwriter3.write("Perm",Dune::VTK::appendedraw);*/

  double Multiplier = fc.infinity_norm(); // Return Maximum

  using Dune::PDELab::Backend::native;
  double S = 0.0;
  for (int i = 0; i < native(ss).size(); i++){
    S += std::abs(native(ss)[i]);
  }
  S /= radius * radius;

  double Q = (S / (Multiplier * tau_y));

  return Q;
}
