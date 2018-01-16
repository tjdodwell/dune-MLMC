/*
  ====== AdaptGrid.hh
  This header file defines functions for adaptivity process
  ------
  Dr Tim Dodwell - University of Exeter - t.dodwell@exeter.ac.uk
  #1 - last updated 9th May 2016.
*/



template <typename MYGRID, typename GFS, typename CC, typename BCTypeParam, typename G, typename U, typename U0>
void AdaptMyGrid(MYGRID& grid, GFS& gfs, CC& cc, BCTypeParam& bctype, G& g, U& p, U0& eta) {

  // ==== Adaptively refine mesh.

  double refinement_factor = 0.4;
  double coursening_factor = 0.0;

  double alpha(refinement_factor);       // refinement fraction
  double eta_alpha(0);     // refinement threshold
  double beta(coursening_factor);        // coarsening fraction
  double eta_beta(0);      // coarsening threshold
  int verbose = 0;

  // <<<10>>> Adapt the grid locally...
  Dune::PDELab::element_fraction( eta, alpha, beta, eta_alpha, eta_beta, verbose );
  Dune::PDELab::mark_grid( grid, eta, eta_alpha, eta_beta ,0 , 100, verbose);
  Dune::PDELab::adapt_grid( grid, gfs, p, 2 );

  Dune::PDELab::constraints(bctype,gfs,cc);
  Dune::PDELab::interpolate(g,gfs,p); // this will overwrite the solution !

}
