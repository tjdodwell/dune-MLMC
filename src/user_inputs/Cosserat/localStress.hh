// -*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>


namespace Dune {
    namespace PDELab {

        template <class PARAM, int nedof, int dim>
        class LocalStress : public NumericalJacobianApplyVolume<LocalStress<PARAM,nedof,dim> >,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<double>,
        public NumericalJacobianVolume<LocalStress<PARAM,nedof,dim> >
        {
        public:
            // pattern assembly flags
            enum { doPatternVolume = true };

            // residual assembly flags
            enum { doAlphaVolume = true };

            LocalStress(PARAM& param_,Dune::FieldVector<double,dim>& y_, double radius_,int intorder_=2)
            : param(param_),y(y_),radius(radius_),intorder(intorder_)
            {}

            // volume integral depending on test and ansatz functions
            template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
            void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
            {
              // This local operator with map from Vh(solution space) to V0 (Space of piecewise constants on a mesh Qh)
              // The operation M*u will return a piecewise constant funciton, containing the maximum value of the failure
              // criterion with each element.

              // Unwrap shape local function spaces - for composite grid function space
              typedef typename LFSU::template Child<0>::Type LFSU_U1;
              typedef typename LFSU::template Child<1>::Type LFSU_U2;
              typedef typename LFSU::template Child<2>::Type LFSU_ROT3;

              const LFSU_U1& lfsu_u1 = lfsu.template child<0>();
              const LFSU_U2& lfsu_u2 = lfsu.template child<1>();
              const LFSU_ROT3& lfsu_rot3 = lfsu.template child<2>();

              const int p1_n = lfsu_rot3.size();
              const int p2_n = lfsu_u1.size();

              // Unwrap element solution into a (column) vector d

              Dune::FieldVector<double,nedof> d(0.0);

              for (int i=0; i < p2_n; ++i){
                  d[i] = x(lfsu_u1,i); // U1
                  d[i + p2_n] = x(lfsu_u2,i); // U2
                  if (i < p1_n){
                    d[i + 2 * p2_n] = x(lfsu_rot3,i); // ROT3
                  }
              }

              typedef typename LFSU_U1::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType JacobianTypeU;
              typedef typename LFSU_ROT3::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType JacobianTypeR;
              typedef typename LFSU_ROT3::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType Range;

              // select quadrature rule
              auto geo = eg.geometry();
              const QuadratureRule<double,dim>& rule = QuadratureRules<double,dim>::rule(geo.type(),intorder);

              double sn = 0.0;

                // Loop over quadrature points
                for (const auto& ip : rule)
                {
                  // Evaluate shape functions at Integration Point
                  std::vector<Range> phi_P1(lfsu_rot3.size());
                  lfsu_rot3.finiteElement().localBasis().evaluateFunction(ip.position(),phi_P1);

                  // Evaluate gradients of shape function at integration point
                  std::vector<JacobianTypeR> js_P1(p1_n);
                  lfsu_rot3.finiteElement().localBasis().evaluateJacobian(ip.position(),js_P1);
                  std::vector<JacobianTypeU> js_P2(p2_n);
                  lfsu_u1.finiteElement().localBasis().evaluateJacobian(ip.position(),js_P2);

                  // Transform gradient to real element
                  const typename EG::Geometry::JacobianInverseTransposed jac = eg.geometry().jacobianInverseTransposed(ip.position());
                  std::vector<Dune::FieldVector<double,dim> > p2_gradphi(p2_n);
                  std::vector<Dune::FieldVector<double,dim> > p1_gradphi(p1_n);

                  for (int i=0; i < p2_n; i++){

                      if (i < p1_n){
                          p1_gradphi[i] = 0.0;
                          jac.umv(js_P1[i][0],p1_gradphi[i]);
                      }
                      p2_gradphi[i] = 0.0;
                      jac.umv(js_P2[i][0],p2_gradphi[i]);
                  }

                  // (A) Compute Elastic part of material

                  // Compute Element Strains from nodal displacements
                  Dune::FieldMatrix<double,6,nedof> B(0.0);
                  for (int i = 0; i < p2_n; i++){
                      B[0][i] = p2_gradphi[i][0]; // E11
                      B[1][i + p2_n] = p2_gradphi[i][1]; // E22
                      B[2][i] = p2_gradphi[i][1];   // E12
                      B[3][i + p2_n] = p2_gradphi[i][0];   // E21
                  }
                  // Rotation Terms
                  for (int i = 0; i < p1_n; i++){
                    B[2][i + 2 * p2_n] = phi_P1[i]; // E12
                    B[3][i + 2 * p2_n] = -phi_P1[i]; // E21
                    B[4][i + 2 * p2_n] = p1_gradphi[i][0]; // K31
                    B[5][i + 2 * p2_n] = p1_gradphi[i][1]; // K32
                  }

                  Dune::FieldVector<double,6> e(0.0), s(0.0); // Initialise Vectors for Cosserat Strains and Stresses respectively

                  B.mv(d,e); // Compute Cosserat Strain at Integration Point

                  Dune::FieldVector<double,2> xg = eg.geometry().global(ip.position());

                  Dune::FieldMatrix<double,6,6> C = param.evaluateMatrix(xg);

                  C.mv(e,s); // Compute Stress

                  double factor = ip.weight() * eg.geometry().integrationElement(ip.position());

                  // Compute Failure criterion at Intergration Point - Quadratic Yield criterion
                  if (std::abs(xg[0] - y[0]) < radius && std::abs(xg[1] - y[1]) < radius){

                  //  std::cout << "Strain" << e << std::endl;
                  //  std::cout << "Stress" << s << std::endl;
                    r.accumulate(lfsv,0,s[0]*factor);
                  } // If Integration Point is with Sub-Region

                } // end for each Intgeration Point


                // Write Max_failure_c to element



            } // end alpha_volume


        private:
            int intorder;
            Dune::FieldVector<double,dim>& y;
            double radius;
            PARAM& param;

        };


}
}
