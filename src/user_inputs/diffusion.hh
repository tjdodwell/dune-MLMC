// -*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-

#ifndef DUNE_PDELAB_DIFFUSION_HH
#define DUNE_PDELAB_DIFFUSION_HH

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

//#include "defaultimp.hh"
//#include "pattern.hh"
#include "flags.hh"
//#include "idefault.hh"

namespace Dune {
    namespace PDELab {

        template<typename LE, typename GV, int dofel>
        class diffuse : public NumericalJacobianApplyVolume<diffuse<LE,GV,dofel>>,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<double>,
        public NumericalJacobianVolume<diffuse<LE,GV,dofel> >
        {
        public:
            // pattern assembly flags
            enum { doPatternVolume = true };

            // residual assembly flags
            enum { doAlphaVolume = true };
            enum { doLambdaVolume = true };

            diffuse (GV& gv_, LE& le_,int intorder_=1,bool Dual = false)
            : le(le_), intorder(intorder_), gv(gv_),isDual(Dual)
            {}

            // volume integral depending on test and ansatz functions
            template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
            void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
            {
                // dimensions
                const int dim = EG::Geometry::mydimension;

                const unsigned int nodel = lfsu.size();

                // domain and range field type
                typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
                typedef typename R::value_type RF;
                typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType JacobianType;
                typedef typename LFSU::Traits::SizeType size_type;

                typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeType;

                // select quadrature rule
                GeometryType gt = eg.geometry().type();
                const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,intorder);

                // Evaluate permeability tensor (assumes it is constant over a single element)




                const typename GV::IndexSet& is(gv.indexSet());

                double Kii = le.evaluateScalar(is.index(eg.entity()));

                Dune::FieldMatrix<double,dim,dim> Kij(0.0);
                for (int i = 0; i < dim; i++) { Kij[i][i] = Kii; }



                // Unwrap solution for element
                FieldVector<double,dofel> p(0.0);
                for (int i=0; i < lfsu.size(); i++){
                    p[i] = x(lfsu,i);
                }



                // Loop over quadrature points
                for (typename Dune::QuadratureRule<DF,dim>::const_iterator it = rule.begin(),endit = rule.end(); it != endit; ++it)
                {

                  // Evaluate Jacobian
                    std::vector<JacobianType> js(nodel);
                    lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);



                    // Transform gradient to real element
                    const typename EG::Geometry::JacobianInverseTransposed jac = eg.geometry().jacobianInverseTransposed(it->position());
                    std::vector<Dune::FieldVector<RF,dim> > gradphi(nodel);

                    for (int i=0; i < lfsu.size(); i++)
                    {
                        gradphi[i] = 0.0;
                        jac.umv(js[i][0],gradphi[i]);
                    }



                    Dune::FieldMatrix<double,dim,dofel> G(0.0);

                    for (int i = 0; i < lfsu.size(); i++){
                      for (int j = 0; j < dim; j++){
                        G[j][i] = gradphi[i][j];
                      }
                    }


                    Dune::FieldVector<double,dim> gradp(0.0),flux(0.0);

                    G.mv(p,gradp); // compute pressure gradient  grad(p) = gradphi * p


                    Kij.mv(gradp,flux); // Compute flux = - Perm * G

                    Dune::FieldVector<double,dofel> residual(0.0);

                    G.mtv(flux,residual); // Compute residual vector res = - G' * Perm * G * p

                    Dune::FieldVector<double,dim> x_global(0.0);

                    // evaluate basis functions
                    std::vector<RangeType> phi(lfsu.size());
                    lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);


                    RF factor = it->weight() * eg.geometry().integrationElement(it->position());

                    for (int i=0; i < lfsu.size(); i++){
                        r.accumulate(lfsv,i,residual[i] * factor);
                    }

                } // end for each quadrature point

            } // end alpha_volume

            template<typename EG, typename LFSV, typename R>
            void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
            {

            // dimensions
            const int dim = EG::Geometry::mydimension;

            // domain and range field type
            typedef typename LFSV::Traits::FiniteElementType::
            Traits::LocalBasisType::Traits::DomainFieldType DF;
            typedef typename R::value_type RF;
            typedef typename LFSV::Traits::FiniteElementType::
            Traits::LocalBasisType::Traits::JacobianType JacobianType;
            typedef typename LFSV::Traits::SizeType size_type;

            typedef typename LFSV::Traits::FiniteElementType::
            Traits::LocalBasisType::Traits::RangeType RangeType;

            // select quadrature rule
            GeometryType gt = eg.geometry().type();
            const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,intorder);


            for (const auto& ip : rule)
            {
              // evaluate basis functions
              std::vector<RangeType> phi(lfsv.size());
              lfsv.finiteElement().localBasis().evaluateFunction(ip.position(),phi);

              // evaluate u at integation points
              Dune::FieldVector<double,dim> xip = eg.geometry().global(ip.position());

              double f =  le.evaluateF(xip,isDual);

              double factor =  ip.weight() * eg.geometry().integrationElement(ip.position());

              for (int i = 0; i < lfsv.size(); i++){
                      r.accumulate(lfsv,i,-phi[i] * f * factor);
              }

            } // end for each quadrature point

            }

        private:

            bool isDual;
            const LE& le;
            const GV& gv;
            int intorder;

        };


        template<typename GV, typename LE, bool dual_problem>
        class error_estimator_residual_based : public NumericalJacobianApplyVolume<error_estimator_residual_based<GV,LE,dual_problem>>,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<double>,
        public NumericalJacobianVolume<error_estimator_residual_based<GV,LE,dual_problem> >
        {
        public:
            // pattern assembly flags
            enum { doPatternVolume = false };
            enum { doPatternSkeleton = false };

            // residual assembly flags
            enum { doAlphaVolume  = true };
            enum { doAlphaSkeleton  = true };

            //! constructor: pass parameter object
            error_estimator_residual_based (GV& gv_,LE& le_)
            : le(le_), gv(gv_)
            {}

            // volume integral depending on test and ansatz functions
            template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
            void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
            {
                // domain and range field type
                typedef typename LFSU::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::DomainFieldType DF;
                typedef typename LFSU::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::RangeFieldType RF;
                typedef typename LFSU::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::RangeType RangeType;
                typedef typename LFSU::Traits::SizeType size_type;

                // dimensions
                const int dim = EG::Geometry::mydimension;
                const int intorder = 2*lfsu.finiteElement().localBasis().order();

                // select quadrature rule
                Dune::GeometryType gt = eg.geometry().type();
                const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

                // loop over quadrature points
                RF sum(0.0);
                for (const auto& ip : rule)
                {
                    // evaluate basis functions
                    std::vector<RangeType> phi(lfsu.size());
                    lfsu.finiteElement().localBasis().evaluateFunction(ip.position(),phi);

                    // evaluate u
                    Dune::FieldVector<double,dim> x_global(0.0);
                    for (size_type i=0; i<lfsu.size(); i++){
                      for (size_type j=0; j<dim; j++){
                        x_global[j] += eg.geometry().corner(i)[j] * phi[i];
                      }
                    }
                    // evaluate right hand side parameter function
                    double f = 1.0;//le.evaluateF(x_global,dual_problem);

                    // integrate f^2
                    RF factor = ip.weight() * eg.geometry().integrationElement(ip.position());
                    sum += f*f * factor;
                }

                // accumulate cell indicator
                DF h_T = diameter(eg.geometry());
                r.accumulate(lfsv,0,h_T*h_T*sum);
            }


            // skeleton integral depending on test and ansatz functions
            // each face is only visited ONCE!
            template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
            void alpha_skeleton (const IG& ig,
                                 const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                 const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                                 R& r_s, R& r_n) const
            {
                // domain and range field type
                typedef typename LFSU::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::DomainFieldType DF;
                typedef typename LFSU::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::RangeFieldType RF;
                typedef typename LFSU::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::JacobianType JacobianType;
                typedef typename LFSU::Traits::SizeType size_type;

                // dimensions
                const int dim = IG::dimension;

                // get cell entities from both sides of the intersection
                auto inside_cell = ig.inside();
                auto outside_cell = ig.outside();

                // evaluate permeability tensors
                const Dune::FieldVector<DF,dim>&
                inside_local = Dune::ReferenceElements<DF,dim>::general(inside_cell.type()).position(0,0);
                const Dune::FieldVector<DF,dim>&
                outside_local = Dune::ReferenceElements<DF,dim>::general(outside_cell.type()).position(0,0);


                // Evaluate elasticity tensor (assumes it is constant over a single element)
                Dune::FieldMatrix<double,dim,dim> K_s(0.0), K_n(0.0);

                const typename GV::IndexSet& is(gv.indexSet());

                le.template evaluateTensor(is.index(inside_cell),K_s);
                le.template evaluateTensor(is.index(outside_cell),K_n);

                // select quadrature rule
                const int intorder = 2*lfsu_s.finiteElement().localBasis().order();
                Dune::GeometryType gtface = ig.geometryInInside().type();
                const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

                // transformation
                typename IG::Entity::Geometry::JacobianInverseTransposed jac;

                // tensor times normal
                const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();
                Dune::FieldVector<RF,dim> Kn_F_s;
                K_s.mv(n_F,Kn_F_s);
                Dune::FieldVector<RF,dim> Kn_F_n;
                K_n.mv(n_F,Kn_F_n);

                // loop over quadrature points and integrate normal flux
                RF sum(0.0);
                for (const auto& ip : rule)
                {
                    // position of quadrature point in local coordinates of elements
                    Dune::FieldVector<DF,dim> iplocal_s = ig.geometryInInside().global(ip.position());
                    Dune::FieldVector<DF,dim> iplocal_n = ig.geometryInOutside().global(ip.position());

                    // evaluate gradient of basis functions
                    std::vector<JacobianType> gradphi_s(lfsu_s.size());
                    lfsu_s.finiteElement().localBasis().evaluateJacobian(iplocal_s,gradphi_s);
                    std::vector<JacobianType> gradphi_n(lfsu_n.size());
                    lfsu_n.finiteElement().localBasis().evaluateJacobian(iplocal_n,gradphi_n);

                    // transform gradients of shape functions to real element
                    jac = inside_cell.geometry().jacobianInverseTransposed(iplocal_s);
                    std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
                    for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);
                    jac = outside_cell.geometry().jacobianInverseTransposed(iplocal_n);
                    std::vector<Dune::FieldVector<RF,dim> > tgradphi_n(lfsu_n.size());
                    for (size_type i=0; i<lfsu_n.size(); i++) jac.mv(gradphi_n[i][0],tgradphi_n[i]);

                    // compute gradient of u
                    Dune::FieldVector<RF,dim> gradu_s(0.0);
                    for (size_type i=0; i<lfsu_s.size(); i++)
                        gradu_s.axpy(x_s(lfsu_s,i),tgradphi_s[i]);
                    Dune::FieldVector<RF,dim> gradu_n(0.0);
                    for (size_type i=0; i<lfsu_n.size(); i++)
                        gradu_n.axpy(x_n(lfsu_n,i),tgradphi_n[i]);

                    // integrate
                    RF factor = ip.weight() * ig.geometry().integrationElement(ip.position());
                    RF jump = (Kn_F_s*gradu_s)-(Kn_F_n*gradu_n);
                    sum += 0.25*jump*jump*factor;
                }

                // accumulate indicator
                // DF h_T = diameter(ig.geometry());
                DF h_T = std::max(diameter(inside_cell.geometry()),diameter(outside_cell.geometry()));
                r_s.accumulate(lfsv_s,0,h_T*sum);
                r_n.accumulate(lfsv_n,0,h_T*sum);
            }


private:
    LE& le;  // two phase parameter class
    const GV& gv;


    template<class GEO>
    typename GEO::ctype diameter (const GEO& geo) const
    {
        typedef typename GEO::ctype DF;
        DF hmax = -1.0E00;
        const int dim = GEO::coorddimension;
        for (int i=0; i<geo.corners(); i++)
        {
            Dune::FieldVector<DF,dim> xi = geo.corner(i);
            for (int j=i+1; j<geo.corners(); j++)
            {
                Dune::FieldVector<DF,dim> xj = geo.corner(j);
                xj -= xi;
                hmax = std::max(hmax,xj.two_norm());
            }
        }
        return hmax;
    }

};

template<typename LE, typename GV, int dofel>
class testError : public NumericalJacobianApplyVolume<testError<LE,GV,dofel>>,
public FullVolumePattern,
public LocalOperatorDefaultFlags,
public InstationaryLocalOperatorDefaultMethods<double>,
public NumericalJacobianVolume<testError<LE,GV,dofel> >
{
public:
    // pattern assembly flags
    enum { doPatternVolume = true };

    // residual assembly flags
    enum { doAlphaVolume = true };

    testError (GV& gv_, LE& le_,int intorder_=4)
    : le(le_), intorder(intorder_), gv(gv_)
    {}

    // volume integral depending on test and ansatz functions
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
    {
        // dimensions
        const int dim = EG::Geometry::mydimension;

        const unsigned int nodel = lfsu.size();

        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename R::value_type RF;
        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef typename LFSU::Traits::SizeType size_type;

        typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeType;


        // Evaluate permeability tensor (assumes it is constant over a single element)

        const typename GV::IndexSet& is(gv.indexSet());

        double Kii = le.evaluateScalar(is.index(eg.entity()));

        Dune::FieldMatrix<double,dim,dim> Kij(0.0);
        for (int i = 0; i < dim; i++) { Kij[i][i] = Kii; }

        // select quadrature rule
        GeometryType gt = eg.geometry().type();
        const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,intorder);

        // Loop over quadrature points
        for (const auto& ip : rule)
        {
          // evaluate basis functions
          std::vector<RangeType> phi(lfsu.size());
          lfsu.finiteElement().localBasis().evaluateFunction(ip.position(),phi);

          // evaluate u at integation points
          double u =0.0;

            for (int i = 0; i < lfsu.size(); i++){
                u += phi[i] * x(lfsu,i);
            }

            Dune::FieldVector<double,dim> xip = eg.geometry().global(ip.position());

            double alpha = 5.0;

            u -= xip[0] * xip[1] * (1.0 - xip[0]) * (1.0  - xip[1]) * std::exp(alpha * xip[0]);

            RF factor = ip.weight() * eg.geometry().integrationElement(ip.position());

            for (int i = 0; i < lfsu.size(); i++){
                  r.accumulate(lfsv,i, u * phi[i] * factor);
            }

        } // end for each quadrature point

        } // end alpha_volume

private:

    const LE& le;
    const GV& gv;
    int intorder;

};

    }
}

template<class GFS, class GV, class X, class COEFF, int dim>
void inline computeFlux(GFS& gfs, GV& gv, X& x, COEFF& le, std::vector<Dune::FieldVector<double,dim>>& flux_all){
/* ===== computeFlux - computes the flux in each element in each direction and returns flux_all a vector of length gv.size(0) containing a vector of length dim.
 ----
 #2 - Last Updated 19th May 2016.
 Dr. Tim Dodwell - t.dodwell@exeter.ac.uk - University of Exeter
 */

  // make local function space
  typedef Dune::PDELab::LocalFunctionSpace<GFS> CLFS;
  CLFS lfsu(gfs);
  typedef Dune::PDELab::LFSIndexCache<CLFS> CLFSCache;
  CLFSCache clfsCache(lfsu);
  std::vector<double> u(lfsu.maxSize());

  typedef typename X::template ConstLocalView<CLFSCache> XView;
  XView xView(x);

  typedef typename GV::IndexSet IndexSet;

  Dune::FieldMatrix<double,dim,dim> Kij(0.0);

  const typename GV::IndexSet& is(gv.indexSet());

  typedef typename CLFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType JacobianType;

  for (const auto& eg : elements(gv)){

    // bind solution x to local element
    lfsu.bind(eg);
    clfsCache.update();
    xView.bind(clfsCache);
    xView.read(u);

    le.template evaluateTensor(is.index(eg),Kij);

    // select quadrature rule
    auto geo = eg.geometry();
    const Dune::QuadratureRule<double,dim>& rule = Dune::QuadratureRules<double,dim>::rule(geo.type(),1);

    // loop over quadrature points
    for (const auto& ip : rule)
    {
      // Evaluate Jacobian
        std::vector<JacobianType> js(lfsu.size());
        lfsu.finiteElement().localBasis().evaluateJacobian(ip.position(),js);

        // Transform gradient to real element
        auto jac = geo.jacobianInverseTransposed(ip.position());
        std::vector<Dune::FieldVector<double,dim> > gradphi(lfsu.size());

        for (int i=0; i < lfsu.size(); i++){
            gradphi[i] = 0.0;
            jac.umv(js[i][0],gradphi[i]);
        }

        Dune::FieldMatrix<double,dim,3> G(0.0);

        for (int i = 0; i < lfsu.size(); i++){
          for (int j = 0; j < dim; j++){
            G[j][i] = gradphi[i][j];
          }
        }
        std::vector<double> gradp(dim), flux(dim);

        G.mv(u,gradp); // compute pressure gradient  grad(p) = gradphi * p
        Kij.mv(gradp,flux); // Compute flux = - Perm * G

        for (int d=0; d < dim; d++){
            flux_all[is.index(eg)][d] = flux[d];
        }

    } // end For each integration point

  } // end For each element


} // End all


void inline evaluateFunction(const Dune::FieldVector<double,2>& xlocal, std::vector<double>& phi){

phi[2] = 4.0 * xlocal[0] * xlocal[1];
phi[1] = 4.0 * (1 - xlocal[0] - xlocal[1]) * xlocal[1];
phi[0] = 4.0 * (1 - xlocal[0] - xlocal[1]) * xlocal[0];

}

void inline evaluateJacobian(const Dune::FieldVector<double,2>& x, std::vector<std::vector<double>>& jac){
for (int i = 0; i < 3; i++){ jac[i].resize(2);}
jac[0][0] = -4.0 * x[1]; jac[0][1] = 4.0 * (1 - x[0] - 2.0 * x[1]);
jac[1][0] = 4.0 * (1 - x[1] - 2.0 * x[0]); jac[1][1] = -4.0 * x[0];
jac[2][0] = 4.0 * x[1]; jac[2][1] = 4.0 * x[0];

}

template<class GFS, class GV, class X, class COEFF, int dim, int dofel>
void inline computeImplicitError(GFS& gfs, GV& gv, X& x, COEFF& le, std::vector<Dune::FieldVector<double,dim>>& flux_all, std::vector<double>& eta_k,  bool isDual = false){
/* ===== computeImplicitError - computes energy norm of the error for each element using an element residual method
 ----
 #3 - 31st August 2016.
 Dr. Tim Dodwell - t.dodwell@exeter.ac.uk - University of Exeter
 ---
 #1 - Implementation of Residual Element Method, as described in 'A posteriori error estimation techniques in practical finite element analysis' by Gratsch & Bathe, Computers & Structures, 2005.

 Problem with implementation as local neuman boundary problem is ill-pose for general boundary conditions. Therefore as expected local element stiffness matrix A is singular, since ker(A) = {[1,1,1]^T} i.e a constant.

 #2 - Issue above solved by solving local neumann boundary problem on a subspace for which the problem is well posed. For this we define error degrees of freedom on each edge of a triangle. Restrict the error (e) to the space spanned by {N1,N2,N3} N1 = 4 * lambda2 * lambda3, N2 = 4 * lambda1 * lambda3, N3 = 4 * lambda1 * lambda2
 where lambda are barycentric coordinates of a triangle.

 #3 - Correction of body force, intoduce boolean variable 'isDual' to control variable  f, calls function of class COEFF called .evalateF(x,isDual)

 */


 typedef typename GV::IndexSet IndexSet;

  Dune::FieldMatrix<double,dim,dim> Kij(0.0);

  const typename GV::IndexSet& is(gv.indexSet());

  for (const auto& eg : elements(gv)){

    // Flux within element

    auto flux_inside = flux_all[is.index(eg)];

    Dune::FieldVector<double,3>  b(0.0); // Initialise right hand side of local error prolem

    

    for (const auto& ig : intersections(gv,eg)){ // Loop through each edge

        if (ig.boundary() == false){ // If element edge is not on the boundary

            auto eg_o = ig.outside(); // Obtain element on the outside

            auto flux_outside = flux_all[is.index(eg_o)];

            // Outward normal to element face
             const Dune::FieldVector<double,dim> n = ig.centerUnitOuterNormal();

            double error_flux = 0.5 * ( (flux_outside * n) - (flux_inside * n));

            Dune::GeometryType gtface = ig.geometryInInside().type();
            const auto& myrule = Dune::QuadratureRules<double,1>::rule(gtface,3);

            for (const auto& ip : myrule){

                double factor = ip.weight() * ig.geometry().integrationElement(ip.position());

                double psi = 4.0 * (1 - ip.position()) * ip.position();

                b[ig.indexInInside()] += psi * error_flux * factor;

            } // end for each integration point

        } // end if interection not on boundary

    } // end of each intersection


    Dune::FieldMatrix<double,3,3> A(0.0);
    Dune::FieldMatrix<double,dim,dim> K(0.0);

    le.template evaluateTensor(is.index(eg), K); // Get Permeability Matrix with element

    auto geo = eg.geometry();

    const Dune::QuadratureRule<double,dim>& rule = Dune::QuadratureRules<double,dim>::rule(geo.type(),2);

    for (const auto& ip : rule){

        std::vector<double> phi(3);

        evaluateFunction(ip.position(),phi);

        // Evaluate Jacobian
        std::vector<std::vector<double>> js(3);
        evaluateJacobian(ip.position(),js);

        // Transform gradient to real element
        auto jac = geo.jacobianInverseTransposed(ip.position());
        std::vector<Dune::FieldVector<double,dim> > gradphi(3);

        for (int i=0; i < 3; i++){
            gradphi[i] = 0.0;
            jac.umv(js[i],gradphi[i]);
        }

        Dune::FieldMatrix<double,dim,dofel> G(0.0);
        Dune::FieldMatrix<double,dofel,dim> GT(0.0);


        for (int i = 0; i < 3; i++){
          for (int j = 0; j < dim; j++){
            G[j][i] = gradphi[i][j];
            GT[i][j] = gradphi[i][j];
          }
        }

        Dune::FieldMatrix<double,dim,3> tmp = K.rightmultiplyany(G);

        double factor = ip.weight() * geo.integrationElement(ip.position());

        Dune::FieldMatrix<double,3,3> Aip = tmp.leftmultiplyany(GT);

        Aip *= factor;

        A += Aip;

        Dune::FieldVector<double,dim> xip = eg.geometry().global(ip.position());

         double f = le.evaluateF(xip,isDual);

        for (int i = 0; i < 3; i++){
          b[i] += f * phi[i] * factor; // Check numbering
        }


    }

    for (const auto& ig : intersections(gv,eg)){ // for each face

      const auto& myrule = Dune::QuadratureRules<double,1>::rule(ig.geometry().type(),1);

      for (const auto& ip : myrule){
        int i = ig.indexInInside();
        auto x = ig.geometry().global(ip.position());
        if (x[0] < 1e-6 || x[0] > 1. - 1e-6 || x[1] < 1e-6 || x[1] > 1. - 1e-6){
            for (int j = 0; j < 3; j++){
              if (i == j){A[i][i] = 1.0; b[i] = 0.0;}
              else {A[i][j] = 0.0; A[j][i] = 0.0;}
            } // end for each j
        } // end if on Dirichlet boundary
     } // End loop through integration points

   } // end loop through each edge

    Dune::FieldVector<double,dofel> e(0.0);

    A.solve(e,b);

    eta_k[is.index(eg)] = e * b; // Energy Norm


  }








} // End all


#endif
