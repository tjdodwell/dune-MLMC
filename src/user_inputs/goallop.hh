// File  goallop.hh contains linear operator for dual problem.

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
class diffuseDUAL : public NumericalJacobianApplyVolume<diffuseDUAL<LE,GV,dofel>>,
public FullVolumePattern,
public LocalOperatorDefaultFlags,
public InstationaryLocalOperatorDefaultMethods<double>,
public NumericalJacobianVolume<diffuseDUAL<LE,GV,dofel> >
{
public:
    // pattern assembly flags
    enum { doPatternVolume = true };

    // residual assembly flags
    enum { doAlphaVolume = true };
    enum { doLambdaVolume = true };

    diffuseDUAL (GV& gv_, LE& le_,int intorder_=1)
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

    // volume integral depending only on test functions
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

  double f =  le.evaluateF(xip,true);

  double factor =  ip.weight() * eg.geometry().integrationElement(ip.position());

  for (int i = 0; i < lfsv.size(); i++){
          r.accumulate(lfsv,i,-phi[i] * f * factor);
  }

} // end for each quadrature point

}


private:

    const LE& le;
    const GV& gv;
    int intorder;

};

}
}
