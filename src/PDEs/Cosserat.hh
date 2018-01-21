#ifndef Cosserat_h
#define Cosserat_h

#include <vector>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>

using namespace std;

namespace Dune {
    namespace PDELab {

template<typename PARAM,int nedof>
class Cosserat:
        public Dune::PDELab::NumericalJacobianApplyVolume<Cosserat<PARAM,nedof>>,
        public Dune::PDELab::NumericalJacobianVolume<Cosserat<PARAM,nedof>>,
        public Dune::PDELab::FullVolumePattern,
        public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
// pattern assembly flags
    enum { doPatternVolume = true };

// residual assembly flags
    enum { doAlphaVolume = true };

Cosserat(PARAM& param_,
          unsigned int intorder_ = 2):
    param(param_), intorder(intorder_)
{}

// Volume Integral Depending on Test and Ansatz Functions
template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
{

    // Galerkin Finite Elements - assumes lfsu == lfsv
    const int dim = 2;

    // Unwrap shape local function spaces
    typedef typename LFSU::template Child<0>::Type LFSU_U1;
    typedef typename LFSU::template Child<1>::Type LFSU_U2;
    typedef typename LFSU::template Child<2>::Type LFSU_ROT3;

    const LFSU_U1& lfsu_u1 = lfsu.template child<0>();
    const LFSU_U2& lfsu_u2 = lfsu.template child<1>();
    const LFSU_ROT3& lfsu_rot3 = lfsu.template child<2>();

    const int p1_n = lfsu_rot3.size();
    const int p2_n = lfsu_u1.size();

    // domain and range field type
    typedef typename LFSU_U1::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU_U1::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainType D;
    typedef typename R::value_type RF;
    typedef typename LFSU_U1::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType JacobianTypeU;
    typedef typename LFSU_U1::Traits::SizeType size_type;

    typedef typename LFSU_ROT3::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType JacobianTypeR;

    typedef typename LFSU_ROT3::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType Range;

    // select quadrature rule

    auto geo =  eg.geometry();
    GeometryType gt = geo.type();
    const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,intorder);

    // Unwrap solution at node, into vector d

    Dune::FieldVector<double,nedof> d(0.0);

    for (size_type i=0; i < p2_n; ++i){
        d[i] = x(lfsu_u1,i); // U1
        d[i + p2_n] = x(lfsu_u2,i); // U2
        if (i < p1_n){
          d[i + 2 * p2_n] = x(lfsu_rot3,i); // ROT3
        }
    }

    double diameter = 7.0e-3;



    // Loop over quadrature points

    for (typename Dune::QuadratureRule<DF,dim>::const_iterator it = rule.begin(),endit = rule.end(); it != endit; ++it)
    {

        // Evaluate shape functions at Integration Point
        std::vector<Range> phi_P1(lfsu_rot3.size());
        lfsu_rot3.finiteElement().localBasis().evaluateFunction(it->position(),phi_P1);

        // Evaluate gradients of shape function at integration point
        std::vector<JacobianTypeR> js_P1(p1_n);
        lfsu_rot3.finiteElement().localBasis().evaluateJacobian(it->position(),js_P1);
        std::vector<JacobianTypeU> js_P2(p2_n);
        lfsu_u1.finiteElement().localBasis().evaluateJacobian(it->position(),js_P2);

        // Transform gradient to real element
        const typename EG::Geometry::JacobianInverseTransposed jac = eg.geometry().jacobianInverseTransposed(it->position());
        std::vector<Dune::FieldVector<RF,dim> > p2_gradphi(p2_n);
        std::vector<Dune::FieldVector<RF,dim> > p1_gradphi(p1_n);

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

        Dune::FieldVector<double,6> e(0.0), s(0.0);

        B.mv(d,e);

        Dune::FieldVector<double,2> xg = eg.geometry().global(it->position());

        Dune::FieldMatrix<double,6,6> C = param.evaluateMatrix(xg);

        C.mv(e,s);

        Dune::FieldVector<double,nedof> res(0.0);

        B.mtv(s,res); // res_u = Bt * sig;

        double factor = it->weight() * eg.geometry().integrationElement(it->position());


            for (size_type i = 0; i < p2_n; i++){
               r.accumulate( lfsu_u1, i, res[i] * factor);
               r.accumulate( lfsu_u2, i, res[i + p2_n] * factor);
            }

            for (size_type i=0; i < p1_n; i++){
                r.accumulate( lfsu_rot3, i, res[i + 2 * p2_n] * factor);
            }


    } // end for each quadrature point


} // end alpha_volume



        private:

            PARAM& param;
            int intorder;
};

    }
}

#endif