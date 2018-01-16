-*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#include "defaultimp.hh"
#include "pattern.hh"
#include "flags.hh"
#include "idefault.hh"

namespace Dune {
namespace PDELab {

template<typename V, typename GV, typename GFS, int dofel>
class errorEstimator : public NumericalJacobianApplyVolume<errorEstimator<V,GV,GFS,dofel>>,
public FullVolumePattern,
public LocalOperatorDefaultFlags,
public InstationaryLocalOperatorDefaultMethods<double>,
public NumericalJacobianVolume<errorEstimator<X,GV,GFS,dofel> >
{
public:
// pattern assembly flags
enum { doPatternVolume = true };

// residual assembly flags
enum { doAlphaVolume = true };

errorEstimator (GV& gv_, V& xc_, GFS& gfs_ int intorder_=2)
: xc(xc_), intorder(intorder_), gv(gv_), gfs(gfs_)
{}

// volume integral depending on test and ansatz functions
template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
{
// dimensions
const int dim = 3;
const int dimw = 3;


const typename GV::IndexSet& is(gv.indexSet());

typedef Dune::PDELab::LocalFunctionSpace<GFS> CLFS;
CLFS clfs(gfs);

typedef Dune::PDELab::LFSIndexCache<CLFS> CLFSCache;
CLFSCache clfsCache(clfs);
std::vector<double> xlc_c(clfs.maxSize());


typedef typename V::template ConstLocalView<CLFSCache> VView;
VView xcView(xc);


// Bind solution x to local element
clfs.bind(*eg);
clfsCache.update();

xcView.bind(clfsCache);
xcView.read(xlc_c);


// extract local function spaces
typedef typename LFSU::template Child<0>::Type LFSU_U1;
typedef typename LFSU::template Child<1>::Type LFSU_U2;
typedef typename LFSU::template Child<2>::Type LFSU_U3;

const LFSU_U1& lfsu_u1 = lfsu.template child<0>();
const LFSU_U2& lfsu_u2 = lfsu.template child<1>();
const LFSU_U3& lfsu_u3 = lfsu.template child<2>();

typedef typename CLFS::template Child<0>::Type LFS1;
typedef typename CLFS::template Child<1>::Type LFS2;
typedef typename CLFS::template Child<2>::Type LFS3;

const LFS1& lfs1 = clfs.template child<0>();
const LFS2& lfs2 = clfs.template child<1>();
const LFS3& lfs3 = clfs.template child<2>();

const unsigned int nodes_per_element = lfsu_u1.size();

// domain and range field type
typedef typename LFSU_U1::Traits::FiniteElementType::
Traits::LocalBasisType::Traits::DomainFieldType DF;
typedef typename R::value_type RF;
typedef typename LFSU_U1::Traits::FiniteElementType::
Traits::LocalBasisType::Traits::JacobianType JacobianType;
typedef typename LFSU_U1::Traits::SizeType size_type;

// select quadrature rule
GeometryType gt = eg.geometry().type();
const QuadratureRule<DF,dim>& rule = QuadratureRules<DF,dim>::rule(gt,intorder);

// Loop over quadrature points
for (typename Dune::QuadratureRule<DF,dim>::const_iterator it = rule.begin(),endit = rule.end(); it != endit; ++it)
{

std::vector<double> phi(nodes_per_element);
lfs1.finiteElement().localBasis().evaluateFunction(it.position(),phi);

Dune::FieldVector<double,3> e(0.0);

for (int i = 0; i < nodes_per_element; i++){
e[0] += (x(lfsu_u1,i) - xlc_c(lfs1,i)) * phi[i];
e[1] += (x(lfsu_u2,i) - xlc_c(lfs2,i)) * phi[i];
e[2] += (x(lfsu_u3,i) - xlc_c(lfs3,i)) * phi[i];
}

double eta = 0.0;

for (int i = 0; i < 3; i++){
eta += std::abs(e[i]);
}


// geometric weight
RF factor = it->weight() * eg.geometry().integrationElement(it->position());


r.accumulate(lfsv,0,eta * factor);


} // end for each quadrature point

} // end alpha_volume


private:

const X& xc;
const GFS& gfs;
const GV& gv;
int intorder;

};

}
}
