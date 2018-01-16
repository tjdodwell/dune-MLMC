// -*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-


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

        class Ltwo : public NumericalJacobianApplyVolume<Ltwo>,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<double>,
        public NumericalJacobianVolume<Ltwo>
        {
        public:
            // pattern assembly flags
            enum { doPatternVolume = true };

            // residual assembly flags
            enum { doAlphaVolume = true };

            Ltwo (int intorder_=1)
            : intorder(intorder_)
            {}

            // volume integral depending on test and ansatz functions
            template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
            void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
            {
                // dimensions
                const int dim = EG::Geometry::mydimension;

                // domain and range field type
                typedef typename LFSU::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::DomainFieldType DF;
                typedef typename R::value_type RF;
                typedef typename LFSU::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::JacobianType JacobianType;
                typedef typename LFSU::Traits::SizeType size_type;

                typedef typename LFSU::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::RangeType RangeType;

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
                  std::vector<double> u(dim);

                  for (int j = 0; j < dim; j++){
                    u[j] = 0.0;
                    for (int i = 0; i < lfsu.size(); i++){
                        u[j] += phi[i] * x(lfsu,i);
                    }
                  }



                    RF factor = ip.weight() * eg.geometry().integrationElement(ip.position());

                    for (int i = 0; i < lfsu.size(); i++){
                          r.accumulate(lfsv,i,u[i] * phi[i] * factor);
                    }

                } // end for each quadrature point

            } // end alpha_volume


        private:
            int intorder;

        };


}
}
