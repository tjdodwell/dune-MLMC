// -*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>


namespace Dune {
    namespace PDELab {

        template <class PARAM, int dim>
        class evalTensor : public NumericalJacobianApplyVolume<evalTensor<PARAM,dim> >,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<double>,
        public NumericalJacobianVolume<evalTensor<PARAM,dim> >
        {
        public:
            // pattern assembly flags
            enum { doPatternVolume = true };

            // residual assembly flags
            enum { doAlphaVolume = true };

            evalTensor(PARAM& param_,int intorder_=2)
            : param(param_),intorder(intorder_)
            {}

            // volume integral depending on test and ansatz functions
            template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
            void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
            {

              // select quadrature rule
              auto geo = eg.geometry();
              const QuadratureRule<double,dim>& rule = QuadratureRules<double,dim>::rule(geo.type(),intorder);

                // Loop over quadrature points
                for (const auto& ip : rule){
                  Dune::FieldVector<double,2> xg = eg.geometry().global(ip.position());
                  double Phi_ip = param.evaluateMisalignment(xg);
                  r.accumulate(lfsv,0,Phi_ip * ip.weight());
                } // end for each Intgeration Point


            } // end alpha_volume


        private:
            int intorder;
            PARAM& param;

        };


}
}
