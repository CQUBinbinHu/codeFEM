#ifndef NUMERICALMESH_EQUATION
#define NUMERICALMESH_EQUATION

#include "../../Common/common.h"
#include "../../Common/math.h"

namespace NumericalMethods {
///
////////////////////////
    typedef void (*MatFun) (MATRIX&, VectorXd);
    typedef void (*VecFun) (VectorXd&,VectorXd);
    typedef VectorXd (*VecTimeFun) (VectorXd, REAL);
    typedef VectorXd (*ICFUN) (VectorXd);
    typedef VectorXd (*REFUN) (VectorXd,REAL);

    class BaseEquation {

        protected :
            int dim;
            int vars;
        public :
            BaseEquation(int dim, int vars) : dim(dim), vars(vars) {};
            ~BaseEquation(){};
    };

    class PDE : public BaseEquation {

        public :
            PDE(int dim,int vars) : BaseEquation(dim,vars) {};
            ~PDE() {};
    };

    class ConservationLaw : public PDE {

        public :
            ConservationLaw(int dim, int vars) : PDE(dim,vars) {};
            ~ConservationLaw(){};
            
            MatFun Flux;
            VecFun Source;

    };
    /////////////////////////
/////////////////////////////
};

#endif