#include "Solver.h"
using namespace NumericalMethods;
/// base functions for finite element
REAL FESolver::baseFunction (int id, Element& ele, POINT po) {

    REAL basefun = 1.;
    switch(id) {
        case 0 : 
            {
                basefun = 1.;
                break;
            };
        case 1 :
            {
                REAL x0 = ele.midPoint[0];
                basefun = 2*(po[0]-x0)/ele.Size();
                break;
            };         
        case 2 :
            {
                REAL x0 = ele.midPoint[0];
                basefun = 2 * (po[0]-x0)/ele.Size();
                basefun = 0.5 * basefun * basefun;
                break;
            };
        default :
            {
                break;
            }
    }
    return basefun;
};
/// derivative of base functions
VectorXd FESolver::baseFunctionDerive (int id,Element& ele, POINT po) {
    VectorXd vec(dim);
    vec[0] = 0.;
    switch(id) {
        case 0 :
            {
                vec[0] = 0;
                break;
            };
        case 1 :
            {
                vec[0] = 2. / ele.Size();
                break;
            };
        case 2 :
            {
                REAL x0 = ele.midPoint[0];
                vec[0] = (2.*(po[0]-x0)/ele.Size()) * (2./ele.Size());
                break;
            };
        default :
            {
                break;
            };
    }
    return vec;
};