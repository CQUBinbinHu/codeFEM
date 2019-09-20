#ifndef NUMERICALMESH_SOLVER
#define NUMERICALMESH_SOLVER

#include "../../Common/math.h"
#include "../../NumericalMesh/core/core.h"
typedef NumericalMesh::Element2D Element;

namespace NumericalMethods {

    class BaseSolver {
        protected:
            int dim;
            int sn;
            MatrixArray solut;
        public :
            BaseSolver(int dim, int sn): dim(dim),sn(sn) {};
            ~BaseSolver(){};

            virtual bool solve() = 0;
            virtual bool writeFile(string filename)const = 0;
            virtual bool readFile(string filename) = 0;
    };
    // Finite Element Solver
    class FESolver : public BaseSolver {
        protected :
            virtual REAL baseFunction(int id,Element&,POINT px);
            virtual VectorXd baseFunctionDerive(int id,Element&, POINT px);
        public :
            FESolver(int dim, int sn ) : BaseSolver(dim,sn) {};
            ~FESolver() {};

            virtual bool writeFile(string filename)const;
            virtual bool readFile(string filename);
    };
    // Finite Difference Solver
    class FDSolver : public BaseSolver {
        public :
            FDSolver(int dim, int sn ) : BaseSolver(dim,sn) {};
            ~FDSolver(){};
    };
};

#endif
