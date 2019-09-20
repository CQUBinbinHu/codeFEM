#ifndef ROBUSTOTIMIZATION_SOLVER
#define ROBUSTOTIMIZATION_SOLVER

#include "Common.h"
#include "BaseSolver.h"
#include "Mesh.h"
#include "../../NumericalMethods/core/core.h"
using namespace NumericalMethods;

namespace  Optimization {
/*
 * finite element solver for isotropic linear elasticity
 */
class SolverElesticity : public BaseFEM{

    private:

    protected :
        const REAL eps = 1e-3;
        INTE EN;
        INTE DISN;
        /// for Mesh
        int NX = 0;
        int NY = 0;
        /// material
        REAL valuePoisson = 0.3;
        REAL valueElesticity = 1.;
        REAL poisson = valuePoisson;
        REAL elesticity = valueElesticity;
        /// Global
        VEC load;
        MatrixXd stiffness;
        MatrixXd DisplaceStiffness;
        MatrixXd DirichletStiffness;
        //local
        MatrixArray stiffnessLocal;
        vector<VEC> loadLocal;
        /// for Dirichlet boundary condition
        vector<bool> ifDirichlet;
        vector<INTE> idDirichlet;
        VEC valDirichlet;
        VEC loadDiri;
        /// for Neumann Boundary condition
        vector<bool> ifNeumann;
        vector<INTE> idNeuman;
        VEC valNeuman;
        /// base function
        ///
        virtual void LocalStiffness();
        virtual void LocalLoad();
        virtual void DirichBound();
        virtual void NeumanBound();
        virtual void GeneralLinearEquation();

        virtual void basefunction(MATRIX&,MATRIX&,VEC);
    public :
        ~SolverElesticity(){};
        SolverElesticity( Mesh& mesh, VEC& disp ) \
            : BaseFEM(mesh,disp){
        
                MD.resize(3,3);
                MD(0,0) = 1.;
                MD(0,1) = poisson;
                MD(0,2) = 0.;
                MD(1,0) = poisson;
                MD(1,1) = 1.;
                MD(1,2) = 0.;
                MD(2,0) = 0.;
                MD(2,1) = 0.;
                MD(2,2) = 0.5*(1-poisson);
                MD *= (elesticity / (1. - poisson*poisson));

                displacement.resize(2*mesh.points.size());
                DISN = displacement.size();
                for(int di=0;di<DISN;di++){
                    displacement[di] = 0.;
                }
                stiffness.resize(DISN,DISN);
                for(int ni=0;ni<DISN;ni++){
                    for(int mi=0;mi<DISN;mi++){
                        stiffness(mi,ni) = 0.;
                    }
                }
        };

        virtual bool solve();
        /// IO
        virtual bool readFile(string filename);
        virtual void set_mesh(Mesh& me){mesh = me;};
        virtual void set_Poisson(REAL poisson){valuePoisson = poisson;};
        virtual void set_DirichBound();
        virtual void set_NeumanBound();
        virtual void set_ShapeFunction();
        void set_Geo(int nx,int ny){ NX=nx;NY=ny; };
        
        vector<INTE> get_idNeuman(){return idNeuman;};
        VEC get_valNeuman(){return valNeuman;};
};
////////////////////////
/// end of namespace
////////////////////////
};
#endif