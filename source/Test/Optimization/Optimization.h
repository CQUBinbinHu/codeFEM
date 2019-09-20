#ifndef OPTIMIZATION_OPTIMIZATION
#define OPTIMIZATION_OPTIMIZATION

#include "Common.h"
#include "Solver.h"
#include "Instance.h"
/*
 * Solve the unconstraint optimization model
 */
namespace Optimization {
    //////////////////
/// start of namespace
/////////////////////////
class SolverOptimization {

private:
    const double eps = 1e-2;
    int PN = 4;
    double x = sqrt(0.);

protected:
    unsigned int EN = 0;
    
    Mesh& mesh;
    MATRIX QuadPoints;
    VEC quadParameters;
    vector<VEC> ShapeFunctionsQuad;
    vector<MATRIX> Grad1ShapeFunctionQuad;
    vector<MATRIX> Grad2ShapeFunctionQuad;

    VEC displacement;
    VEC level;
    VEC levelDescent;
    /// get information about Shape Functions and gradients 
    void getShapeFunction(VEC& shapeInfo,VEC point);
    void getShapeFunctionGrad(MATRIX& shapeInfo,VEC point);
    void getShapeFunctionGrad2(MATRIX& shapeInfo,VEC point);
    void getDescent();
    void SearchLine(double& delta);
    double getObjection(VEC& valLevel);
    double HeasiveFun(double);
    double DeltaFun(double);
    SolverElesticity& FEM;
public:
    ~SolverOptimization(){};
    SolverOptimization(Mesh& i_mesh,SolverElesticity& fem):mesh(i_mesh),FEM(fem){
        EN = mesh.elements.size();
        level.resize(mesh.points.size());
        levelDescent.resize(mesh.points.size());
    };

    bool Initialization();
    bool solve();

    /// IO set
    bool set_level();
    bool set_displace();
    bool set_Objection();
    /// IO get
    bool get_level();
    bool get_displace();
    bool get_Objection();
};
    //////////////////
/// end of namespace
/////////////////////////
}
#endif
