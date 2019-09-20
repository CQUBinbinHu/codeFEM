#ifndef OPTIMIZATION_BASESOLVER
#define OPTIMIZATION_BASESOLVER

#include "Common.h"
#include "Mesh.h"

using namespace NumericalMesh;
using namespace Optimization;


namespace  Optimization {

class BaseFEM {
protected:
    const int Dim = 2;
    int NumShape = 4;
    Mesh& mesh;
    VEC& displacement; // unknown

    MatrixXd MD;
    VEC quadParameters;
    vector<MATRIX> JacobMats;
    vector<MATRIX> ShapeMats;
    vector<MATRIX> JacobMatsNode;
    vector<MATRIX> ShapeMatsNode;

    bool ifStress = false;
    bool ifStrain = false;
    vector<MATRIX> StressMats;
    vector<MATRIX> StrainMats;
    vector<MATRIX> StressMatsNode;
    vector<MATRIX> StrainMatsNode;
public:
    BaseFEM(Mesh& me,VEC& disp): mesh(me),displacement(disp){};

    /// IO
    virtual bool writeFile(string)const;
    virtual bool readFile (string) = 0;
    virtual void set_NumShape(int num){ NumShape = num; };
    virtual void get_Stress_Quad(MATRIX&);
    virtual void get_Strain_Quad(MATRIX&);
    virtual void get_Stress_Node(MATRIX&);
    virtual void get_Strain_Node(MATRIX&);
    virtual void get_NodeMatInfo();

    vector<MATRIX> get_JacobMats(){return JacobMats;};
    vector<MATRIX> get_ShapeMats(){return ShapeMats;};
};

}
#endif
