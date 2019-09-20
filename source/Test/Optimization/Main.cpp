#include "Common.h"
#include "Solver.h"
#include "PreProcess.h"
#include "Mesh.h"
#include "Instance.h"

using namespace std;
using namespace Optimization;
//using namespace OpenMesh;
using namespace NumericalMesh;

int main()
{
    /// for different version of debug examples;
    int INSTANCE_MODE = 0;
    /// for quad
    int NX = 1,NY = 1;
    /// the Number of interpolation functions
    int VN = 8;
    const int Dim = 2;
    string fileMesh;
    /// the Number of interpolation functions
    /// switch for different instances;
    switch (INSTANCE_MODE){
        case 0 : {
    ////////////////////////////////////////
    // trival Instance for linear-element //
    ////////////////////////////////////////
    NX = 1,NY = 1;
    VN = 4;
    fileMesh = "../data/meshTrivial.inp";
            break;
        };
        case 1 : {
    /////////////////////////////////////
    // trival Instance for quad-element//
    /////////////////////////////////////
    NX = 1,NY = 1;
    VN = 8;
    fileMesh = "../data/quadratic.inp";
            break;
        };
        default :
        {
            break;
        };
    }
    //string fileMesh = "../data/quadratic.inp";
    //string fileMesh = "../data/meshTrivial.inp";
    //string fileMesh = "../data/stress-concentrate.inp";
    //string fileMesh = "../data/stress-concen-2.inp";

    /// get Rectangle Mesh
    cout<< "---- Generate Mesh ----"<<endl;
    Optimization::Mesh quadMesh;
    GenerateMesh generator(Dim);
    generator.set_NumShape(VN);
    generator.get_MeshFromFile(fileMesh,quadMesh);
    cout<< "---- Generate Mesh : END ----"<<endl;
    cout<< "---- Finite Element Method ----"<<endl;
    ////////////////////
    // Finite Element //
    ////////////////////
    VEC displacement;
    //INSTANCE_Man LinearFE(quadMesh,displacement);
    //INSTANCE2 LinearFE(quadMesh,displacement);
    SolverElesticity LinearFE(quadMesh,displacement);

    LinearFE.set_Geo(NX,NY);
    LinearFE.set_NumShape(VN);
    LinearFE.set_DirichBound();
    LinearFE.set_NeumanBound();
    LinearFE.set_ShapeFunction();
    LinearFE.solve();
    /// logout
    MATRIX strainNode;
    MATRIX stressNode;
    MATRIX strainQuad;
    MATRIX stressQuad;
    LinearFE.get_NodeMatInfo();
    LinearFE.get_Strain_Node(strainNode);
    LinearFE.get_Stress_Node(stressNode);
    LinearFE.get_Strain_Quad(strainQuad);
    LinearFE.get_Stress_Quad(stressQuad);
    /// output file
    string logfile = "output.dat";
    cout<< "---- solver END; write to file: "<< logfile<< " ----"<<endl;
    LinearFE.writeFile(logfile);
    // Solver
    //Optimization::SolverOptimization optimal;
    //optimal.solve();
    /// Debug : End()
    std::cout << "---- End Program ----" << std::endl;
    return 0;
}