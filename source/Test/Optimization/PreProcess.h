#ifndef ROBUSTOPTIMIZATION_PREPROCESS
#define ROBUSTOPTIMIZATION_PREPROCESS

#include "Common.h"
#include "Mesh.h"
using namespace NumericalMesh;

class GenerateMesh {

    private:
        int Dim;
        vector<POINT> Points;
        int NVert = 4;
        vector<int> VecNx;
    public:
        GenerateMesh(int dim) :Dim(dim) {
            /// debug
            cout<< "GenerateMesh()"<<endl;
        };

        ~GenerateMesh(){};

        /// IO
        bool get_Mesh(Mesh&);
        bool get_MeshFromFile(string,Mesh&);
        bool set_Geometry(vector<POINT> geometry,vector<int> Nxy);
        bool set_GeometryFromFile(string);
        void set_NumShape(int num){ NVert = num; };
};

#endif


