#ifndef NUMERICALMESH_CORE_MESH
#define NUMERICALMESH_CORE_MESH

#include "../../Common/math.h"
#include "../../Common/common.h"

namespace NumericalMesh {
    class Element2D {
        private:
            REAL size;
        public:
        ~Element2D(){};
        INTE id;
        vector<INTE> vertex;
        POINT2 midPoint;
        REAL size;
        vector<POINT2> nvec;

        MATRIX quadPoints;
        VEC quadParameters;

        vector<MATRIX> quadPointsBound;
        vector<VEC> quadPointsBound;
        vector<int> NeiEleID;
    };

    class BaseMesh {
        public:
        ~BaseMesh(){};
        vector<POINT2> points;
    };
    class FEMesh : public BaseMesh {
        public :
        ~FEMesh(){};

        //
        vector<INTE> innerEleID;
        vector<INTE> boundEleID;
        vector<Element2D> elements;
    };
};

#endif