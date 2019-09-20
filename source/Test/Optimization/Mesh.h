#ifndef OPTIMIZATION_MESH
#define OPTIMIZATION_MESH

#include "Common.h"
#include "../../NumericalMesh/core/core.h"

using namespace NumericalMesh;
typedef NumericalMesh::Element2D Element;

/// namespace Optimization
namespace Optimization{
/// start of namespace
////////////////////////////////////////////////////////
class Mesh : public NumericalMesh::FEMesh {
    public :
        ~Mesh(){}
        vector<MATRIX> contralPoints;
};
/// end of namespace
////////////////////////////////////////////////////////
};
#endif