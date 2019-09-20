#ifndef ROBUSTOPTIMIZATION_INSTANCE
#define ROBUSTOPTIMIZATION_INSTANCE

#include "Solver.h"

class INSTANCE_Man : public SolverElesticity{

public:
    INSTANCE_Man( Mesh& mesh,VEC& vecPlace) \
        : SolverElesticity(mesh,vecPlace){};
    virtual bool set_BoundFromFile(string filename);
    virtual void set_DirichBound();
    virtual void set_NeumanBound();

};

/*
 * sysmetric boundary bondary condition
 */
class INSTANCE2: public SolverElesticity{

public:
    INSTANCE2( Mesh& mesh,VEC& vecPlace) \
        : SolverElesticity(mesh,vecPlace){};
    virtual void set_DirichBound();
    virtual void set_NeumanBound();
};

#endif

