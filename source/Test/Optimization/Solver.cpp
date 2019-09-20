#include "Solver.h"
using namespace Optimization;
//#include "../../core/NumericalMethods/NumericalMethods/Tools.h"
void SolverElesticity::set_ShapeFunction(){

    int GN = 3;
    NumerIntegral GaussInte(GN,'L');

    VEC x,y;
    x.resize(2);
    y.resize(2);

    x[0] = -1.;
    x[1] = +1.;
    y[0] = -1.;
    y[1] = +1.;

    MATRIX nodeMat(2,2*NumShape);
    switch(NumShape){
    case 4 : {

        nodeMat(0,0) = -1.;
        nodeMat(0,1) = +1.;
        nodeMat(0,2) = +1.;
        nodeMat(0,3) = -1.;
        nodeMat(1,0) = -1.;
        nodeMat(1,1) = -1.;
        nodeMat(1,2) = +1.;
        nodeMat(1,3) = +1.;

               break;
           }
    case 8 : {
        nodeMat(0,0) = -1.;
        nodeMat(0,1) =  0.;
        nodeMat(0,2) = +1.;
        nodeMat(0,3) = +1.;
        nodeMat(0,4) = +1.;
        nodeMat(0,5) =  0.;
        nodeMat(0,6) = -1.;
        nodeMat(0,7) = -1.;
        nodeMat(1,0) = -1.;
        nodeMat(1,1) = -1.;
        nodeMat(1,2) = -1.;
        nodeMat(1,3) =  0.;
        nodeMat(1,4) = +1.;
        nodeMat(1,5) = +1.;
        nodeMat(1,6) = +1.;
        nodeMat(1,7) =  0.;
                break;
            }
    }

    POINT contralPoint0(2);
    POINT contralPoint1(2);
    contralPoint0[0] = -1.;
    contralPoint0[1] = -1.;
    contralPoint1[0] = +1.;
    contralPoint1[1] = +1.;

    MATRIX quadPoints;
        /*
         * use equivalent parameter commutation
         */
    GaussInte.GQPoints(quadPoints, \
                contralPoint0, \
                contralPoint1);

    GaussInte.GQWeight(quadParameters, \
                contralPoint0, \
                contralPoint1);

    quadParameters *= 4.;

    MATRIX matJacob(2,NumShape);
    MATRIX matShape(2,2*NumShape);
    for(int ip=0;ip<quadPoints.cols();ip++){
        basefunction(matShape,matJacob,quadPoints.col(ip));
        ShapeMats.push_back(matShape);
        JacobMats.push_back(matJacob);
    }

    MATRIX matJacobN(2,NumShape);
    MATRIX matShapeN(2,2*NumShape);
    for(int ip=0;ip<nodeMat.cols();ip++){
        basefunction(matShapeN,matJacobN,nodeMat.col(ip));
        ShapeMatsNode.push_back(matShapeN);
        JacobMatsNode.push_back(matJacobN);
    }
}

bool SolverElesticity::solve(){

    /// Local Stiffness Matrix
    LocalStiffness();
    /// Local Load
    LocalLoad();
    /// Global Stiffness and Load; 
    GeneralLinearEquation();
    /// Neumann Boundary condition
    /// Neuman boundary condition must be applied at first
    NeumanBound();
    /// Dirichlet boundary condition;
    DirichBound();
    /* Solve the linear equations;
     * Applying the SVD decomposition within Eigen3
     *  for linear least squares
     */
    //VEC displacementX=DisplaceStiffness.bdcSvd(ComputeThinU | ComputeThinV).solve(loadDiri);
    VEC displacementX=DisplaceStiffness.ldlt().solve(loadDiri);

    /// get the displacement
    INTE iterDirich = 0;
    INTE iterDisp = 0;
    for (int iter=0;iter<displacement.size();iter++){
        if(ifDirichlet[iter]){
            displacement[iter] = valDirichlet[iterDirich];
            iterDirich++;
        }else{
            displacement[iter] = displacementX[iterDisp];
            iterDisp++;
        }
    }
    return true;
};

/*
 * set the Dirichlet boundary conditions
 */
void SolverElesticity::DirichBound(){

    MatrixXd DirichletStiffness0(DISN,idDirichlet.size());
    MatrixXd DisplaceStiffness0(DISN,DISN-idDirichlet.size());
    DirichletStiffness.resize(DISN-idDirichlet.size(),idDirichlet.size());
    DisplaceStiffness.resize(DISN-idDirichlet.size(),DISN-idDirichlet.size());

    INTE iterDirich = 0;
    INTE iterDisplace = 0;
    for(int iter = 0;iter<displacement.size();iter++){
        if(ifDirichlet[iter]){
            DirichletStiffness0.col(iterDirich) = stiffness.col(iter);
            iterDirich++;
        }else{
            loadDiri[iterDisplace] = load[iter];
            DisplaceStiffness0.col(iterDisplace) = stiffness.col(iter);
            iterDisplace++;
        }
    }
    iterDirich = 0;
    iterDisplace = 0;
    for(int iter = 0;iter<displacement.size();iter++){
        if(ifDirichlet[iter]){
            DirichletStiffness.row(iterDirich) = DirichletStiffness0.row(iter);
            iterDirich++;
        }else{
            DisplaceStiffness.row(iterDisplace) = DisplaceStiffness0.row(iter);
            iterDisplace++;
        }
    }
    loadDiri += - DirichletStiffness*valDirichlet;
};

/*
 * Neuman Bopundary conditions;
 * the right hand of the equation
 */
void SolverElesticity::NeumanBound(){

    for(int iter=0;iter<valNeuman.size();iter++){
        int Neuman=idNeuman[iter];
        load[Neuman] += valNeuman[iter];
    }
};

/*
 *  get the global information through local information
 */
void SolverElesticity::GeneralLinearEquation(){
    load.resize(DISN);
    for(int li=0;li<DISN;li++){
        load[li] = 0.;
    }
    for(ULINTE iterEle=0;iterEle<mesh.elements.size();iterEle++){
        vector<INTE> itGlobal;
        for(ULINTE itLocal=0;itLocal<mesh.elements[iterEle].vertex.size();itLocal++){
            itGlobal.push_back(mesh.elements[iterEle].vertex[itLocal]);
        }

        vector<INTE> itGlobal2(2*itGlobal.size());
        for(ULINTE itrow=0;itrow<itGlobal.size();itrow++){
            itGlobal2[2*itrow]=2*itGlobal[itrow];
            itGlobal2[2*itrow+1]=2*itGlobal[itrow]+1;
        }
        
        for(ULINTE itrow=0;itrow<itGlobal2.size();itrow++){
            for(ULINTE itcol=0;itcol<itGlobal2.size();itcol++){
                stiffness(itGlobal2[itrow],itGlobal2[itcol]) \
                    += stiffnessLocal[iterEle](itrow,itcol);
            }
            load[itGlobal2[itrow]] += loadLocal[iterEle][itrow];
        }
    }
};

/*
 * Calculate local stiffness matrix and the load vector
 */
void SolverElesticity::LocalStiffness(){
    set_ShapeFunction();
    int NS = 2*NumShape;
    for(vector<Element>::iterator it=mesh.elements.begin();it<mesh.elements.end();it++){
        MatrixXd stiffi(NS,NS);
        VEC loadi(NS);
        for(int ri=0;ri<stiffi.rows();ri++){
            for(int ci=0;ci<stiffi.cols();ci++){
                stiffi(ri,ci) = 0.;
            }
            loadi(ri) = 0.;
        }
        ///
        MATRIX MB(3,NS);
        for(int ri=0;ri<MB.rows();ri++){
            for(int ci=0;ci<MB.cols();ci++){
                MB(ri,ci) =0.;
            }
        }
        MATRIX MK(NS,NS);
        for(int ri=0;ri<MK.rows();ri++){
            for(int ci=0;ci<MK.cols();ci++){
                MK(ri,ci) =0.;
            }
        }
        MATRIX shape;
        MATRIX JacoMat(Dim,Dim);
        for(int ip=0;ip<quadParameters.size();ip++){

            JacoMat = JacobMats[ip] * mesh.contralPoints[(*it).id];
            shape = JacoMat.inverse() * JacobMats[ip];
            for(int ns=0;ns<NumShape;ns++){
                MB(0,2*ns  ) = shape(0,ns);
                MB(1,2*ns+1) = shape(1,ns);
                MB(2,2*ns  ) = shape(1,ns);
                MB(2,2*ns+1) = shape(0,ns);
            }
            MATRIX stress = MD*MB;
            MK = MB.transpose()*stress;

            VEC fvec(2);
            fvec[0] = 0.;
            fvec[1] = 0.;
            VEC VP = ShapeMats[ip].transpose()*fvec;

            double lambda = JacoMat.determinant();
            stiffi += quadParameters[ip]*lambda*MK;
            loadi += quadParameters[ip]*lambda*VP;
            StrainMats.push_back(MB);
            StressMats.push_back(stress);
        }
        stiffnessLocal.push_back(stiffi);
        loadLocal.push_back(loadi);
    }
};
/*
 * local load
 */
void SolverElesticity::LocalLoad(){
/// the local load calculation is included in the LocalStiffness;
};

void SolverElesticity::set_DirichBound(){

    for(int ri = 0;ri<NY+1;ri++){
        idDirichlet.push_back(2*((NX+1)*ri));
        idDirichlet.push_back(2*((NX+1)*ri)+1);
    }
    valDirichlet.resize(idDirichlet.size());
    for(int ni=0;ni<displacement.size();ni++){
        ifDirichlet.push_back(false);
    }
    for(unsigned long ni=0;ni<idDirichlet.size();ni++){
        ifDirichlet[idDirichlet[ni]] = true;
        valDirichlet[ni] = 0.;
    }

    loadDiri.resize(DISN-idDirichlet.size());
    for(int li=0;li<loadDiri.size();li++){
        loadDiri[li] = 0.;
    }
}

void SolverElesticity::set_NeumanBound(){
    for (int ri = 0;ri<NY+1;ri++){
        idNeuman.push_back(2*((NX+1)*(ri+1)-1));
    }
    valNeuman.resize(idNeuman.size());

    for(int ni=0;ni<displacement.size();ni++){
        ifNeumann.push_back(false);
    }
    for(unsigned long ni=0;ni<idNeuman.size();ni++){
        ifNeumann[idNeuman[ni]] = true;
        valNeuman[ni] = 1.;
    }
}

void SolverElesticity::basefunction(MATRIX& matShape, MATRIX& matJacob, VEC xy){

    VEC x,y;
    x.resize(2);
    y.resize(2);

    x[0] = -1.;
    x[1] = +1.;
    y[0] = -1.;
    y[1] = +1.;

        /*
         * use equivalent parameter commutation
         */
    matJacob.resize(2,NumShape);
    matShape.resize(2,2*NumShape);
    for(int ri=0;ri<2;ri++){
        for(int ci=0;ci<NumShape;ci++){
            matJacob(ri,ci) = 0.;
        }
    }
    for(int ri=0;ri<2;ri++){
        for(int ci=0;ci<2*NumShape;ci++){
            matShape(ri,ci) =0.;
        }
    }

        double N1,N2,N3,N4;

        N1 = 0.25* (1. + xy[0]*x[0])*(1. + xy[1]*y[0]);
        N2 = 0.25* (1. + xy[0]*x[1])*(1. + xy[1]*y[0]);
        N3 = 0.25* (1. + xy[0]*x[1])*(1. + xy[1]*y[1]);
        N4 = 0.25* (1. + xy[0]*x[0])*(1. + xy[1]*y[1]);

        switch(NumShape){
        case 4 : {
            matShape(0,0) = N1;
            matShape(0,2) = N2;
            matShape(0,4) = N3;
            matShape(0,6) = N4;
            matShape(1,1) = matShape(0,0);
            matShape(1,3) = matShape(0,2);
            matShape(1,5) = matShape(0,4);
            matShape(1,7) = matShape(0,6);

            matJacob(0,0) = 0.25* (x[0])*(1. + xy[1]*y[0]);
            matJacob(0,1) = 0.25* (x[1])*(1. + xy[1]*y[0]);
            matJacob(0,2) = 0.25* (x[1])*(1. + xy[1]*y[1]);
            matJacob(0,3) = 0.25* (x[0])*(1. + xy[1]*y[1]);
            matJacob(1,0) = 0.25* (y[0])*(1. + xy[0]*x[0]);
            matJacob(1,1) = 0.25* (y[0])*(1. + xy[0]*x[1]);
            matJacob(1,2) = 0.25* (y[1])*(1. + xy[0]*x[1]);
            matJacob(1,3) = 0.25* (y[1])*(1. + xy[0]*x[0]);

        break;
                 }
        case 8 : {
            matShape(0, 8) = 0.5 * (1. - xy[0]*xy[0])*(1. - xy[1]);
            matShape(0,10) = 0.5 * (1. - xy[1]*xy[1])*(1. + xy[0]);
            matShape(0,12) = 0.5 * (1. - xy[0]*xy[0])*(1. + xy[1]);
            matShape(0,14) = 0.5 * (1. - xy[1]*xy[1])*(1. - xy[0]);

            matShape(0,0) = N1 - 0.5*matShape(0, 8) -0.5*matShape(0,14);
            matShape(0,2) = N2 - 0.5*matShape(0, 8) -0.5*matShape(0,10);
            matShape(0,4) = N3 - 0.5*matShape(0,10) -0.5*matShape(0,12);
            matShape(0,6) = N4 - 0.5*matShape(0,12) -0.5*matShape(0,14);

            matShape(1,1) = matShape(0,0);
            matShape(1,3) = matShape(0,2);
            matShape(1,5) = matShape(0,4);
            matShape(1,7) = matShape(0,6);
            matShape(1,9) = matShape(0,8);
            matShape(1,11) = matShape(0,10);
            matShape(1,13) = matShape(0,12);
            matShape(1,15) = matShape(0,14);

            matJacob(0,4) = - xy[0]*(1. - xy[1]);
            matJacob(1,4) = - 0.5 * (1. - xy[0]*xy[0]);
            matJacob(0,5) =   0.5 * (1. - xy[1]*xy[1]);
            matJacob(1,5) = - xy[1]*(1. + xy[0]);
            matJacob(0,6) = - xy[0]*(1. + xy[1]);
            matJacob(1,6) =   0.5 * (1. - xy[0]*xy[0]);
            matJacob(0,7) = - 0.5 * (1. - xy[1]*xy[1]);
            matJacob(1,7) = - xy[1]*(1. - xy[0]);

            matJacob(0,0) = 0.25* (x[0])*(1. + xy[1]*y[0]) \
                            -0.5* matJacob(0,4) \
                            -0.5* matJacob(0,7);

            matJacob(1,0) = 0.25* (y[0])*(1. + xy[0]*x[0]) \
                            -0.5* matJacob(1,4) \
                            -0.5* matJacob(1,7) ;

            matJacob(0,1) = 0.25* (x[1])*(1. + xy[1]*y[0]) \
                            -0.5* matJacob(0,4) \
                            -0.5* matJacob(0,5) ;

            matJacob(1,1) = 0.25* (y[0])*(1. + xy[0]*x[1]) \
                            -0.5* matJacob(1,4) \
                            -0.5* matJacob(1,5);

            matJacob(0,2) = 0.25* (x[1])*(1. + xy[1]*y[1]) \
                            -0.5* matJacob(0,5) \
                            -0.5* matJacob(0,6) ;

            matJacob(1,2) = 0.25* (y[1])*(1. + xy[0]*x[1]) \
                            -0.5* matJacob(1,5) \
                            -0.5* matJacob(1,6) ;

            matJacob(0,3) = 0.25* (x[0])*(1. + xy[1]*y[1]) \
                            -0.5* matJacob(0,6) \
                            -0.5* matJacob(0,7) ;

            matJacob(1,3) = 0.25* (y[1])*(1. + xy[0]*x[0]) \
                            -0.5* matJacob(1,6) \
                            -0.5* matJacob(1,7) ;

                     break;
                 }
        default :{
                     cout << "Error 1 : wrong parameter; NumShape "<<endl;
                     break;
                 }
        }
}

bool SolverElesticity::readFile(string ){
    return true;
}
