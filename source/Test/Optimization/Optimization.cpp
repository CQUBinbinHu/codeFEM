#include "Optimization.h"
using namespace Optimization;

bool SolverOptimization::Initialization(){

    /// Finite Element Method
    FEM.set_NumShape(PN);
    FEM.set_DirichBound();
    FEM.set_NeumanBound();
    /// get QuadPoints
    int GN = 3;
    NumerIntegral GaussInte(GN,'L');
    POINT contralPoint0(2);
    POINT contralPoint1(2);
    contralPoint0[0] = -1.;
    contralPoint0[1] = -1.;
    contralPoint1[0] = +1.;
    contralPoint1[1] = +1.;

        /*
         * use equivalent parameter commutation
         */
    GaussInte.GQPoints(QuadPoints, \
                contralPoint0, \
                contralPoint1);

    GaussInte.GQWeight(quadParameters, \
                contralPoint0, \
                contralPoint1);

    quadParameters *= 4.;
    
    /// get ShapeFunctionsQuad
    /// get Grad1ShapeFunctionsQuad
    /// get Grad2ShapeFunctionsQuad
    for(int qi = 0;qi<quadParameters.size();qi++){
        VEC shapeInfo;
        getShapeFunction(shapeInfo,QuadPoints.col(qi));
        ShapeFunctionsQuad.push_back(shapeInfo);
        MATRIX GradShapeInfo;
        getShapeFunctionGrad(GradShapeInfo,QuadPoints.col(qi));
        Grad1ShapeFunctionQuad.push_back(shapeInfo);
        getShapeFunctionGrad2(GradShapeInfo,QuadPoints.col(qi));
        Grad2ShapeFunctionQuad.push_back(shapeInfo);
    }
    return true;
};
void SolverOptimization::getDescent(){

    MATRIX matDescent(EN*quadParameters.size(),level.size());
    VEC levelDescentQuad(EN*quadParameters.size());

    int iter = 0;
    MATRIX strain,stress;
    FEM.get_Strain_Quad(strain);
    FEM.get_Stress_Quad(stress);
    VEC levelVec(PN);
    for(vector<Element>::iterator it=mesh.elements.begin();it<mesh.elements.end();it++){
        for(int vi=0;vi<PN;vi++){
            levelVec[vi] = level[(*it).vertex[vi]];
        }
        MATRIX UM(PN,2);
        for(int ip=0;ip<PN;ip++){
            UM(ip,0) = displacement[2*(*it).vertex[ip]];
            UM(ip,1) = displacement[2*(*it).vertex[ip]+1];
        }
        for(int ip=0;ip<quadParameters.size();ip++){

            matDescent(iter,0) = 0.;
            for(int ci=0;ci<level.size();ci++){
                matDescent(iter,ci) = 0.;
            }
            for(U_LONG vi=0;vi<(*it).vertex.size();vi++){
                matDescent(iter,(*it).vertex[vi]) = ShapeFunctionsQuad[ip][vi];
            }
            /// get descent level
            MATRIX Jacob = Grad1ShapeFunctionQuad[ip] * mesh.contralPoints[(*it).id];
            
            double LDQ = 0.;
            LDQ -= 100.;
            double En = (stress.col(iter).transpose()*strain.col(iter));
            LDQ += En * Jacob.determinant();

            VEC Q1(2),Q2(2),Q3(2);
            /// get Gradient of Phi
            VEC GradPhi(2);
            GradPhi = Jacob.reverse() * Grad1ShapeFunctionQuad[ip] * levelVec;
            /// get Gradient2 of Phi
            MATRIX GradPhi2Xi(2,2);
            GradPhi2Xi(0,0) = 0.;
            for(int vi=0;vi<PN;vi++){
                GradPhi2Xi(0,0) += Grad2ShapeFunctionQuad[ip](0,ip)*levelVec[ip];
                GradPhi2Xi(0,1) += Grad2ShapeFunctionQuad[ip](2,ip)*levelVec[ip];
                GradPhi2Xi(1,0) += Grad2ShapeFunctionQuad[ip](2,ip)*levelVec[ip];
                GradPhi2Xi(1,1) += Grad2ShapeFunctionQuad[ip](1,ip)*levelVec[ip];
            }
            ///
            MATRIX JacobMat2(3,2);
            JacobMat2 = Grad2ShapeFunctionQuad[ip]*mesh.contralPoints[(*it).id];
            MATRIX deltGradPhi(2,2);
            deltGradPhi(0,0) = JacobMat2(0,0)*GradPhi[0] + JacobMat2(0,1)*GradPhi[1];
            deltGradPhi(1,1) = JacobMat2(1,0)*GradPhi[0] + JacobMat2(1,1)*GradPhi[1];
            deltGradPhi(0,1) = JacobMat2(2,0)*GradPhi[0] + JacobMat2(2,1)*GradPhi[1];
            deltGradPhi(1,0) = JacobMat2(2,0)*GradPhi[0] + JacobMat2(2,1)*GradPhi[1];
            MATRIX GradPhi2(2,2);
            GradPhi2 = Jacob.reverse() * (GradPhi2 - deltGradPhi) \
                     * Jacob.reverse().transpose();
            /// stress Matrix
            MATRIX ST(2,2);
            ST(0,0) = stress(0,iter);
            ST(1,1) = stress(1,iter);
            ST(0,1) = stress(2,iter);
            ST(1,0) = stress(2,iter);
            /// gradient of U
            MATRIX DU(2,2);
            DU = Grad1ShapeFunctionQuad[ip] * UM ;
            /// get Q1
            Q1 = -2. * DU * ST * GradPhi;
            Q3 = -2. * ST * GradPhi2 * UM.transpose() \
                 * ShapeFunctionsQuad[ip].transpose();
            ///////////////////////////////////////////////////////////////////
            //
            //  Q2 : = -2. * u . Ae(\nabla u) . \nabla Phi
            //
            ///////////////////////////////////////////////////////////////////

            double GradPhiNorm = std::sqrt(GradPhi.transpose()*GradPhi);
            double L0 = (GradPhi.transpose()*(Q1+Q2+Q3));
            LDQ -= (1./GradPhiNorm)* L0;
            levelDescentQuad(iter) = LDQ * Jacob.determinant();
            /// need to accomplish
            iter++;
        }
    } 
    /// get levelDescent using linear least squares
    levelDescent = matDescent.bdcSvd(ComputeThinU | ComputeThinV).solve(levelDescentQuad);
};
double SolverOptimization::getObjection(VEC& newLevel){

    FEM.solve();
    double obj = 0.;
    VEC levelVec(PN);
    for(vector<Element>::iterator it=mesh.elements.begin();it<mesh.elements.end();it++){
        for(int vi=0;vi<PN;vi++){
            levelVec[vi] = newLevel[(*it).vertex[vi]];
        }
        for(int ip=0;ip<quadParameters.size();ip++){
            MATRIX JacoMat = Grad1ShapeFunctionQuad[ip] * mesh.contralPoints[(*it).id];
            double levelPo = ShapeFunctionsQuad[ip].transpose()*levelVec; 
            obj += HeasiveFun(levelPo) * JacoMat.determinant();
        }
    }
    obj *= 100.;
    for(int ni=0;ni<FEM.get_valNeuman().size();ni++){
        obj += FEM.get_valNeuman()[ni]*displacement[FEM.get_idNeuman()[ni]];
    }

    return obj;
};
void SolverOptimization::getShapeFunction(VEC& shapeInfo,VEC point){

    shapeInfo.resize(PN);
    switch(PN){
    case 4 : {
             MATRIX xis(4,2);
             xis(0,0) = -1.;
             xis(1,0) = +1.;
             xis(2,0) = +1.;
             xis(3,0) = -1.;
             xis(0,1) = -1.;
             xis(1,1) = -1.;
             xis(2,1) = +1.;
             xis(3,1) = +1.;
             for(int ip=0;ip<4;ip++){
                shapeInfo[ip] = 0.25*(1.+xis(ip,0)*point[0])*(1.+xis(ip,1)*point[1]);
             }
                break;
             }
    case 8 : {
                 break;
             }
    default :{
                 break;
             }
    };
};
void SolverOptimization::getShapeFunctionGrad(MATRIX& shapeInfo,VEC point){

    shapeInfo.resize(2,PN);
    switch(PN){
    case 4 : {
             MATRIX xis(4,2);
             xis(0,0) = -1.;
             xis(1,0) = +1.;
             xis(2,0) = +1.;
             xis(3,0) = -1.;
             xis(0,1) = -1.;
             xis(1,1) = -1.;
             xis(2,1) = +1.;
             xis(3,1) = +1.;
             for(int ip=0;ip<4;ip++){
                shapeInfo(0,ip) = 0.25*xis(ip,0)*(1.+xis(ip,1)*point[1]);
                shapeInfo(1,ip) = 0.25*xis(ip,1)*(1.+xis(ip,0)*point[0]);
             }
                break;
             }
    case 8 : {
                 break;
             }
    default :{
                 break;
             }
    };
};
void SolverOptimization::getShapeFunctionGrad2(MATRIX& shapeInfo,VEC point){

    shapeInfo.resize(2,PN);
    switch(PN){
    case 4 : {
             MATRIX xis(4,2);
             xis(0,0) = -1.;
             xis(1,0) = +1.;
             xis(2,0) = +1.;
             xis(3,0) = -1.;
             xis(0,1) = -1.;
             xis(1,1) = -1.;
             xis(2,1) = +1.;
             xis(3,1) = +1.;
             for(int ip=0;ip<4;ip++){
                shapeInfo(0,ip) = 0.;
                shapeInfo(1,ip) = 0.;
                shapeInfo(2,ip) = 0.25*xis(ip,0)*xis(ip,1);
             }
                break;
             };
    case 8 : {
                 break;
             };
    default :{
                 break;
             };
    };
};
void SolverOptimization::SearchLine(double& delta){
    double Ros = 0.1;
    vector<double> searchRo;
    double valObject0 = 0.;
    double valObject1 = 0.;
    double valObject2 = 0.;

    double r1,r2;
    double ra,rb;
    // finite element
    valObject0 = getObjection(level);
    valObject1 = valObject0;
    getDescent();

    searchRo.push_back(Ros);
    VEC levelNew = level + Ros * levelDescent;
    // finite element
    valObject2 = getObjection(levelNew);
    /// update Ro
    if(valObject2 < valObject1){
        while (valObject2 < valObject1 ){
            Ros *= 2.;
            searchRo.push_back(Ros);
            levelNew = level + Ros * levelDescent;
            valObject1 = valObject2;
            valObject2 = getObjection(levelNew);
        }
        vector<double>::iterator ro = searchRo.end();
        ra = (*(ro--));
        rb = (*ro);
    }else{
        ra = 0;
        rb = Ros; 
    }
    /// employing 1-dim search using 0.618 method
    while(abs(ra-rb)>eps){
        r1 = ra + 0.382 * (rb-ra);
        r2 = ra + 0.618 * (rb-ra);
        levelNew = level + ra * levelDescent;
        valObject1 = getObjection(levelNew);
        levelNew = level + rb * levelDescent;
        valObject2 = getObjection(levelNew);
        if(valObject1>valObject2){
            ra = r1;
        }else{
            rb = r2;
        }
    }
    level += ra * levelDescent;
    delta = (valObject0 - valObject1)/valObject0;
};
bool SolverOptimization::solve(){

    double delta = 1.;
    while (delta>eps){
        /// update
        getDescent();
        SearchLine(delta);
    }
    return true;
};

double SolverOptimization::HeasiveFun(double){
    return 0.;
};