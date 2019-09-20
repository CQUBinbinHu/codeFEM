#include <sstream>
#include "PreProcess.h"
using namespace std;
using namespace Optimization;
typedef ULINTE U_LONG;
/// get the Mesh form file of the outputs of ABAQUS 
bool GenerateMesh::get_MeshFromFile(string filename,Mesh& mesh){
    /*
     * Read from file for meshgrid;
     */
    vector<POINT2> nodes;
    vector<VectorXi> elements;
    POINT2 coord;
    VectorXi ele;
    ele.resize(1+NVert);
    string type;
    string newline;
    string typeNode("*Node\r");
    string typeEle("*Element\r");
    string typeEnd("*End\r");
    string typeNode1("*Node");
    string typeEle1("*Element");
    string typeEnd1("*End");

    char GET_CHAR;
    int GET_INT;

    ifstream infile(filename);
    std::getline(infile,type);
    if (type == typeNode || type==typeNode1){
        while(infile.eof()==false){
            std::getline(infile,newline);
            if(newline == typeEle || newline==typeEle1){
                break;
            }else{
                stringstream lineStream(newline);
                lineStream >> GET_INT;
                lineStream >> GET_CHAR;
                lineStream >> coord[0];
                lineStream >> GET_CHAR;
                lineStream >> coord[1];
                nodes.push_back(coord);
            }
        }
    }
    if(newline == typeEle||newline==typeEle1){
        while (infile.eof()==false){
            getline(infile,newline);
            if(newline == typeEnd || newline == typeEnd1){
                break;
            }else{
                stringstream lineStream(newline);
                for(int ni=0;ni<ele.size();ni++){
                    lineStream >> ele[ni];
                    lineStream >> GET_CHAR;
                }
                elements.push_back(ele);
            }
        }
    }
    /*
     * Get the mesh from file;
     */
    //int GN = 2;
    //NumerIntegral GaussInte(GN,'L');
    //
    mesh.points = nodes;

    for(U_LONG ei=0;ei<elements.size();ei++){
        Element eleSample;
        eleSample.id = elements[ei][0]-1;
        eleSample.vertex.resize(NVert);
        for (int vi=0;vi<NVert;vi++){
            eleSample.vertex[vi] = elements[ei][vi+1]-1;
        };
        MATRIX contralPoints(2,NVert);
        for(int vi=0;vi<NVert;vi++){
            contralPoints.col(vi) = mesh.points[eleSample.vertex[vi]];
        }
        /*
         * use equivalent parameter commutation
         */
        mesh.elements.push_back(eleSample);
        mesh.contralPoints.push_back(contralPoints.transpose());
    }
    infile.close();
    return true;
};
/*
 * set the Geometry of the domain
 */
bool GenerateMesh::set_GeometryFromFile(string filename){
    return true;
};
/*
 *  Generate uniform quadracture element mesh
 */
inline int indexPoint(int xi,int yi,int Nx){return xi+yi*Nx+yi;};
bool GenerateMesh::set_Geometry(vector<POINT> geo,vector<int> Nxy) {
   return true;
}
bool GenerateMesh::get_Mesh(Mesh& mesh){
    return true;
};