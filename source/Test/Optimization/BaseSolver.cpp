#include "BaseSolver.h"

/// get the strain on quadracture points
void BaseFEM::get_Strain_Quad(MATRIX& strain){
    int NS = 2*NumShape;
    ifStrain = true;
    strain.resize(3,mesh.elements.size()*quadParameters.size());
    int iter=0;
    for(vector<Element>::iterator ei=mesh.elements.begin();ei<mesh.elements.end();ei++){
        VEC vertex(NS);
        for(int vi=0;vi<NumShape;vi++){
            vertex(2*vi) = displacement[2*(*ei).vertex[vi]];
            vertex(2*vi+1) = displacement[2*(*ei).vertex[vi]+1];
        }
        for(int ip=0;ip<quadParameters.size();ip++){
            strain.col(iter) = StrainMats[quadParameters.size()*(*ei).id+ip]*vertex;
            iter++;
        }
    };
    /// out to file
    ofstream outfile("strain_quad.dat");
    for(int ci = 0;ci<strain.cols();ci++){
        for(int ri =0;ri<strain.rows();ri++){
            outfile << strain(ri,ci) << '\t' ;
        }
        outfile << endl;
    }
    outfile.close();
}
/// get the strain on nodes
void BaseFEM::get_Strain_Node(MATRIX& strain){
    int NS = 2*NumShape;
    ifStrain = true;
    strain.resize(3,mesh.elements.size()*NumShape);
    int iter=0;
    for(vector<Element>::iterator ei=mesh.elements.begin();ei<mesh.elements.end();ei++){
        VEC vertex(NS);
        for(int vi=0;vi<NumShape;vi++){
            vertex(2*vi) = displacement[2*(*ei).vertex[vi]];
            vertex(2*vi+1) = displacement[2*(*ei).vertex[vi]+1];
        }
        for(int ip=0;ip<NumShape;ip++){
            strain.col(iter) = StrainMatsNode[NumShape*(*ei).id+ip]*vertex;
            iter++;
        }
    };
    /// out to file
    ofstream outfile("strain_quad.dat");
    for(int ci = 0;ci<strain.cols();ci++){
        for(int ri =0;ri<strain.rows();ri++){
            outfile << strain(ri,ci) << '\t' ;
        }
        outfile << endl;
    }
    outfile.close();
}
/// get the stress on quadracture points
void BaseFEM::get_Stress_Quad(MATRIX& stress){
    ifStress = true;
    stress.resize(3,mesh.elements.size()*quadParameters.size());
    int iter=0;
    for(vector<Element>::iterator ei=mesh.elements.begin();ei<mesh.elements.end();ei++){
        VEC vertex(2*NumShape);
        for(int vi=0;vi<NumShape;vi++){
            vertex(2*vi) = displacement[2*(*ei).vertex[vi]];
            vertex(2*vi+1) = displacement[2*(*ei).vertex[vi]+1];
        }
        for(int ip=0;ip<quadParameters.size();ip++){
            stress.col(iter) = StressMats[quadParameters.size()*(*ei).id+ip]*vertex;
            iter++;
        }
    };
    /// out to file
    ofstream outfile("stress_quad.dat");
    for(int ci = 0;ci<stress.cols();ci++){
        for(int ri =0;ri<stress.rows();ri++){
            outfile<< stress(ri,ci) << '\t' ;
        }
        outfile << endl;
    }
    outfile.close();
}
/// get the stress on Nodes
void BaseFEM::get_Stress_Node(MATRIX& stress){
    ifStress = true;
    stress.resize(3,mesh.elements.size()*NumShape);
    int iter=0;
    for(vector<Element>::iterator ei=mesh.elements.begin();ei<mesh.elements.end();ei++){
        VEC vertex(2*NumShape);
        for(int vi=0;vi<NumShape;vi++){
            vertex(2*vi) = displacement[2*(*ei).vertex[vi]];
            vertex(2*vi+1) = displacement[2*(*ei).vertex[vi]+1];
        }
        for(int ip=0;ip<NumShape;ip++){
            stress.col(iter) = StressMatsNode[NumShape*(*ei).id+ip]*vertex;
            iter++;
        }
    };
    /// out to file
    ofstream outfile("stress_node.dat");
    for(int ci = 0;ci<stress.cols();ci++){
        for(int ri =0;ri<stress.rows();ri++){
            outfile<< stress(ri,ci) << '\t' ;
        }
        outfile << endl;
    }
    outfile.close();
}
/// get the Matrix information on each Nodes
void BaseFEM::get_NodeMatInfo(){
    int NS = 2*NumShape;
    for(vector<Element>::iterator it=mesh.elements.begin();it<mesh.elements.end();it++){
        /// strain Matrix
        MATRIX MB(3,NS);
        for(int ri=0;ri<MB.rows();ri++){
            for(int ci=0;ci<MB.cols();ci++){
                MB(ri,ci) =0.;
            }
        }
        /// stiffness Matrix
        MATRIX MK(NS,NS);
        for(int ri=0;ri<MK.rows();ri++){
            for(int ci=0;ci<MK.cols();ci++){
                MK(ri,ci) =0.;
            }
        }
        /// for construct the strain Mat
        MATRIX shape;
        MATRIX JacoMat(Dim,Dim);
        for(int ip=0;ip<NumShape;ip++){
            JacoMat = JacobMatsNode[ip] * mesh.contralPoints[(*it).id];
            shape = JacoMat.inverse() * JacobMatsNode[ip];
            for(int ns=0;ns<NumShape;ns++){
                MB(0,2*ns  ) = shape(0,ns);
                MB(1,2*ns+1) = shape(1,ns);
                MB(2,2*ns  ) = shape(1,ns);
                MB(2,2*ns+1) = shape(0,ns);
            }
            /// stress Matrix
            MATRIX stress = MD*MB;
            StrainMatsNode.push_back(MB);
            StressMatsNode.push_back(stress);
        }
    }
}
/// write result into file
bool BaseFEM::writeFile(string filename) const {
    ofstream logFile(filename);
    logFile << "displacement of each nodes:" <<endl;
    for (ULINTE ni=0;ni<mesh.points.size();ni++){
        logFile << setw(12) <<dec<< ni;
        logFile << setw(12) <<setprecision(3) << displacement[2*ni];
        logFile << setw(12) <<setprecision(3) << displacement[2*ni+1];
        logFile << endl;
    }
    logFile.close();
    return true;
}
