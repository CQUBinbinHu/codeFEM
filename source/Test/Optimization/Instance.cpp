#include "Instance.h"

bool INSTANCE_Man::set_BoundFromFile(string filename){
    ifstream infile(filename);

    infile.close();
    return true;
};

void INSTANCE_Man::set_DirichBound(){

    vector<int> vecDrich;
    vecDrich.push_back(0);
    vecDrich.push_back(3);
    vecDrich.push_back(33);
    vecDrich.push_back(34);
    vecDrich.push_back(35);
    vecDrich.push_back(36);
    vecDrich.push_back(37);
    vecDrich.push_back(38);
    vecDrich.push_back(39);
    vecDrich.push_back(40);
    vecDrich.push_back(41);
    for(int ri = 0;ri<vecDrich.size();ri++){
        idDirichlet.push_back(2*vecDrich[ri]);
        idDirichlet.push_back(2*vecDrich[ri]+1);
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

};

void INSTANCE_Man::set_NeumanBound(){

    vector<int> vecNeuman;
    vecNeuman.push_back(1);
    vecNeuman.push_back(2);
    vecNeuman.push_back(15);
    vecNeuman.push_back(16);
    vecNeuman.push_back(17);
    vecNeuman.push_back(18);
    vecNeuman.push_back(19);
    vecNeuman.push_back(20);
    vecNeuman.push_back(21);
    vecNeuman.push_back(22);
    vecNeuman.push_back(23);

    for (int ri = 0;ri<vecNeuman.size();ri++){
        idNeuman.push_back(2*vecNeuman[ri]);
    }
    valNeuman.resize(idNeuman.size());

    for(int ni=0;ni<displacement.size();ni++){
        ifNeumann.push_back(false);
    }
    for(unsigned long ni=0;ni<idNeuman.size();ni++){
        ifNeumann[idNeuman[ni]] = true;
        valNeuman[ni] = 1.;
    }
    valNeuman[0] = 0.5;
    valNeuman[idNeuman.size()-1] = 0.5;

};
////////////////////////////////////////////////////////////////
void INSTANCE2::set_DirichBound(){
    for(int ri = 0;ri<NY+1;ri++){
        idDirichlet.push_back(2*((NX+1)*ri));
        //idDirichlet.push_back(2*((mesh.NX+1)*ri)+1);
    }
    for(int ri = 0;ri<NX+1;ri++){
        idDirichlet.push_back(2*ri+1);
        //idDirichlet.push_back(2*((mesh.NX+1)*ri)+1);
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
};

void INSTANCE2::set_NeumanBound(){

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
    valNeuman[0] = 0.5;
    valNeuman[idNeuman.size()-1] = 0.5;
};
