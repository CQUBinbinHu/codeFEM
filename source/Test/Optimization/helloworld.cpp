#include <iostream>
#include "../../Common/math.h"
#include <ctime>

using namespace std;

int main(){

    // helloworld
    cout<< "code helloworld!"<<endl;
    // for runtime checkout
    clock_t startTime,endTime;
    // initializaiton of Matrix
    ULINTE SN = 10000;
    vector<vector<REAL>> vec1;
    vector<vector<REAL>> vec2;
    vector<vector<REAL>> vec3;
    // use matrix in Eigen
    MATRIX mat1(SN,SN);
    MATRIX mat2(SN,SN);
    MATRIX mat3(SN,SN);
    for(int mi=0;mi<SN;mi++){
        for(int ni=0;ni<SN;ni++){
            mat1(mi,ni) = 0.;
            mat2(mi,ni) = 0.;
            mat3(mi,ni) = 0.;
        }
    }
    for(ULINTE vi=0;vi<SN;vi++){
        vector<REAL> vec0(SN);
        vec1.push_back(vec0);
        vec2.push_back(vec0);
        vec3.push_back(vec0);
    }
    // matrix operation #1
    startTime = clock();
    for(ULINTE mi=0;mi<SN;mi++){
        for(ULINTE ni = 0; ni<SN; ni++){
            vec3[mi][ni] = 0.;
            for(ULINTE vi=0;vi<SN;vi++){
                vec3[mi][ni] += vec1[mi][vi]*vec2[vi][ni];
            }
        }
    }
    endTime = clock();
    cout<< "runtime of operation #1 : " \
        <<(double)(endTime-startTime)/CLOCKS_PER_SEC<< "s"<<endl;
    // matrix operation #2
    startTime = clock();
    for(ULINTE mi=0;mi<SN;mi++){
        for(ULINTE ni = 0; ni<SN; ni++){
            vec3[ni][mi] = 0.;
            for(ULINTE vi=0;vi<SN;vi++){
                vec3[ni][mi] += vec1[ni][vi]*vec2[vi][mi];
            }
        }
    }
    endTime = clock();
    cout<< "runtime of operation #2 : " \
        <<(double)(endTime-startTime)/CLOCKS_PER_SEC<< "s"<<endl;
    // matrix operation #3
    startTime = clock();
    mat3 = mat1*mat2;
    endTime = clock();
    cout<< "runtime of operation #3 : " \
        <<(double)(endTime-startTime)/CLOCKS_PER_SEC<< "s"<<endl;

    return 0;
}