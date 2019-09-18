#ifndef COMMON_MATH
#define COMMON_MATH

#include "common.h"

typedef int INTE;
typedef double REAL;

typedef VectorXd POINT;
typedef Vector2d POINT2;
typedef Vector3d POINT3;
typedef VectorXd VEC;
typedef MatrixXd MATRIX;
typedef vector<MATRIX> MatrixArray;
typedef vector<vector<MATRIX>> Tens4;

const REAL PI = 3.1415926535;
template <class Tem>
vector<Tem> operator +(vector<Tem> a,vector<Tem> b){
    vector<Tem> vec(a.size());
    for(unsigned int ai=0;ai<a.size();ai++)
        vec[ai]=a[ai]+b[ai];
    return vec;
};

template <class Tem>
vector<Tem> operator -(vector<Tem> a,vector<Tem> b){
    vector<Tem> vec(a.size());
    for(unsigned int ai=0;ai<a.size();ai++)
        vec[ai]=a[ai]-b[ai];
    return vec;
};

template <class Tem>
vector<Tem> operator*(REAL k,vector<Tem> a){
    vector<Tem> vec(a.size());
    for(unsigned int i=0;i<a.size();i++)
        vec[i]=k*a[i];
    return vec;
};

template <class Tem>
Tem SumVec(vector<Tem> vec){
    for(unsigned int vi=1;vi<vec.size();vi++)
        vec[0]=vec[0]+vec[vi];
    return vec[0];
}

#endif