#include "Tools.h"

void NumerIntegral::GLPoints(MATRIX& GP,POINT& p1,POINT& p2){
  GP.resize(Dim,Gvec.size());
  for(unsigned int ni=0;ni<Gvec.size();ni++){
    POINT p0;
    REAL par1=(1. - Gvec[ni])/2.;
    REAL par2=(1. + Gvec[ni])/2.;
    p0=par1*p1+par2*p2;
    GP.col(ni)=p0;
  }

};

void NumerIntegral::GLWeight(VEC& wvec,POINT&,POINT&){
  //wvec=WGvec;
  wvec.resize(WGvec.size());
  for(int i=0;i<wvec.rows();i++){
    wvec(i)=WGvec[i];
  }
};

void NumerIntegral::GQPoints(MATRIX& GP,POINT& p1,POINT& p2){
    GP.resize(2,Gvec.size()*Gvec.size());
    POINT gpoint(2);
    int iter = 0;
    for(unsigned int ni=0;ni<Gvec.size();ni++){
        for(unsigned int mi=0;mi<Gvec.size();mi++){
            gpoint[0]=(1. - Gvec[mi])/2.*p1[0]+(1. + Gvec[mi])/2.*p2[0];
            gpoint[1]=(1. - Gvec[ni])/2.*p1[1]+(1. + Gvec[ni])/2.*p2[1];
            GP.col(iter)=gpoint;
            iter++;
        }
    }
};

void NumerIntegral::GQWeight(VEC& wvec,POINT&,POINT&){
  //wvec=WGvec;
  wvec.resize(WGvec.size()*WGvec.size());
  int iter =0;
  for(unsigned long i=0;i<WGvec.size();i++){
      for (unsigned long j=0;j<WGvec.size();j++){
          wvec(iter)=WGvec[i]*WGvec[j];
          iter++;
      }
  }
};
