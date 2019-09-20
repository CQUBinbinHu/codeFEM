#ifndef NUMERICALMETHODS_TOOLS
#define NUMERICALMETHODS_TOOLS

#include "../../Common/common.h"
#include "../../Common/math.h"

namespace NumericalMethods{
class NumerIntegral
{
  private:

    int Dim = 1;

    const double G10=0;
    const double WG10=2.*0.5;
    
    const double G21=0.5773502691896258;
    const double WG21=1.*0.5;

    const double G31=0;
    const double WG31=8./9.*0.5;
    const double G32=0.7745966692414834;
    const double WG32=5./9.*0.5;

    const double G41=0.3399810435848563;
    const double G42=0.8611363115940526;
    const double WG41=0.6521451548625462*0.5;
    const double WG42=0.34785484513745385*0.5;

    vector<double>Gvec;
    vector<double>WGvec;
  public:
    NumerIntegral(int dim):Dim(dim){};
    NumerIntegral(int np,char TypeS){
      switch(TypeS){
        case 'L' :
          switch(np){
            case 1 :
              Gvec.push_back(G10);
              WGvec.push_back(WG10);
              break;
            case 2 :
              Gvec.push_back(-G21);
              Gvec.push_back(+G21);
              WGvec.push_back(WG21);
              WGvec.push_back(WG21);
              break;
            case 3 :
              Gvec.push_back(-G32);
              Gvec.push_back(G31);
              Gvec.push_back(+G32);
              WGvec.push_back(WG32);
              WGvec.push_back(WG31);
              WGvec.push_back(WG32);
              break;
            case 4 :
              Gvec.push_back(-G42);
              Gvec.push_back(-G42);
              Gvec.push_back(+G42);
              Gvec.push_back(+G42);
              WGvec.push_back(WG42);
              WGvec.push_back(WG41);
              WGvec.push_back(WG41);
              WGvec.push_back(WG42);
              break;
            default:
              cout<<"GaussPar:: Error02: wrong NP";
              break;
          }
          break;
        default:
          cout<<"Error in GaussPar";
      }
    };
    ~NumerIntegral(){};

    /// Gauss points for line element
    void GLPoints(MATRIX&,POINT&,POINT&);
    /// weight of Gauss points for line element
    void GLWeight(VEC&,POINT&,POINT&);
    /// gauss points for Quadracture Element
    void GQPoints(MATRIX&,POINT&,POINT&);
    /// weight of gauss points for Quadracture Element
    void GQWeight(VEC&,POINT&,POINT&);
};
///////////////////////////
/// end of namespace
///////////////////////////
};
#endif