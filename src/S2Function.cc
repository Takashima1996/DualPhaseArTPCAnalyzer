#include "S2SingleFunc.hh"
#include "S2Function.hh"

namespace artpc
{

S2Function::S2Function(const std::vector<double>& parameter)
{
    funcParameter_=parameter;
    s2Num_ = parameter.size() / numPar_;
}

double S2Function::operator()(double x)  const
{
    std::vector<double>par(funcParameter_.size());

    double f = funcParameter_.at(0);
    
    for(int i=0;i<s2Num_;i++){
        std::vector<double>par(numPar_);
        std::copy(funcParameter_.begin() + numPar_ * i + 1,
            funcParameter_.begin() + numPar_ * (i+1), par.begin());
        
        f += S2SingleFunc(x, par);
        par.clear();
    }

     return f;
}

double S2Function::operator()(double *x, double *par) const
{
    double f = par[0];

    for(int i=0;i<s2Num_;i++){
        std::vector<double>vPar;
        for(int j=0;j<numPar_;j++){
            vPar.push_back(par[numPar_*i+1+j]);
        }
        f += S2SingleFunc(x[0], vPar);
        vPar.clear();
    }

   return f;
}

}   /* namespace artpc */
