#include <iostream>
#include <cmath>

#include "S2Function.hh"
#include "S2LossFCN.hh"

namespace artpc
{

S2LossFCN::S2LossFCN()
{
}

S2LossFCN::S2LossFCN(std::vector<double>& data, std::vector<double>& time, std::vector<double>& error, std::vector<int>& flag, const int s2Number)
{
    fData_=data;
    fTime_=time;
    fError_=error;
    flagR_=flag;
    fErrorDef_=1.0;
    s2Num_ = s2Number;
}

double S2LossFCN::operator()(const std::vector<double>& par)  const
{
    double chi2 = 0;

    S2Function s2Func(par);
    for(int i=0; i<(int)fTime_.size(); i++){
        if(flagR_.at(i)==0)  continue;
        chi2 += (s2Func(fTime_.at(i))-fData_.at(i))*(s2Func(fTime_.at(i))-fData_.at(i))/(fError_.at(i)*fError_.at(i));
    }

    return chi2;
}

}   /* namespace artpc */
