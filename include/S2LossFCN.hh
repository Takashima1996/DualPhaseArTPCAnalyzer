#ifndef COMPTONSOFT_S2LossFCN_H
#define COMPTONSOFT_S2LossFCN_H 

#include <vector>

#include "Minuit2/FCNBase.h"

namespace artpc
{

class S2LossFCN : public ROOT::Minuit2::FCNBase
{
public:
    S2LossFCN();
    S2LossFCN(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<int>&, const int);
    ~S2LossFCN(){};

    virtual double Up() const {return fErrorDef_;}
    virtual double operator()(const std::vector<double>&) const;

    std::vector<double> Parameter() const {return parameter_;}
    void SetErrorDef(const double value){fErrorDef_ = value;}

private:
    int s2Num_ = 0;
    double fErrorDef_ = 1.0;
    std::vector<double>parameter_;
    std::vector<double>fData_;
    std::vector<double>fTime_;
    std::vector<double>fError_;
    std::vector<int>flag_;
    std::vector<int>flagR_;
};

}   /* namespace artpc */

#endif
