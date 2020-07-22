#include <vector>
#include "S2FitParHolder.hh"

namespace artpc
{

void S2FitParHolder::setA(std::vector<double> &inputA)
{
    A_.resize(s2Number_);
    A_ = inputA;
}

void S2FitParHolder::setT(std::vector<double> &inputT)
{
    T_.resize(s2Number_);
    T_ = inputT;
}


void S2FitParHolder::sett0(std::vector<double> &inputt0)
{
    t0_.resize(s2Number_);
    t0_ = inputt0;
}

void S2FitParHolder::setSigma(std::vector<double> &inputSigma)
{
    sigma_.resize(s2Number_);
    sigma_ = inputSigma;
}

void S2FitParHolder::setTau1(std::vector<double> &inputTau1)
{
    tau1_.resize(s2Number_);
    tau1_ = inputTau1;
}

void S2FitParHolder::setTau2(std::vector<double> &inputTau2)
{
    tau2_.resize(s2Number_);
    tau2_ = inputTau2;
}

void S2FitParHolder::setS2RatioFast(std::vector<double> &inputP)
{
    p_.resize(s2Number_);
    p_ = inputP;
}

}   /* namespace artpc */
