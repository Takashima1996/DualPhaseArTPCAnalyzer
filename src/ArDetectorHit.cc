#include <vector>
#include "ArDetectorHit.hh"

namespace artpc
{

void ArDetectorHit::setAmplitude(std::vector<double> Amp)
{
    amplitude_.resize(Amp.size()); 
    amplitude_= Amp;
}

void ArDetectorHit::setDriftTime(std::vector<double> driftT)
{
    driftT_.resize(driftT.size());
    driftT_ = driftT;
}

void ArDetectorHit::setGasTime(std::vector<double> gasT)
{
    gasT_.resize(gasT.size());
    gasT_ = gasT;
}

void ArDetectorHit::setS2RatioFast(std::vector<double> s2RatioFast)
{
    s2RatioFast_.resize(s2RatioFast.size());
    s2RatioFast_ = s2RatioFast;
}

void ArDetectorHit::setTau1(std::vector<double> tau1)
{
    tau1_.resize(tau1.size());
    tau1_ = tau1;
}

void ArDetectorHit::setTau2(std::vector<double> tau2)
{
    tau2_.resize(tau2.size());
    tau2_ = tau2;
}

void ArDetectorHit::setSigma(std::vector<double> sigma)
{
    sigma_.resize(sigma.size());
    sigma_ = sigma;
}

void ArDetectorHit::setS1S2Components(std::vector<double> s1s2Components)
{
    s1s2Components_.resize(s1s2Components.size());
    s1s2Components_ = s1s2Components;
}

}   /* namespace artpc */
