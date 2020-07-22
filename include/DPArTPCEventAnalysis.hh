#ifndef COMPTONSOFT_DPArTPCEventAnalysis_H
#define COMPTONSOFT_DPArTPCEventAnalysis_H

#include <iostream>
#include <vector>
#include <memory>
#include <ctime>

#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TString.h>

#include "ArTPCData_sptr.hh"
#include "PulseSignalArray.hh"
#include "S2FitParHolder.hh"

namespace artpc
{

class DPArTPCEventAnalysis
{
public:
    DPArTPCEventAnalysis(ArTPCData_sptr);
    ~DPArTPCEventAnalysis(){};

    void makeHistogram();
    void s1s2Integral();
    int pulseDetect();
    void s2FitInitPars(const int s2Num);
    void s2Fit(const int s2Num, bool limited);
    void s2Refit(const int s2Num);
    void flagCalculate(const int s2Num, const bool limited);
    bool s1BadEventSelect();
    bool s2AutoFit();
    void s2Integral();
    void s2Chi2();
    void s2BIC(const double lambda);
    void setErrorMinos(bool bminos)   {bminos_ = bminos;}

    double Chi2_norm(){return chi2_norm_;}
    double NBins()  {return nBins_;}
    double FCNMin() {return FCNmin_;}
    std::vector<double> Parameter() {return parameter_;}
    std::vector<double> S1S2Components();
    std::vector<double> S2Chi2();
    std::vector<double> S2BIC()  {return vBIC_;}
    std::vector<double> BestParameter()  {return bestParameter_;}
    std::vector<int> FlagR()    {return flagR_;}
    double S2Value();

    void passFitResult(S2FitParHolder&);

private:
    int runID_ = 0;
    int fileID_ = 0;
    int eventID_ = 0;
    int timeWindow_ = 0;
    time_t timeStamp_ = 0;
    std::vector<double> trigTime_;

    double s1Fast_=0.0;
    double s1Slow_=0.0;
    double s1Total_=0.0;
    double s1Ratio_=0.0;
    double s2Total_=0.0;
    bool bminos_ = false;

    const double FADC_Sampling_ = 4.0e-3; //250MHz
    const double epsilon_ = 1.0e-8;
    const double minT_ = -10.0;
    const double maxT_ = 130.0;
    const int nch_ = 7 * 2;

    const double s1Range_ = 10.0; //range for s1analysis[us]
    const int toBin_ = 10;
    const int Numpar_ = 9;
    const int maxS2Num_ = 4;
    int nBins_ = 0;
    double FCNmin_ = 0.0;
    double chi2_norm_ = 0.0;


    std::vector<double> chi2_;
    std::vector<double> initPars_;
    std::vector<int> flagR_;
    std::vector<double>vChi2_;
    std::vector<double>vBIC_;

    std::unique_ptr<TH1D> histogramSum_;
    std::unique_ptr<TH1D> histogramSumR_;
    std::vector<double>error_;
    std::vector<double>errorR_;
    std::vector<double> parameter_;
    std::vector<std::vector<double>> parMemory_;
    std::vector<double> bestParameter_;
    const std::vector<double> step_{0.01,0.01,0.01,0.01,0.01,0.01,0.001,0.01,0.01,0.01};

    std::unique_ptr<PulseSignalArray> phaseSignalArray_;
};

}   /* namespace artpc */

#endif
