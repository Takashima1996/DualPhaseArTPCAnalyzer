#ifndef COMPTONSOFT_ArDetectorHit_H
#define COMPTONSOFT_ArDetectorHit_H 

#include <ctime>
#include <vector>
#include "ArTPCData_sptr.hh"

namespace artpc
{

class ArDetectorHit
{
public:
    ArDetectorHit(){};
    ~ArDetectorHit(){};

    void analyze(ArTPCData_sptr);

    int RunID() const {return runID_;}
    int FileID()    const {return fileID_;}
    int EventID()   const {return eventID_;}
    time_t TimeStamp() const  {return timeStamp_;}
    std::vector<double> TrigTime()  const  {return trigTime_;}
    int S2Number()  const {return s2Number_; }
    int NBins() const {return nBins_;}
    double FCNMin() const {return FCNmin_;}
    double y0() const {return y0_;}
    std::vector<double> Amplitude() const {return amplitude_; }
    std::vector<double> DriftTime() const {return driftT_;}
    std::vector<double> GasTime() const {return gasT_;}
    std::vector<double> S2RatioFast() const {return s2RatioFast_;}
    std::vector<double> Tau1() const {return tau1_;}
    std::vector<double> Tau2() const  {return tau2_;}
    std::vector<double> Sigma() const {return sigma_;}
    std::vector<double> S1S2Components() const {return s1s2Components_;}
    bool BadFlag() const {return badFlag_; }

    void setRunID(int runID)    {runID_ = runID;}
    void setFileID(int fileID)  {fileID_ = fileID;}
    void setEventID(int eventID)    {eventID_ = eventID;}
    void setTrigTime(std::vector<double> trigTime)   {trigTime_ = trigTime;}
    void setTimeStamp(time_t timeStamp) {timeStamp_ = timeStamp;}
    void setS2Number(int S2Num) {s2Number_ = S2Num; }
    void setNBins(int nBins)   {nBins_ = nBins;}
    void setFCNMin(double FCNmin)   {FCNmin_ = FCNmin;}
    void sety0(double y0)   {y0_ = y0;}
    void setAmplitude(std::vector<double>);
    void setDriftTime(std::vector<double>);
    void setGasTime(std::vector<double>);
    void setS2RatioFast(std::vector<double>);
    void setTau1(std::vector<double>);
    void setTau2(std::vector<double>);
    void setSigma(std::vector<double>);
    void setS1S2Components(std::vector<double>);
    void setBadFlag(bool badFlag)    {badFlag_ = badFlag;}

private:
    int runID_ = 0;
    int fileID_ = 0;
    int eventID_=0;
    std::vector<double> trigTime_;
    time_t timeStamp_ = 0;
    int  s2Number_=0;
    int nBins_ = 0;
    double FCNmin_ = 0.0;
    double y0_=0.0;
    bool badFlag_ = false;
    std::vector<double> driftT_;
    std::vector<double> amplitude_;
    std::vector<double> gasT_;
    std::vector<double> s2RatioFast_;
    std::vector<double> tau1_;
    std::vector<double> tau2_;
    std::vector<double> sigma_;
    std::vector<double> s1s2Components_;
};

}   /* namespace artpc */

#endif
