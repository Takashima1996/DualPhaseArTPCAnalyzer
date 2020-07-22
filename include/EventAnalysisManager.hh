#ifndef COMPTONSOFT_EventAnalysisManager_H
#define COMPTONSOFT_EventAnalysisManager_H

#include <ctime>
#include <vector>
#include "ArTPCData_sptr.hh"
#include "ArDetectorHit_sptr.hh"

namespace artpc
{

class EventAnalysisManager
{
public:
    EventAnalysisManager(ArTPCData_sptr);
    ~EventAnalysisManager(){};

    void run();
    void passDetectorHit(ArDetectorHit_sptr);

private:
    ArTPCData_sptr arTPCData_;
    bool badFlag_ = false;
    int runID_= 0;
    int fileID_ = 0;
    int eventID_ = 0;
    int s2Number_ = 0;
    time_t timeStamp_ = 0;
    std::vector<double> trigTime_;
    int nBins_ = 0;
    double  FCNmin_ = 0.0;
    const int maxS2Number_ = 4;
    
    double y0_ = 0.0;
    std::vector<double> s1s2Components_;
    std::vector<double> driftT_;
    std::vector<double> amplitude_;
    std::vector<double> gasT_;
    std::vector<double> s2RatioFast_;
    std::vector<double> tau1_;
    std::vector<double> tau2_;
    std::vector<double> sigma_;
};

}   /* namespace artpc */

#endif
