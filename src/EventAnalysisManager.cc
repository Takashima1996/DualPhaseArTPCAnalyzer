#include <vector>

#include "ArTPCData_sptr.hh"
#include "DPArTPCEventAnalysis.hh"
#include "ArDetectorHit.hh"
#include "ArDetectorHit_sptr.hh"
#include "S2FitParHolder.hh"
#include "EventAnalysisManager.hh"

namespace artpc
{

EventAnalysisManager::EventAnalysisManager(ArTPCData_sptr arTPCData)
{
    arTPCData_ = arTPCData;
    runID_ = arTPCData_->RunID();
    fileID_ = arTPCData_->FileID();
    eventID_ = arTPCData_->EventID();
    timeStamp_ = arTPCData_->TimeStamp();
    trigTime_ = arTPCData_->TrigTime();
}

void EventAnalysisManager::run()
{
    DPArTPCEventAnalysis analyzer(arTPCData_);
    analyzer.s1s2Integral();
    s2Number_ = analyzer.pulseDetect();
    s1s2Components_.resize(s2Number_);
    s1s2Components_ = analyzer.S1S2Components();
    
    if((s2Number_>0) && (s2Number_<=maxS2Number_)){
        S2FitParHolder s2FitParHolder;
        analyzer.s2AutoFit();
        analyzer.passFitResult(s2FitParHolder);
        badFlag_ = analyzer.s1BadEventSelect();

        driftT_.resize(s2Number_);
        amplitude_.resize(s2Number_);
        gasT_.resize(s2Number_);
        s2RatioFast_.resize(s2Number_);
        tau1_.resize(s2Number_);
        tau2_.resize(s2Number_);
        sigma_.resize(s2Number_);

        y0_ = s2FitParHolder.y0();
        driftT_ = s2FitParHolder.DriftT();
        amplitude_ = s2FitParHolder.Amplitude();
        gasT_ = s2FitParHolder.GasT();
        s2RatioFast_ = s2FitParHolder.S2RatioFast();
        tau1_ = s2FitParHolder.Tau1();
        tau2_ = s2FitParHolder.Tau2();
        sigma_ = s2FitParHolder.Sigma();
        nBins_ = s2FitParHolder.NBins();
        FCNmin_ = s2FitParHolder.FCNMin();
    }
}

void EventAnalysisManager::passDetectorHit(ArDetectorHit_sptr arDetectorHit)
{
    arDetectorHit->setRunID(runID_);
    arDetectorHit->setFileID(fileID_);
    arDetectorHit->setEventID(eventID_);
    arDetectorHit->setTrigTime(trigTime_);
    arDetectorHit->setTimeStamp(timeStamp_);
    arDetectorHit->setS2Number(s2Number_);
    arDetectorHit->sety0(y0_);
    arDetectorHit->setAmplitude(amplitude_);
    arDetectorHit->setDriftTime(driftT_);
    arDetectorHit->setGasTime(gasT_);
    arDetectorHit->setS2RatioFast(s2RatioFast_);
    arDetectorHit->setTau1(tau1_);
    arDetectorHit->setTau2(tau2_);
    arDetectorHit->setSigma(sigma_);
    arDetectorHit->setFCNMin(FCNmin_);
    arDetectorHit->setNBins(nBins_);
    arDetectorHit->setS1S2Components(s1s2Components_);
    arDetectorHit->setBadFlag(badFlag_);
}

}   /* namespace artpc  */
