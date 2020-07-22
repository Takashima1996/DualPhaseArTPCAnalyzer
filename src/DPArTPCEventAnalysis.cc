#include <cmath>
#include <numeric>
#include <boost/format.hpp>
#include <memory>

#include <TTree.h>
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnPrint.h"

#include "PulseSignalArray.hh"
#include "DPArTPCEventAnalysis.hh"
#include "ArTPCData_sptr.hh"
#include "S2LossFCN.hh"
#include "S2Function.hh"
#include "S2FitParHolder.hh"

namespace artpc
{

DPArTPCEventAnalysis::DPArTPCEventAnalysis(ArTPCData_sptr arTPCData)
{
    phaseSignalArray_.reset(new PulseSignalArray(arTPCData));
    phaseSignalArray_->setPulseHeight();
    phaseSignalArray_->errorCalculate();

    timeWindow_ = phaseSignalArray_->TimeWindow();
    runID_ = phaseSignalArray_->RunID();
    eventID_ = phaseSignalArray_->EventID();
    fileID_ = phaseSignalArray_->FileID();
    timeStamp_ = phaseSignalArray_->TimeStamp();
    trigTime_.resize(nch_);
    trigTime_ = phaseSignalArray_->TrigTime();
}

void DPArTPCEventAnalysis::makeHistogram()
{
    phaseSignalArray_->setPulseHeight();
    phaseSignalArray_->errorCalculate();
    errorR_ = phaseSignalArray_->ErrorR();
    error_ = phaseSignalArray_->Error();

    std::vector<double>dataArray = phaseSignalArray_->getPulseSignalArray();
    std::vector<double>errorArray = phaseSignalArray_->Error();
    std::vector<double>errorArrayR = phaseSignalArray_->ErrorR();

    histogramSum_.reset();
    histogramSum_.reset(new TH1D("sum","",timeWindow_, minT_, maxT_));
    histogramSum_->SetStats(0);
    histogramSum_->GetXaxis()->SetTitle("Time[#mus]");
    histogramSum_->GetYaxis()->SetTitle("Photoelectron[p.e.]");
    for(int ibin=0;ibin<timeWindow_;ibin++)
    {
        histogramSum_->SetBinContent(ibin+1,dataArray.at(ibin));
        histogramSum_->SetBinError(ibin+1, errorArray.at(ibin));
    }

    histogramSumR_.reset((TH1D*)histogramSum_->Clone("sum_rebin"));
    histogramSumR_->Rebin(toBin_);
    
    for(int ibin=0;ibin<timeWindow_ / toBin_;ibin++){
        histogramSumR_->SetBinError(ibin+1, errorArrayR.at(ibin));
    }
}

//hyperparameters are binElapse & level
int DPArTPCEventAnalysis::pulseDetect()
{
    makeHistogram();
    std::unique_ptr<TH1D> histogramRR;
    histogramRR.reset((TH1D*)histogramSumR_->Clone("RR"));
    histogramRR->Rebin(toBin_);

    const int bin_1us = histogramRR->FindBin(-1.0);
    double stdRR = 0.0;
    for(int ibin=1; ibin<=bin_1us; ibin++){
        double val = histogramRR->GetBinContent(ibin);
        stdRR += val * val;
    }
    stdRR = std::sqrt(stdRR / bin_1us);

    const double level = 10.0 * stdRR;

    const int bin5us = histogramRR->FindBin(5.0);
    const int bin110us = histogramRR->FindBin(110.0);
    int s2number = 0;
    bool flag = false;  //true: in pulse false: out of pulse
    int lastPulseBin = 1;
    const int binElapse = 3.0 / (FADC_Sampling_ * toBin_ * toBin_);
    for(int ibin=bin5us; ibin< bin110us; ibin++){
        double val1 = histogramRR->GetBinContent(ibin);
        double val2 = histogramRR->GetBinContent(ibin+1);
        double val3 = histogramRR->GetBinContent(ibin+2);
        double val4 = histogramRR->GetBinContent(ibin+3);
        bool sign = (val1>level)&&(val2>level)&&(val3>level)&&(val4>level);
        bool rising = (val1<val2)&&(val2<val3)&&(val3<val4);
        bool elapse = (ibin>(lastPulseBin + binElapse));

       if((flag==false)&&sign&&rising){
           s2number += 1;
           flag = true;
           lastPulseBin = ibin;
       }    else if(flag && elapse && rising){
           s2number += 1;
           lastPulseBin = ibin;
       }    else if(flag && elapse && (val1<level) && (val2<level) && (val3<level)){
           flag=false;
       }
    }

    //from here, eliminating an event in which S2 occurs after 110us
    for(int ibin=bin110us+1; ibin<= timeWindow_/(toBin_*toBin_); ibin++){
        double val = histogramRR->GetBinContent(ibin);
        if(val > level){
            s2number = 0;
            break;
        }
    }

    return s2number;
}

void DPArTPCEventAnalysis::s1s2Integral()
{
    makeHistogram();

    const int S1EndBin = histogramSum_->GetXaxis()->FindBin(s1Range_/2.0);
    std::unique_ptr<TH1D> temp;
    temp.reset((TH1D*)histogramSum_->Clone("narrow"));
    const int s1RangeBin = static_cast<int>(s1Range_ / FADC_Sampling_);

    std::vector<double>photoElectron(s1RangeBin+1,0);

    temp->GetXaxis()->SetRange(1, S1EndBin);
    const int Maxbin= temp->GetMaximumBin();

    std::vector<double>tempArray = phaseSignalArray_->getPulseSignalArray();

    std::copy(tempArray.begin() + Maxbin - 1 - s1RangeBin/2,
    tempArray.begin()+Maxbin-1+s1RangeBin/2, photoElectron.begin());

    const int bin004us = 0.04 / FADC_Sampling_;
    const int bin018us = 0.18 / FADC_Sampling_;
    const int bin3us = 3.0 / FADC_Sampling_;
    s1Fast_ = std::accumulate(photoElectron.begin() + s1RangeBin/2 - bin004us, 
        photoElectron.begin() + s1RangeBin/ 2 + bin018us, 0.0);
    s1Slow_ = std::accumulate(photoElectron.begin() + s1RangeBin/2 + bin018us+1,
        photoElectron.begin() + s1RangeBin/2 + bin3us, 0.0);
    s1Total_ = s1Fast_ + s1Slow_;
    s1Ratio_ = s1Slow_ / s1Total_;

    const int bin110us = 110.0 / FADC_Sampling_;
    
    const int binS2Start = histogramSum_->FindBin(3.0);
    const int binS2End = histogramSum_->FindBin(110.0);
    s2Total_ = histogramSum_->Integral(binS2Start, binS2End);
}

bool DPArTPCEventAnalysis::s1BadEventSelect()
{
    s1s2Integral();

    if((s1Total_>600.0)&&(s1Total_<3000.0) && (0.58<s1Ratio_) && (s1Ratio_<0.65)){
        const double peIntegral = histogramSum_->Integral(1, timeWindow_);
        if(peIntegral < 200000.0)  return false;
    }

    return true;
}

void DPArTPCEventAnalysis::s2FitInitPars(const int s2Num)
{
    initPars_.clear();
    makeHistogram();

    for(int i=0;i<(s2Num * Numpar_+1);i++){
        initPars_.push_back(0.0);
    }

    const int bin5us = histogramSumR_->FindBin(5.0);
    const int bin110us = histogramSumR_->FindBin(110.0);
    std::unique_ptr<TH1D> copy((TH1D*)histogramSumR_->Clone("copy"));
    copy->GetXaxis()->SetRange(bin5us,bin110us);

    std::vector<double> MaxValBin, MaxValBinCopy, MaxVal;
    std::vector<int> overlap;

    for(int i=0; i<s2Num;i++){
        MaxValBin.push_back(copy->GetMaximumBin());
        MaxVal.push_back(copy->GetMaximum());
        const double MaxTime = (MaxValBin.at(i) - 1) * FADC_Sampling_ * toBin_ + minT_;
        const int binStart = copy->FindBin(MaxTime - 3.0);
        const int binEnd = copy->FindBin(MaxTime + 5.0);
        for(int l = binStart; l<=binEnd; l++){
            copy->SetBinContent(l,0);
        }
    }
    MaxValBinCopy = MaxValBin;
    std::sort(MaxValBinCopy.begin(), MaxValBinCopy.end());

    const double critLenBin = 10.0 / (FADC_Sampling_ * toBin_);
    //
    
    if(s2Num == 1)  overlap.push_back(0);
    if(s2Num == 2){
        if((MaxValBinCopy.at(1)-MaxValBinCopy.at(0))  < critLenBin)
        {
            overlap.push_back(1);
            overlap.push_back(2);
        }   else 
        {
            overlap.push_back(0);
            overlap.push_back(0);
        }
    }
    if (s2Num > 2) {
        if((MaxValBinCopy.at(1) - MaxValBinCopy.at(0)) < critLenBin){
            overlap.push_back(1);
        }   else{
            overlap.push_back(0);
        }
        for(int i=1;i<(s2Num-1); i++){
            if(overlap.at(i-1)==0){
                if((MaxValBinCopy.at(i+1) - MaxValBinCopy.at(i)) < critLenBin){
                    overlap.push_back(1);
                }   else{
                    overlap.push_back(0);
                }
            }
            else if(overlap.at(i-1)==1){
                if(((MaxValBinCopy.at(i+1) - MaxValBinCopy.at(i))) < critLenBin){
                    overlap.push_back(3);
                }   else{
                    overlap.push_back(2);
                }
            }
            else if(overlap.at(i-1)==2){
                if((MaxValBinCopy.at(i+1) - MaxValBinCopy.at(i)) < critLenBin){
                    overlap.push_back(1);
                }   else{
                    overlap.push_back(0);
                }
            }
            else if(overlap.at(i-1)==3){
                if((MaxValBinCopy.at(i+1) - MaxValBinCopy.at(i)) < critLenBin){
                    overlap.push_back(3);
                }   else{
                    overlap.push_back(2);
                }
            }

        }
        if((overlap.at(s2Num-2)==1)||(overlap.at(s2Num-2)==3)){
            overlap.push_back(2);
        }   else{
            overlap.push_back(0);
        }
    }

    initPars_.at(0);
    for(int i=0;i<s2Num;i++){
        const double MaxTime = MaxValBin.at(i) * FADC_Sampling_ * toBin_ + minT_;
        const int binStart = copy->FindBin(MaxTime-5.0);
        const int binStartShort = copy->FindBin(MaxTime- 1.0);
        const int binEnd = copy->FindBin(MaxTime + 8.0);
        const int binEndShort = copy->FindBin(MaxTime + 3.0);

        double A=0.0;
        if(overlap.at(i) == 0){
            A = histogramSumR_->Integral(binStart, binEnd) * FADC_Sampling_ * toBin_;
        }   else if(overlap.at(i) == 1){
            A = histogramSumR_->Integral(binStart, binEndShort) * FADC_Sampling_ * toBin_;
        }   else if(overlap.at(i) == 2){
            A = histogramSumR_->Integral(binStartShort, binEnd) * FADC_Sampling_ * toBin_;
        }   else{
            A = histogramSumR_->Integral(binStartShort, binEndShort) * FADC_Sampling_ * toBin_;
        }
        
        initPars_.at(1+i*Numpar_)=A;
        initPars_.at(2+i*Numpar_)=1.1;
        initPars_.at(3+i*Numpar_)=MaxTime-1.0;
        initPars_.at(4+i*Numpar_)=0.3;
        initPars_.at(5+i*Numpar_)=0.011;
        initPars_.at(6+i*Numpar_)=0.20;
        initPars_.at(7+i*Numpar_)=3.151;
        initPars_.at(8+i*Numpar_)=MaxTime - 5.0;
        initPars_.at(9+i*Numpar_)=MaxTime + 8.0;
    }
}

void DPArTPCEventAnalysis::flagCalculate(const int s2Num, bool limited=true)
{
    flagR_.clear();
    std::vector<int> zerosR(timeWindow_ / toBin_, 0);
    flagR_ = zerosR;

    std::vector<double>initStartBins, initEndBins;
    for(int i=0;i<s2Num;i++){
        initStartBins.push_back(
            histogramSumR_->FindBin(initPars_.at(8 + i * Numpar_)));
        initEndBins.push_back(
            histogramSumR_->FindBin(initPars_.at(9 + i * Numpar_)));
    }

    std::sort(initStartBins.begin(), initStartBins.end());
    std::sort(initEndBins.begin(), initEndBins.end());

    for(int i=initStartBins.at(0); i<=initEndBins.at(s2Num-1); i++){
        flagR_.at(i-1) = 1;
    }
    if(limited==true){
        for(int j=0;j<(s2Num-1);j++){
            if((initEndBins.at(j)+1) < initStartBins.at(j+1)){
                for(int l=(initEndBins.at(j)+1); l<initStartBins.at(j+1); l++){
                    flagR_.at(l-1) = 0;
                }
            }
        }
    }
}


void DPArTPCEventAnalysis::s2Fit(const int s2Num, const bool limited=true)
{
    const int npar = Numpar_ * s2Num + 1;
    if(limited==true){
        s2FitInitPars(s2Num);
    }
    if(limited==false) {
        initPars_.clear();
        initPars_ = parameter_;
    }
    flagCalculate(s2Num, limited);
    
    std::vector<double>dataR, timeR;

    for(int i=0;i<histogramSumR_->GetNbinsX();i++){
        dataR.push_back(histogramSumR_->GetBinContent(i+1));
        timeR.push_back(minT_ + FADC_Sampling_*toBin_*(i+0.5));
    }

    S2LossFCN fFCN(dataR, timeR, errorR_, flagR_, s2Num);
    ROOT::Minuit2::MnUserParameters mnPars;

    mnPars.Add("y0",initPars_.at(0),step_.at(0),-5.0,5.0);
    for(int i=0;i<s2Num;i++){
        mnPars.Add((boost::format("A_%1%") %(i+1)).str().c_str(), initPars_.at(Numpar_*i+1), step_.at(1));
        mnPars.Add((boost::format("T_%1%") %(i+1)).str().c_str(),initPars_.at(Numpar_*i+2),step_.at(2));
        mnPars.Add((boost::format("t0_%1%") %(i+1)).str().c_str(),initPars_.at(Numpar_*i+3),step_.at(3),initPars_.at(Numpar_*i+3)-5, initPars_.at(Numpar_*i+3)+5);
        mnPars.Add((boost::format("sigma_%1%") %(i+1)).str().c_str(),initPars_.at(Numpar_*i+4));
        mnPars.Add((boost::format("tau1_%1%") %(i+1)).str().c_str(),initPars_.at(Numpar_*i+5));
        mnPars.Add((boost::format("p_%1%") %(i+1)).str().c_str(),initPars_.at(Numpar_*i+6),step_.at(6));
        mnPars.Add((boost::format("tau2_%1%") %(i+1)).str().c_str(),initPars_.at(Numpar_*i+7),step_.at(7));
        mnPars.Add((boost::format("FirStart_%1%") %(i+1)).str().c_str(),initPars_.at(Numpar_*i+8));
        mnPars.Add((boost::format("FitEnd_%1%") %(i+1)).str().c_str(),initPars_.at(Numpar_*i+9));

        //mnPars.Fix((boost::format("T_%1%") %(i+1)).str().c_str());
        mnPars.Fix((boost::format("sigma_%1%") %(i+1)).str().c_str());
        mnPars.Fix((boost::format("tau1_%1%") %(i+1)).str().c_str());
        mnPars.Fix((boost::format("tau2_%1%") %(i+1)).str().c_str());
        mnPars.Fix((boost::format("FirStart_%1%") %(i+1)).str().c_str());
        mnPars.Fix((boost::format("FitEnd_%1%") %(i+1)).str().c_str());
    }
    std::cout<<"Initial parameter:\n";
    for(int i=0;i<initPars_.size();i++){
        std::cout<<initPars_.at(i)<<std::endl;
    }

    ROOT::Minuit2::MnMigrad migrad(fFCN, mnPars);
    ROOT::Minuit2::FunctionMinimum min = migrad();
    ROOT::Minuit2::MnUserParameterState state = min.UserState();
    std::cout<< min << std::endl;

    if(bminos_==true){
        ROOT::Minuit2::MnMinos minos(fFCN, min);
        for(int i=0; i<initPars_.size()/Numpar_ ;i++){
            std::pair<double, double>eA=minos(i*Numpar_ + 1);
            std::cout<<"A"<<i+1<<"error : "<<eA.first<<" "<<eA.second<<std::endl;
            std::pair<double, double>et0=minos(i*Numpar_ + 3);
            std::cout<<"t0"<<i+1<<"error : "<<et0.first<<" "<<et0.second<<std::endl;
        }
    }

    parameter_.clear();
    for(int ipar = 0;ipar<npar; ipar++){
        parameter_.push_back(state.Value(ipar));
    }
    FCNmin_ = min.Fval();
    nBins_ = std::accumulate(flagR_.begin(), flagR_.end(), 0);
    chi2_norm_ = FCNmin_ / nBins_;
}

bool DPArTPCEventAnalysis::s2AutoFit()
{
    int s2Number = pulseDetect();
    if(s2Number==0) return false;
    s2Refit(s2Number);
    return true;
}

//Using all range 
void DPArTPCEventAnalysis::s2Refit(const int s2Num)
{
    s2Fit(s2Num);
    s2Fit(s2Num, false);
}

void DPArTPCEventAnalysis::s2Chi2()
{
    vChi2_.clear();
    parMemory_.clear();
    const bool limited=false;
    for(int k=1; k<=maxS2Num_; k++){
        s2Refit(k);
        vChi2_.push_back(FCNmin_);
        parMemory_.push_back(parameter_);
    }
    for(int k=0;k<maxS2Num_;k++){
        std::cout<<k+1<<":"<<vChi2_.at(k)<<std::endl;
    }
}

void DPArTPCEventAnalysis::s2BIC(const double lambda)
{
    vBIC_.clear();
    s2Chi2();
    nBins_ = std::accumulate(flagR_.begin(), flagR_.end(), 0);

    int bestNum = 1;
    for(int i=0;i<maxS2Num_;i++){
        const double BIC = 
            vChi2_.at(i) + lambda * (4*(i+1)+1) * std::log(nBins_);
        vBIC_.push_back(BIC);
        std::cout<<"BIC(k="<<i+1<<"):"<<BIC<<std::endl;
        if((i>0) && (vBIC_.at(i)<vBIC_.at(bestNum-1)))   bestNum = i+1;
    }

    bestParameter_.clear();
    bestParameter_ = parMemory_.at(bestNum-1);
}

std::vector<double> DPArTPCEventAnalysis::S1S2Components()
{
    std::vector<double>components{s1Fast_, s1Slow_, s2Total_};

    return components;
}

void DPArTPCEventAnalysis::passFitResult(S2FitParHolder& s2ParHolder)
{

    int s2num = parameter_.size() / Numpar_;
    std::vector<double>vec_A, vec_t0, vec_sigma, vec_T,  vec_tau1, vec_tau2, vec_p;
    for(int ipar=0; ipar<s2num; ipar++){
        vec_A.push_back(parameter_.at(1 + ipar*Numpar_));
        vec_T.push_back(parameter_.at(2 + ipar*Numpar_));
        vec_t0.push_back(parameter_.at(3 + ipar*Numpar_));
        vec_sigma.push_back(parameter_.at(4 + ipar*Numpar_));
        vec_tau1.push_back(parameter_.at(5 + ipar*Numpar_));
        vec_p.push_back(parameter_.at(6 + ipar*Numpar_));
        vec_tau2.push_back(parameter_.at(7 + ipar*Numpar_));
    }

    s2ParHolder.setRunID(runID_);
    s2ParHolder.setFileID(fileID_);
    s2ParHolder.setEventID(eventID_);
    s2ParHolder.setS2Number(s2num);

    s2ParHolder.setNBins(nBins_);
    s2ParHolder.setFCNMin(FCNmin_);
    s2ParHolder.sety0(parameter_.at(0));
    s2ParHolder.setA(vec_A);
    s2ParHolder.sett0(vec_t0);
    s2ParHolder.setSigma(vec_sigma);
    s2ParHolder.setT(vec_T);
    s2ParHolder.setTau1(vec_tau1);
    s2ParHolder.setTau2(vec_tau2);
    s2ParHolder.setS2RatioFast(vec_p);
}

}   /* namespace artpc */
