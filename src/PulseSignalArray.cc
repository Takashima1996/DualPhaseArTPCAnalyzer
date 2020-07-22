#include <iostream>
#include <cmath>
#include <numeric>
#include <boost/format.hpp>
#include <memory>

#include <TFile.h>
#include <TF1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>

#include "ArTPCData_sptr.hh"
#include "PulseSignalArray.hh"

namespace artpc
{

PulseSignalArray::PulseSignalArray(ArTPCData_sptr InputData)
{
    runID_ = InputData->RunID();
    fileID_ = InputData->FileID();
    eventID_ =InputData->EventID();
    timeStamp_ = InputData->TimeStamp();
    trigTime_.resize(nch_);
    trigTime_ = InputData->TrigTime();
    timeWindow_ = InputData->TimeWindow();
    tempTP_ = InputData->TriggerPosition();
    FADCCount_.resize(nch_);
    FADCCount_ = InputData->FADCCount();
    std::vector<double> zerosDouble((int)timeWindow_,0.0);
    std::vector<double> zerosDoubleR((int)timeWindow_/toBin_,0.0);
    std::vector<double> zerosDoubleCh(nch_, 0.0);

    for(int ibin=0;ibin<(int)timeWindow_;ibin++){
        timeSeries_.push_back(FADC_Sampling_ * (ibin+1/2) + minT_);
    }
    for(int ibin=0;ibin<(int)timeWindow_/toBin_;ibin++){
        timeSeriesR_.push_back(FADC_Sampling_ * toBin_ * (ibin+1/2) + minT_);
    }
    pulseHeightSum_ = zerosDouble;
    error_ = zerosDouble;
    pedestalMean_ = zerosDoubleCh;
    pedestalSD_ = zerosDoubleCh;
    pedestal_ = zerosDoubleCh;
    pulseHeightSumR_ = zerosDoubleR; 

    for(int ich=0;ich<nch_;ich++){
        pulseHeight_.push_back(zerosDouble);
    }
}

void PulseSignalArray::setPulseHeight()
{
    //Calculate pedestal
    for(int ich=0;ich<nch_;ich++){
        pedestal_.at(ich)=0.0;
        for(int ibin=0;ibin<(tempTP_ - Ped_to_TrigPos_);ibin++){
            pedestal_.at(ich) += FADCCount_.at(ich)->at(ibin);
        }
        pedestal_.at(ich) /= (tempTP_ - Ped_to_TrigPos_);
    }

    for(int ich =0;ich<nch_;ich++){
        for(int ibin=0;ibin<timeWindow_;ibin++){
            pulseHeight_.at(ich).at(ibin) = FADCCount_.at(ich)->at(ibin) - pedestal_.at(ich);
            pulseHeight_.at(ich).at(ibin) /= FADCtoPE_.at(ich);
        }
    }

    for(int ibin=0;ibin<timeWindow_;ibin++){
        pulseHeightSum_.at(ibin) = 0.0;
        for(int ich=0;ich<nch_;ich++){
            pulseHeightSum_.at(ibin) += pulseHeight_.at(ich).at(ibin);
        }
    }

    for(int i=0;i<timeWindow_/toBin_;i++){
        pulseHeightSumR_.at(i) = 
            std::accumulate(pulseHeightSum_.begin()+i*toBin_, 
                pulseHeightSum_.begin()+(i+1)*toBin_-1, 0.0);
    }
}

void PulseSignalArray::pedestalAnalysis()
{
    std::unique_ptr<TH1D> pedHist(new TH1D("h","h",20,-2,2));

    for(int ibin=0;ibin<(tempTP_-Ped_to_TrigPos_)/toBin_;ibin++){
        pedHist->Fill(pulseHeightSumR_.at(ibin));
    }
    sumSDR_ = pedHist->GetStdDev();
}

void PulseSignalArray::errorCalculate()
{
    pedestalAnalysis();

    for(int ibin=0;ibin<timeWindow_;ibin++){
        if(std::sqrt(pulseHeightSum_.at(ibin)) > sumSD_){
            error_.at(ibin) = std::sqrt(pulseHeightSum_.at(ibin));
        }   else if(std::sqrt(pulseHeightSum_.at(ibin))<sumSD_){
            error_.at(ibin) = sumSD_;
        }
    }

    if(errorR_.size()==0){
        std::vector<double> zeros(timeWindow_/toBin_, 0.0);
        errorR_ = zeros;
    }

    for(int ibin=0;ibin<timeWindow_/toBin_;ibin++){
        //if(std::sqrt(pulseHeightSumR_.at(ibin)) > sumSDR_*3.0){
        if(std::sqrt(pulseHeightSumR_.at(ibin)) > sumSDR_){
            errorR_.at(ibin) = std::sqrt(pulseHeightSumR_.at(ibin));
        }   else {
            errorR_.at(ibin) = sumSDR_;
        }
    }
}


std::vector<int> PulseSignalArray::EventInfo()
{
    std::vector<int> info{runID_, fileID_, eventID_};
    return info;
}

//for test
void PulseSignalArray::makeHistogram()
{
    setPulseHeight();

    dataHistogram_.reset(new TH1D("hist","", timeWindow_, minT_, maxT_));
    for(int ibin=0;ibin<timeWindow_;ibin++){
        dataHistogram_->SetBinContent(ibin+1, pulseHeightSum_.at(ibin));
    }
    dataHistogram_->SetStats(0);
    dataHistogram_->GetXaxis()->SetTitle("Time[#mus]");
    dataHistogram_->GetYaxis()->SetTitle("Photoelectron[p.e.]");
}

void PulseSignalArray::saveHistogram()
{
    makeHistogram();
    dataHistogram_->SaveAs("test.root");
}

}   /* namespace artpc */
