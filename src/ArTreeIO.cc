#include <iostream>
#include <vector>
#include <TTree.h>
#include <ctime>

#include "ArTreeIO.hh"
#include "S2FitParHolder_sptr.hh"
#include "ArDetectorHit_sptr.hh"

namespace artpc
{

void ArTreeIO::defineBranches()
{
    artree_.reset(new TTree("artree","artree"));
    artree_->Branch("run", &runID_, "run/I");
    artree_->Branch("file", &fileID_, "file/I");
    artree_->Branch("eventid", &eventID_, "eventid/I");
    artree_->Branch("TimeStamp", &timeStamp_, "TimeStamp/L");
    artree_->Branch("TrigTime", &trigTime_);
    artree_->Branch("num_hits",&num_hits_, "num_hits/I");
    artree_->Branch("nbins", &nBins_, "nbins/I");
    artree_->Branch("FCNmin", &FCNmin_, "FCNmin/D");
    artree_->Branch("y0", &y0_, "y0/D");
    artree_->Branch("Amplitude", &amp_);
    artree_->Branch("DriftTime", &driftT_);
    artree_->Branch("GasTime", &gasT_);
    artree_->Branch("Sigma", &sigma_);
    artree_->Branch("Tau1", &tau1_);
    artree_->Branch("Tau2", &tau2_);
    artree_->Branch("S2RatioFast", &s2RatioFast_);
    artree_->Branch("S1RatioSlow", &s1RatioSlow_);
    artree_->Branch("S1Total", &s1Total_, "S1Total/D");
    artree_->Branch("S1Slow", &s1Slow_, "S1Slow/D");
    artree_->Branch("S1Fast", &s1Fast_, "S1Fast/D");
    artree_->Branch("S2Total", &s2Total_, "S2Total/D");
    artree_->Branch("BadEvent", &badFlag_, "BadEvent/O");
}

void ArTreeIO::setBranchAddresses()
{
    artree_->SetBranchAddress("run", &runID_);
    artree_->SetBranchAddress("file", &fileID_);
    artree_->SetBranchAddress("eventid", &eventID_);
    artree_->SetBranchAddress("TrigTime", &trigTime_);
    artree_->SetBranchAddress("TimeStamp", &timeStamp_);
    artree_->SetBranchAddress("num_hits", &num_hits_);
    artree_->SetBranchAddress("nbins", &nBins_);
    artree_->SetBranchAddress("FCNmin", &FCNmin_);
    artree_->SetBranchAddress("y0", &y0_);
    artree_->SetBranchAddress("Amplitude", &amp_b_);
    artree_->SetBranchAddress("DriftTime", &driftT_b_);
    artree_->SetBranchAddress("GasTime", &gasT_b_);
    artree_->SetBranchAddress("Sigma", &sigma_b_);
    artree_->SetBranchAddress("Tau1", &tau1_b_);
    artree_->SetBranchAddress("Tau2", &tau2_b_);
    artree_->SetBranchAddress("S2RatioFast", &s2RatioFast_);
    artree_->SetBranchAddress("S1RatioSlow", &s1RatioSlow_);
    artree_->SetBranchAddress("S1Total", &s1Total_);
    artree_->SetBranchAddress("S1Slow", &s1Slow_);
    artree_->SetBranchAddress("S1Fast", &s1Fast_);
    artree_->SetBranchAddress("S2Total", &s2Total_);
    artree_->SetBranchAddress("BadEvent", &badFlag_);
}

void ArTreeIO::fillHits(const std::vector<ArDetectorHit_sptr>& hits)
{
    const int NumHits = hits.size();
    num_hits_ = NumHits;

    for(int i=0; i<NumHits; i++){
        const ArDetectorHit_sptr& hit = hits.at(i);
        runID_ = hit->RunID();
        fileID_ = hit->FileID();
        eventID_ = hit->EventID();
        trigTime_ = hit->TrigTime();
        timeStamp_ = long(hit->TimeStamp());
        num_hits_ = hit->S2Number();
        nBins_ = hit->NBins();
        FCNmin_ = hit->FCNMin();
        y0_ = hit->y0();
        amp_ = hit->Amplitude();
        driftT_ = hit->DriftTime();
        gasT_ = hit->GasTime();
        sigma_ = hit->Sigma();
        tau1_ = hit->Tau1();
        tau2_ = hit->Tau2();
        s2RatioFast_ = hit->S2RatioFast();
        std::vector<double> s1s2comp = hit->S1S2Components();
        s1Fast_ = s1s2comp[0];
        s1Slow_ = s1s2comp[1];
        s1Total_ = s1Fast_ + s1Slow_;
        s1RatioSlow_ = s1Slow_ / s1Total_;
        s2Total_ = s1s2comp[2];
        badFlag_ = hit->BadFlag();

        artree_->Fill();
    }
    artree_->SaveAs("ArHittree.root");
}

}   /* namespace artpc */
