#ifndef COMPTONSOFT_ArTreeIO_H
#define COMPTONSOFT_ArTreeIO_H

#include <vector>
#include <ctime>
#include <memory>
#include <TTree.h>
#include "ArDetectorHit_sptr.hh"

namespace artpc
{

class ArTreeIO
{
public:
    ArTreeIO(){};
    ~ArTreeIO(){};

    void defineBranches();
    void setBranchAddresses();
    void fillHits(const std::vector<ArDetectorHit_sptr>&);
    void retrieveHits();

private:
    int runID_ = 0;
    int fileID_ = 0;
    int eventID_=0;
    std::vector<double> trigTime_;
    long timeStamp_ = 0;
    int num_hits_=0;
    int nBins_ = 0;
    double FCNmin_ = 0.0;
    double y0_ = 0.0;
    std::vector<double> amp_;
    std::vector<double> driftT_;
    std::vector<double> gasT_;
    std::vector<double> sigma_;
    std::vector<double> tau1_;
    std::vector<double> tau2_;
    std::vector<double> s2RatioFast_;
    bool badFlag_ = false;

    std::vector<double>* amp_b_=0;
    std::vector<double>* driftT_b_=0;
    std::vector<double>* gasT_b_=0;
    std::vector<double>* sigma_b_=0;
    std::vector<double>* tau1_b_=0;
    std::vector<double>* tau2_b_=0;
    double s1Fast_ = 0.0;
    double s1Slow_ = 0.0;
    double s1Total_ = 0.0;
    double s1RatioSlow_ = 0.0;
    double s2Total_ = 0.0;
    
    std::unique_ptr<TTree> artree_;
};

} /* namespace artpc */

#endif
