#ifndef COMPTONSOFT_S2FitParHolder_H
#define COMPTONSOFT_S2FitParHolder_H

#include <vector>

namespace artpc
{

class S2FitParHolder
{
public:
    S2FitParHolder(){};
    ~S2FitParHolder(){};

    void setS2Number(int inputS2Number) {s2Number_ = inputS2Number;}
    void setRunID(int inputRunID) {runID_ = inputRunID;}
    void setFileID(int inputFileID) {fileID_ = inputFileID;}
    void setEventID(int inputEventID)   {eventID_ = inputEventID;}
    void sety0(double inputy0)  {y0_ = inputy0;}
    void setA(std::vector<double>&);
    void setT(std::vector<double>&);
    void sett0(std::vector<double>&);
    void setSigma(std::vector<double>&);
    void setTau1(std::vector<double>&);
    void setTau2(std::vector<double>&);
    void setS2RatioFast(std::vector<double>&);
    void setNBins(int nBins)    {nBins_ = nBins;}
    void setFCNMin(double FCNmin)   {FCNmin_ = FCNmin;}

    int S2Number() {return s2Number_;}
    int RunID() {return runID_;}
    int FileID() {return fileID_;}
    int EventID()   {return eventID_;}
    double y0()  {return y0_;}
    std::vector<double> Amplitude() {return A_;}
    std::vector<double> GasT() {return T_;}
    std::vector<double> DriftT() {return t0_;}
    std::vector<double> Sigma() {return sigma_;}
    std::vector<double> Tau1() {return tau1_;}
    std::vector<double> Tau2() {return tau2_;}
    std::vector<double> S2RatioFast() {return p_;}
    int NBins() {return nBins_;}
    double FCNMin() {return FCNmin_;}

private:
    int s2Number_=0;
    int runID_=0;
    int fileID_=0;
    int eventID_=0;
    int nBins_ = 0;

    double FCNmin_ = 0.0;
    double y0_ = 0.0;
    std::vector<double>A_;
    std::vector<double>T_;
    std::vector<double>t0_;
    std::vector<double>sigma_;
    std::vector<double>tau1_;
    std::vector<double>tau2_;
    std::vector<double>p_;
    std::vector<double>fitStart_;
    std::vector<double>fitEnd_;
};

}

#endif

