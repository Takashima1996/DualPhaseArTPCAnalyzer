#ifndef COMPTONSOFT_PulseSignalArray_H
#define COMPTONSOFT_PulseSignalArray_H

#include <vector>
#include <ctime>

#include <TH1.h>
#include <TString.h>

#include "ArTPCData_sptr.hh"

namespace artpc
{

class PulseSignalArray
{
private:
    int runID_ = 0;
    int fileID_ = 0;
    int eventID_ = 0;
    time_t timeStamp_ = 0;
    std::vector<double> trigTime_;

    double initialTime_ = 0.0;

    const double FADC_Sampling_ = 4.0e-3; //250MHz
    const int nch_ = 7 * 2;
    const int Ped_to_TrigPos_ = 250; 
    int CurrentEvent_ = -1;
    const int toBin_ = 10;
    int timeWindow_;
    int tempTP_ = 0;
    
    const double minT_ = -10.0;
    const double maxT_ = 130.0;
    double sumSD_ = 0.0;
    double sumSDR_ = 0.0;
    double pedestalSDSum_=0.0; 

    std::unique_ptr<TH1D> dataHistogram_;

    std::vector<std::vector<double>> pulseHeight_;
    std::vector<std::vector<long>*> FADCCount_;
    std::vector<double> timeSeries_;
    std::vector<double> timeSeriesR_;
    std::vector<double> pulseHeightSum_;
    std::vector<double> pulseHeightSumR_;
    std::vector<double> error_;
    std::vector<double> errorR_;
    std::vector<double> pedestal_;
    std::vector<double> pedestalMean_;
    std::vector<double> pedestalSD_;

    const std::vector<double>FADCtoPE_{50.0,50.0,50.0,50.0,50.0,50.0,40.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0};

public:
    PulseSignalArray(ArTPCData_sptr);
    ~PulseSignalArray(){};

    void makeTimeSeries();
    void changeEvent();
    void pedestalCalculate();
    void setPulseHeight();
    void pedestalAnalysis();
    void errorCalculate();
    void makeHistogram();
    void saveHistogram();

    std::vector<double> Error()  {return error_;}
    std::vector<double> ErrorR() {return errorR_;}
    std::vector<double> getPulseSignalArray()   {return pulseHeightSum_;}
    std::vector<double> getPulseSignalArrayR()   {return pulseHeightSumR_;}
    std::vector<double> TimeSeries()   {return timeSeries_;}
    std::vector<double> TimeSeriesR()   {return timeSeriesR_;}
    std::vector<int> EventInfo();
    int TimeWindow()    {return timeWindow_;}
    int RunID()  {return runID_;}
    int FileID() {return fileID_;}
    int EventID()    {return eventID_;}
    time_t TimeStamp()  {return timeStamp_;}
    std::vector<double> TrigTime() {return trigTime_;}
    double SumSDR() {return sumSDR_;}

};

} /*namespace artpc */

#endif
