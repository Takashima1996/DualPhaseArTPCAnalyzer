#ifndef COMPTONSOFT_ArTPCData_H
#define COMPTONSOFT_ArTPCData_H

#include <vector>
#include <memory>
#include <ctime>

using std::vector;

namespace artpc
{

class ArTPCData
{
public:
    ArTPCData() {FADCCount_.resize(nch_);}
    ~ArTPCData(){};

    void setFADCCount(vector<vector<long>*> &RawData)   
    {
        FADCCount_ = RawData;
    }
    void setTimeWindow(int window)    {timeWindow_ = window;}
    void setTriggerPosition(int trigpos)    {triggerPosition_ = trigpos;}
    void setRunID(int runID)   {runID_ = runID;}
    void setFileID(int fileID) {fileID_ = fileID;}
    void setEventID(int eventID)    {eventID_ = eventID;}
    void setTimeStamp(time_t timeStamp) {timeStamp_ = timeStamp;}
    void setTrigTime(vector<double> trigTime)   {trigTime_ = trigTime;}

    vector<vector<long>*> FADCCount() {return FADCCount_;}
    int RunID()  {return runID_;}
    int FileID() {return fileID_;}
    int EventID()    {return eventID_;}
    int TriggerPosition()   {return triggerPosition_;}
    int TimeWindow()    {return timeWindow_;}
    time_t TimeStamp()  {return timeStamp_;}
    vector<double> TrigTime()   {return trigTime_;}

private:
    const int nch_ = 7 * 2;
    vector<vector<long>*> FADCCount_;
    int runID_ = 0;
    int fileID_ = 0;
    int eventID_ = 0;
    int timeWindow_ = 0;
    int triggerPosition_ = 0;
    time_t timeStamp_ = 0;
    vector<double> trigTime_;
};

} /* namespace artpc */

#endif
