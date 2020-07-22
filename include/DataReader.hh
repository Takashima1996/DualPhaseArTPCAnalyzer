#ifndef COMPTONSOFT_DataReader_H
#define COMPTONSOFT_DataReader_H

#include <boost/format.hpp>
#include <cmath>
#include <vector>
#include <memory>
#include <TString.h>
#include <TChain.h>
#include <ctime>
#include "ArTPCData_sptr.hh"

namespace artpc
{

using std::vector;

typedef struct {
  // General Information
  uint32_t FUNCTION_ID;
  time_t TIME_STAMP;

  // Run Configuration
  uint32_t NEVENT;
  uint32_t NEVENT_IN_BANK;
  uint32_t NBANK_IN_FILE;
  uint32_t TIME_WINDOW;
  uint32_t TRIGGER_POSITION;
  uint32_t TRIGGER_MASK;
  uint32_t SELF_TRIGGER_BIT;
  uint32_t COINCIDENCE_LEVEL;
  uint32_t COINCIDENCE;
  char COINCIDENCE_TABLE_FILE[64];
  uint32_t COINCIDENCE_GATE_LENGTH;
  uint32_t HE_SUPPRESS_BIT;
  uint32_t THRESHOLD[16 * 2];
  uint32_t HE_THRESHOLD[16 * 2];
  uint32_t FREQ;
  uint32_t PEAKING_TIME;
  uint32_t GAP_TIME;
  uint32_t CHANNEL_INV_BIT;

  // Board Information
  uint32_t MODULE_ID[2];
  uint32_t SERIAL_NUMBER[2];
  uint32_t HARDWARE_VER[2];
  uint32_t FARMWARE_VER_ADCG[4 * 2];

  // Dummy
  uint32_t DUMMY[32];

} ConfData;

typedef struct {
  uint32_t RUN_NUMBER;
  uint32_t FILE_NUMBER;
  uint32_t EVENT_NUMBER;
  uint32_t EVENT_SIZE;
  uint32_t CHANNEL_MASK;
  uint32_t EVENT_COUNTER;
  uint32_t TRIGGER_TIME_TAG_FIRST[32];
  uint32_t TRIGGER_TIME_TAG_SECOND[32];
  uint32_t TRIGGER_TIME_TAG_THIRD[32];
  uint32_t MAW_TEST_FLAG;
  uint32_t STATUS_FLAG;
  uint32_t DUMMY_F[32];
  uint32_t DUMMY_S[32];
} HeaderData;

class DataReader
{
public:
    DataReader(TString InputRootFile);
    ~DataReader(){};

    void setArTPCData(ArTPCData_sptr);
    void setEvent(const int event);

private:
    int runID_ = 0;
    int fileID_ = 0;
    int eventID_ = 0;
    int nch_ = 7 * 2;
    std::shared_ptr<TChain> chain_;
    vector<double> trigTime_;
    vector<vector<long>*> FADCRawCount_;
    ConfData conf_;
    HeaderData header_;
};

DataReader::DataReader(TString InputRootFile)
{
    chain_.reset(new TChain("rawevt"));
    chain_->Add(InputRootFile);
    chain_->SetBranchAddress("conf", &conf_);
    chain_->SetBranchAddress("header", &header_);
    chain_->GetEntry(0);

    std::vector<long>* zerosLong = 0;
    FADCRawCount_.clear();
    for(int ich=0;ich<nch_;ich++){
        FADCRawCount_.push_back(zerosLong);
        double TrigTime=(header_.TRIGGER_TIME_TAG_THIRD[ich] *
            std::pow(2,32)+header_.TRIGGER_TIME_TAG_SECOND[ich]
            * std::pow(2,16)+header_.TRIGGER_TIME_TAG_FIRST[ich])*4.0e-9;
        trigTime_.push_back(TrigTime);
    }
    for(int ich=0;ich<nch_;ich++){
        chain_->SetBranchAddress((boost::format("ch%1%") %(ich+1)).str().c_str(),
            &(FADCRawCount_.at(ich)));
    }
    chain_->GetEntry(0);
}

void DataReader::setArTPCData(ArTPCData_sptr arTPCData)
{
    arTPCData->setRunID(header_.RUN_NUMBER);
    arTPCData->setFileID(header_.FILE_NUMBER);
    arTPCData->setEventID(eventID_);
    arTPCData->setTimeWindow(conf_.TIME_WINDOW);
    arTPCData->setTriggerPosition(conf_.TRIGGER_POSITION);
    arTPCData->setTimeStamp(conf_.TIME_STAMP);
    arTPCData->setTrigTime(trigTime_);

    arTPCData->setFADCCount(FADCRawCount_);
}

void DataReader::setEvent(const int event)
{
    eventID_ = event;
    chain_->GetEntry(event);
}

} /* namespace artpc */

#endif
