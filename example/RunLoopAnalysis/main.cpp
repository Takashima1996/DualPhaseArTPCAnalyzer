#include <TString.h>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

#include "ArDetectorHit.hh"
#include "ArDetectorHit_sptr.hh"
#include "ArTPCData.hh"
#include "ArTPCData_sptr.hh"
#include "ArTreeIO.hh"
#include "DataReader.hh"
#include "EventAnalysisManager.hh"

int main(int argc, char *argv[]) {
  time_t time1, time2;
  time(&time1);
  int runNum = 192;

  artpc::ArTreeIO arTreeIO;
  std::vector<artpc::ArDetectorHit_sptr> arHitArray;

  for (int ifile = 1; ifile <= 1; ifile++) {
    TString Inputfile =
        (boost::format("../../../ANKOK_RUN-%1%_%2%.root") % runNum % ifile)
            .str();
    artpc::DataReader Reader(Inputfile);
    for (int event = 0; event < 10; event++) {
      std::cout << ifile << " " << event << std::endl;
      artpc::ArTPCData_sptr arTPCData(new artpc::ArTPCData());
      Reader.setEvent(event);
      Reader.setArTPCData(arTPCData);

      artpc::EventAnalysisManager eventAnalysisManager(arTPCData);
      eventAnalysisManager.run();
      artpc::ArDetectorHit_sptr arDetectorHit(new artpc::ArDetectorHit());
      eventAnalysisManager.passDetectorHit(arDetectorHit);
      arHitArray.push_back(arDetectorHit);
    }
  }

  artpc::ArTreeIO treeIO;
  treeIO.defineBranches();
  treeIO.fillHits(arHitArray);

  time(&time2);
  std::cout << "Elapsed time: " << time2 - time1 << "second" << std::endl;

  return 0;
}
