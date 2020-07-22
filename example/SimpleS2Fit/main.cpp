#include <TString.h>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <ctime>
#include <fstream>
#include <vector>

#include "ArTPCData.hh"
#include "ArTPCData_sptr.hh"
#include "DPArTPCEventAnalysis.hh"
#include "DataReader.hh"
#include "PulseSignalArray.hh"
#include "S2FitParHolder.hh"
#include "Visualizer.hh"

int main(int argc, char *argv[]) {
  time_t time1, time2;
  time(&time1);
  int runNum = 192;
  int fileNum = 1;

  int eventNum = boost::lexical_cast<int>(argv[1]);
  int s2Num = boost::lexical_cast<int>(argv[2]);

  TString Inputfile =
      (boost::format("../../../ANKOK_RUN-%1%_%2%.root") % runNum % fileNum)
          .str();

  artpc::ArTPCData_sptr arTPCData(new artpc::ArTPCData());

  artpc::DataReader Reader(Inputfile);
  Reader.setEvent(eventNum);
  Reader.setArTPCData(arTPCData);

  artpc::DPArTPCEventAnalysis analyzer(arTPCData);
  analyzer.setErrorMinos(true);
  analyzer.s2Fit(s2Num, true);

  artpc::S2FitParHolder parHolder;
  analyzer.passFitResult(parHolder);

  artpc::PulseSignalArray phaseSignalArray(arTPCData);
  phaseSignalArray.setPulseHeight();
  artpc::Visualizer visualizer;
  visualizer.plotFit(phaseSignalArray, parHolder);

  time(&time2);
  std::cout << "Elapsed time: " << time2 - time1 << "second" << std::endl;

  return 0;
}
