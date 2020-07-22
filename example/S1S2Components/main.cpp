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

int main(int argc, char *argv[]) {
  time_t time1, time2;
  time(&time1);

  int runNum = 192;
  int fileNum = 1;
  int eventNum = boost::lexical_cast<int>(argv[1]);

  TString Inputfile =
      (boost::format("../../../ANKOK_RUN-%1%_%2%.root") % runNum % fileNum)
          .str();

  artpc::ArTPCData_sptr arTPCData(new artpc::ArTPCData());
  artpc::DataReader Reader(Inputfile);
  Reader.setEvent(eventNum);
  Reader.setArTPCData(arTPCData);

  // Analyze part
  artpc::DPArTPCEventAnalysis analyzer(arTPCData);
  analyzer.s1s2Integral();

  std::vector<double> PulseComponents = analyzer.S1S2Components();
  std::cout << "S1Fast: " << PulseComponents.at(0) << std::endl;
  std::cout << "S1Slow: " << PulseComponents.at(1) << std::endl;
  std::cout << "S2Total: " << PulseComponents.at(2) << std::endl;

  time(&time2);
  std::cout << "Elapsed time: " << time2 - time1 << "second" << std::endl;

  return 0;
}
