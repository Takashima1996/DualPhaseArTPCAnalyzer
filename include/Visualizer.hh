#ifndef COMPTONSOFT_Visualizer_H
#define COMPTONSOFT_Visualizer_H

#include "PulseSignalArray.hh"
#include "S2FitParHolder.hh"

namespace artpc
{

class Visualizer
{
public:
    Visualizer(){};
    ~Visualizer(){};

    void plotFit(PulseSignalArray&, S2FitParHolder&);

private:
};

}   /* namespace artpc */

#endif
