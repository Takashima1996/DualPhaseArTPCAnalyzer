#ifndef COMPTONSOFT_S2FUNCTION_H
#define COMPTONSOFT_S2FUNCTION_H

#include <cmath>
#include <vector>
#include <boost/math/special_functions/erf.hpp>

namespace artpc
{

class S2Function
{
public:
    S2Function(){};
    S2Function(const std::vector<double>& parameter);
    ~S2Function(){};
    double operator()(double x) const;
    double operator()(double* x, double* par) const;

    std::vector<double> Parameter()  {return funcParameter_;}

private:
    std::vector<double>funcParameter_;
    double s2Num_ = 0;
    const int numPar_ = 9;
};

}   /* namespace artpc */

#endif
