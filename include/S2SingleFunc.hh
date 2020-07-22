#ifndef COMPTONSOFT_S2SingleFunc_H
#define COMPTONSOFT_S2SingleFunc_H

#include <vector>
#include <cmath>
#include <boost/math/special_functions/erf.hpp>

namespace artpc
{

double S2SingleFunc(double x, std::vector<double>& par)
{
      double f=
      (par[0] / (2.0 * par[1])) *
      ((boost::math::erf((x - par[2]) / (std::sqrt(2.0) * par[3])) -
        std::exp(-(x - par[2]) / par[4] +
            par[3] * par[3] / (2.0 * par[4] * par[4]) +
            std::log(boost::math::erfc((par[3] * par[3] - par[4] * (x - par[2])) /
                               (std::sqrt(2) * par[3] * par[4])))) -
        boost::math::erf((x - par[1] - par[2]) / (std::sqrt(2.0) * par[3])) +
        std::exp(-(x - par[1] - par[2]) / par[4] +
                 par[3] * par[3] / (2.0 * par[4] * par[4]) +
                 std::log(boost::math::erfc(
                     (par[3] * par[3] - par[4] * (x - par[1] - par[2])) /
                     (std::sqrt(2) * par[3] * par[4]))))) *
           par[5] +
       (boost::math::erf((x - par[2]) / (std::sqrt(2.0) * par[3])) -
        std::exp(-(x - par[2]) / par[6] +
                 par[3] * par[3] / (2.0 * par[6] * par[6])) *
            boost::math::erfc((par[3] * par[3] - par[6] * (x - par[2])) /
                      (std::sqrt(2) * par[3] * par[6])) -
        boost::math::erf((x - par[1] - par[2]) / (std::sqrt(2.0) * par[3])) +
        std::exp(-(x - par[1] - par[2]) / par[6] +
                 par[3] * par[3] / (2.0 * par[6] * par[6])) *
            boost::math::erfc((par[3] * par[3] - par[6] * (x - par[2] - par[1])) /
                      (std::sqrt(2) * par[3] * par[6]))) *
           (1 - par[5]));

      return f;
}

} /* namespace artpc */

#endif
