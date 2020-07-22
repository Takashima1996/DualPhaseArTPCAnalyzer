#include <vector>
#include <memory>
#include <iostream>

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TSpline.h>

#include "S2Function.hh"
#include "PulseSignalArray.hh"
#include "S2FitParHolder.hh"
#include "Visualizer.hh"

namespace artpc
{

void Visualizer::plotFit(PulseSignalArray& phaseSignalArray, S2FitParHolder& s2FitParHolder)
{
    phaseSignalArray.setPulseHeight();
    phaseSignalArray.errorCalculate();
    std::vector<double> exData(phaseSignalArray.getPulseSignalArrayR());
    std::vector<double> exError(phaseSignalArray.ErrorR());
    std::vector<double> timeArray(phaseSignalArray.TimeSeriesR());
    int t_size = exData.size();
    std::vector<double> timeError(t_size, 0.0);
    int s2Number = s2FitParHolder.S2Number();
    std::vector<double>par;

    std::vector<double> A(s2FitParHolder.Amplitude());
    std::vector<double> T(s2FitParHolder.GasT());
    std::vector<double> t0(s2FitParHolder.DriftT());
    std::vector<double> sigma(s2FitParHolder.Sigma());
    std::vector<double> tau1(s2FitParHolder.Tau1());
    std::vector<double> tau2(s2FitParHolder.Tau2());
    std::vector<double> p(s2FitParHolder.S2RatioFast());
    
    par.push_back(s2FitParHolder.y0());
    for(int i=0; i<s2Number; i++){
        par.push_back(A.at(i));
        par.push_back(T.at(i));
        par.push_back(t0.at(i));
        par.push_back(sigma.at(i));
        par.push_back(tau1.at(i));
        par.push_back(p.at(i));
        par.push_back(tau2.at(i));
        par.push_back(0.0);
        par.push_back(0.0);
    }

    S2Function func(par);
    std::vector<double> PredictData;

    for(int i=0;i < t_size; i++){
        PredictData.push_back(func(timeArray.at(i)));
    }
    
    TCanvas canvas1("c1","",10,10,300,300);
    TGraphErrors tg1(t_size, &(timeArray.at(0)), &(exData.at(0)), 
                                &(timeError.at(0)), &(exError.at(0)));
    tg1.GetXaxis()->SetTitle("Time[#mus]");
    tg1.GetYaxis()->SetTitle("Photoelectron Number");
    tg1.GetXaxis()->SetLimits(-10.0,130.0);
    tg1.SetTitle("");
    tg1.SetMarkerColor(kBlue);
    tg1.Draw("APZ");
    TGraph tg2(t_size, &(timeArray.at(0)), &(PredictData.at(0)));
    TSpline3 sp("spline", &tg2);
    sp.SetLineColor(kRed);
    sp.Draw("same");

    canvas1.SaveAs("FitResult.pdf");
}

}   /* namespace artpc */
