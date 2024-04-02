#include <iostream>
#include <iomanip>
#include <cmath>

#include <TFile.h>
#include <vector>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TMath.h>

#include <TGraphErrors.h>

#include <RooRealVar.h>
#include <RooChebychev.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooGlobalFunc.h>


//Purpose of this scripts:
//Do fitting for scatter plots inside fitCLG1_compare.root--- 



//linear fitting---
void LinearFit(TGraphErrors* graph, double minX, double maxX) {
    //Perform linear fitting with chi-squared method
    TF1 *fit = new TF1("linear_fit", "[0] + [1]*x", minX, maxX);
    graph->Fit(fit, "Q", "", minX, maxX);//"Q" for quiet mode and chi-squared method
    
    //Get the fitted parameters
    double slope = fit->GetParameter(1);
    double slopeError = fit->GetParError(1); // Uncertainty of slope parameter
    double intercept = fit->GetParameter(0);
    double interceptError = fit->GetParError(0); // Uncertainty of intercept parameter


    //fit goodness---
    double chiSquare = fit->GetChisquare();//chi^2
    int ndf = fit->GetNDF();//ndf
    double chiNDF = chiSquare / ndf;

    //Print the results
    
    std::cout << "\nk = " << std::setprecision(6) << slope << " ± " << slopeError 
              << ", b = " << std::setprecision(3) << intercept << " ± " << interceptError 
              << ",  chi^2/ndf = " << chiSquare << "/" << ndf << " = " << chiNDF << "\n" << std::endl;


    //Set the fitted curve color to match the scatter plot color
    fit->SetLineColor(graph->GetMarkerColor());

    //Draw the fitted line on the canvas
    fit->Draw("same");
}





//zoom in y and y error values of graph---
void scaleGraph(TGraphErrors* graph, double factor) {
    for (int i = 0; i < graph->GetN(); ++i) {
        double x_val, y_val, y_err;
        graph->GetPoint(i, x_val, y_val);
        y_err = graph->GetErrorY(i);
        // Multiply y-values and y-errors by the factor
        graph->SetPoint(i, x_val, factor * y_val);
        graph->SetPointError(i, factor * y_err);
    }
}




//Chebyshev Fitting---
void ChebyshevFit(TGraphErrors* graph, double minX, double maxX) {

//    TF1 *fit = new TF1("Chebyshev_fit", "[0]*1 + [1]*x + [2]*(2*x*x - 1)", minX, maxX);
    TF1 *fit = new TF1("Chebyshev_fit", "[0]*1 + [1]*x + [2]*(2*x*x - 1) + [3]*(4*x*x*x - 3*x)", minX, maxX);


    graph->Fit(fit, "Q", "", minX, maxX);//"Q" for quiet mode and chi-squared method
    
    //Retrieve the parameters
    double param0 = fit->GetParameter(0);
    double param1 = fit->GetParameter(1);
    double param2 = fit->GetParameter(2);
    double param3 = fit->GetParameter(3);

    //fit goodness---
    double chiSquare = fit->GetChisquare();//chi^2
    int ndf = fit->GetNDF();//ndf
    double chiNDF = chiSquare / ndf;

    //Print the results
    std::cout<< std::fixed;
//    std::cout << "p0 = " << std::setprecision(3) << param0 << ", p1 = " << std::setprecision(6) << param1 << ", p2 = " << std::setprecision(6) << param2 << ",   chi^2/ndf = " << chiSquare << "/" << ndf << " = " << chiNDF <<std::endl;

    std::cout << "p0 = " << std::setprecision(3) << param0 << ", p1 = " << std::setprecision(6) << param1 << ", p2 = " << std::setprecision(6) << param2 << ", p3 = " << std::setprecision(6) << param3 << ",   chi^2/ndf = " << chiNDF <<std::endl;

    //Set the fitted curve color to match the scatter plot color
    fit->SetLineColor(graph->GetMarkerColor());

    //Draw the fitted line on the canvas
    fit->Draw("same");
}









//===The MAIN function=========================================================
void compare_fit(){

    string file_path = "/Users/shuaixiangzhang/Work/current/FNAL_Work2024/rsl_analyses/v4_analysis/results/fit_Develop/cathode/";
    TFile *file = TFile::Open((file_path + "fitCLG1_compare.root").c_str(), "UPDATE");

    TCanvas *canvas = (TCanvas*)file->Get("fitBias");    

    //Delete added canvas due to last execution---
    gDirectory->cd(file->GetName());
    gDirectory->Delete("linear_fit;*");

    //clone the canvas---
    TCanvas *new_canvas = (TCanvas*)canvas->Clone();
    new_canvas->SetName("linear_fit");
    new_canvas->SetTitle("Plots with fitted lines");

    

    TIter next(canvas->GetListOfPrimitives());
    TObject *obj;
    while ((obj = next())) {
        // Check if the object is a TGraphErrors
        if (obj->InheritsFrom(TGraphErrors::Class())) {

            TGraphErrors* graph = (TGraphErrors*)obj;
            new_canvas->cd();

            LinearFit(graph, 50.0, 230.0);
//            ChebyshevFit(graph, 50.0, 200.0);
        }
    }


    new_canvas->Write();

    file->Close();


}

