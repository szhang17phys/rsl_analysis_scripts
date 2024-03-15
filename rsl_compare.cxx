#include <iostream>

#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>


void readFile(const std::string& fileName, double distances[], double mpvConv[], double sigConv[], int maxArraySize) {
    //Open the input file for reading
    std::ifstream inputFile(fileName);

    //Check if the file is open
    if(!inputFile.is_open()) {
        std::cerr << "Error: Cannot open input file '" << fileName << "' for reading" << std::endl;
        return;
    }

    //Temporary variables to store data from each line
    double distance, mpv, sig;
    int index = 0;
    
    //Read data from each line and store them into arrays
    while (inputFile >> distance >> mpv >> sig && index < maxArraySize) {
        distances[index] = distance;
        mpvConv[index] = mpv;
        sigConv[index] = sig;
        index++;
    }

    //Close the input file
    inputFile.close();
}






//======MAIN FUNCTION=======================================
void rsl_compare(){

    const int num = 60; //num of slices---

    std::string inputFile99 = "../results/fit_Develop/rsl99_fit.txt";



    double dis99[num];
    double mpv99[num];
    double sig99[num];

    readFile(inputFile99, dis99, mpv99, sig99, num);

    //TEST------------------
    for(int i=0; i<60; ++i){
        std::cout<<"mpv: "<<mpv99[i]<<std::endl;
    }


}





