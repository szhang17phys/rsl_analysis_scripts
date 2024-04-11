#include <iostream>
#include "../lib/XArapucaInfo.h"
#include "../lib/distance.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <string.h>
#include <stdio.h>
#include <vector>

#include <cstdlib>


using namespace std;

void test(){
    Point3D opch = {0, 0, 0};
    Point3D topCRT = {400, 588.2, 100};
    Point3D botCRT = {-300, -588.2, 180};

    m_solid_angle(opch, topCRT, botCRT);

}