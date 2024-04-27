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
    Point3D v1 = {302, -417, 149};
    Point3D v2 = {380, 588, 130};
    Point3D v3 = {-400, -588, 100};


    cout<<"\nAverage Solid angle: "<< distance_point_line(v1, v2, v3)<<endl;
    cout<<"\nAverage Solid angle: "<< solid_angle2(v1, v2, v3)<<endl;

}