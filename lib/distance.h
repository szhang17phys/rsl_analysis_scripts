#ifndef DISTANCE_H
#define DISTANCE_H
//This script is used to calculate the distance between 3D point and line segment---
//Author: Shuaixiang (Shu)---
//Date: Oct 8, 2023---

struct Point3D{
    double x, y, z;
};

double dotProduct(const Point3D& v1, const Point3D& v2){
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double magnitudeSquared(const Point3D& v){
    return v.x * v.x + v.y * v.y + v.z * v.z;
}


//distance between point and line=====================
double distance_point_line(const Point3D& point, const Point3D& segmentStart, const Point3D& segmentEnd){
    // Vector representing the line segment---
    Point3D segmentVector = {segmentEnd.x-segmentStart.x, segmentEnd.y-segmentStart.y, segmentEnd.z-segmentStart.z};

    //Vector from the start of the segment to the point---
    Point3D pointVector = {point.x-segmentStart.x, point.y-segmentStart.y, point.z-segmentStart.z};

    //Calculate the projection of pointVector onto segmentVector---
    double t = dotProduct(pointVector, segmentVector) / magnitudeSquared(segmentVector);

    //Check if the projection point is outside the line segment---
    if(t < 0.0){
        //The point is closest to the start of the segment---
        return std::sqrt(magnitudeSquared(pointVector));
    } 
    else if(t > 1.0){
        //The point is closest to the end of the segment---
        Point3D endToPoint = {point.x-segmentEnd.x, point.y-segmentEnd.y, point.z-segmentEnd.z};
        return std::sqrt(magnitudeSquared(endToPoint));
    } 
    else{
        //The point is closest to a point on the segment---
        Point3D projection = {
            segmentStart.x + t*segmentVector.x,
            segmentStart.y + t*segmentVector.y,
            segmentStart.z + t*segmentVector.z
        };
        Point3D pointToProjection = {point.x-projection.x, point.y-projection.y, point.z-projection.z};
        return std::sqrt(magnitudeSquared(pointToProjection));
    }

}





//output cos of angle of two vectors=======================
double cos(const Point3D& v1, const Point3D& v2){
    double dot = 0.0; //dot product
    double mag1 = 0.0; //length of v1
    double mag2 = 0.0;

    dot = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    mag1 = std::sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
    mag2 = std::sqrt(v2.x * v2.x + v2.y * v2.y + v2.z * v2.z);

    return dot/(1.0 * mag1 * mag2);
}


//output square of vector==================
double square(const Point3D& v1){
    return (v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
}




//distance, the solid angle this time===================================
//top is top crt, bot is bottom crt, opch is center of opch
double solid_angle(const Point3D& opch, const Point3D& top, const Point3D& bot){

    //specify points num from track-----------------------
    int num = 1024;
    //----------------------------------------------------

    //Store points of the track---
    Point3D points[num];
    double step = 1.0 / (num + 1);

    Point3D trueStart = {0.0, 0.0, 0.0};//near top crt---
    Point3D trueEnd = {0.0, 0.0, 0.0};

    trueStart.y = 417.6;
    trueEnd.y = -417.6;

    double ratio = (top.y - trueStart.y) / (top.y - bot.y);

    trueStart.x = top.x - ratio * (top.x - bot.x);
    trueStart.z = top.z - ratio * (top.z - bot.z);
    trueEnd.x = bot.x + ratio * (top.x - bot.x);
    trueEnd.z = bot.z + ratio * (top.z - bot.z);    

    for(int i=0; i<num; ++i){
        points[i].x = trueStart.x - (i+1) * step * (trueStart.x - trueEnd.x);
        points[i].y = trueStart.y - (i+1) * step * (trueStart.y - trueEnd.y);
        points[i].z = trueStart.z - (i+1) * step * (trueStart.z - trueEnd.z);
    }

    Point3D vectors[num];//from opch center to point---
    for(int i=0; i<num; ++i){
        vectors[i].x = points[i].x - opch.x;
        vectors[i].y = points[i].y - opch.y;
        vectors[i].z = points[i].z - opch.z;        
    }

    Point3D opchV = {0.0, 0.0, 0.0};//normal line of opch---
    if(opch.y>400 && opch.x>-300){//opch 0, 2, 16, 22
        opchV.x = 0.0;
        opchV.y = -1.0;
        opchV.z = 0.0;
    }
    if(opch.y<-400 && opch.x>-300){//opch 1, 3, 17, 23
        opchV.x = 0.0;
        opchV.y = 1.0;
        opchV.z = 0.0;
    }
    if(opch.x<50 && opch.x>-50){//cathode xa
        opchV.x = 1.0;
        opchV.y = 0.0;
        opchV.z = 0.0;
    }
    if(opch.z<-300){//ground pmts
        opchV.x = 1.0;
        opchV.y = 0.0;
        opchV.z = 0.0;
    }    
    if(opch.z<-50 && opch.x>-300){//-z high pmts
        opchV.x = 0.0;
        opchV.y = 0.0;
        opchV.z = 1.0;
    }
    if(opch.z>350 && opch.x>-300){//+z high pmts
        opchV.x = 0.0;
        opchV.y = 0.0;
        opchV.z = -1.0;
    }

    //Define solid angle---
    double omg =  0.0;
    double tmp = 0.0;
    for(int i=0; i<num; ++i){
        tmp = 3600 * cos(opchV, vectors[i]) / square(vectors[i]);
        if(opch.x<50 && opch.x>-50){//cathode xa double-side active
            tmp = std::abs(tmp);
        }
        omg += tmp;
    }
    omg = 1000.0 * omg / num;//make omg looking-better


    //test---
//    for(int i=0; i<num; ++i){
//        cout<<"Point "<<i<<": ("<<points[i].x<<", "<<points[i].y<<", "<<points[i].z<<")"<<endl;
//    }


    return omg;

}








#endif
