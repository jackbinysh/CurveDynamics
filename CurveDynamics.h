#include "CurveDynamics_Constants.h"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <time.h>
using namespace std;

struct CurvePoint
{
    double xcoord;   //position vector x coord
    double ycoord;   //position vector y coord
    double zcoord;   //position vector z coord
    double tx;       //grad vector x coord
    double ty;       //grad vector y coord
    double tz;       //grad vector z coord
    double nx;       //grad vector x coord
    double ny;       //grad vector y coord
    double nz;       //grad vector z coord
    double bx;       //grad vector x coord
    double by;       //grad vector y coord
    double bz;       //grad vector z coord
    double vx;       //grad vector x coord
    double vy;       //grad vector y coord
    double vz;       //grad vector z coord
    double vdotn;       //grad vector y coord
    double vdotb;       //grad vector x coord
    double vdott;       //grad vector x coord
    double curvature;        // curvature
    double torsion;        // torsion
    double length;   //length of line
};

struct Curve
{
    std::vector<CurvePoint> Curve; // the actual data of the curve
    bool InSync;   // Says whether or not the geometric data (curvature torision, velocities) are in sync with the position data 
    double length;   //total lengthh of line
    double xavgpos;  // average position of the knot
    double yavgpos;
    double zavgpos;
};
// Main functions
void UpdateState(Curve& Curve);
void UpdateVelocityFromGeometry(Curve& Curve);
void UpdatePositionFromVelocity(Curve& Curve);
void UpdateGeometryFromPosition(Curve& Curve);
void PrintCurve( double t, const Curve& Curve);
void InitialiseCurve(string CurveFilename,Curve& Curve);
// little inline guys
inline int incp(int i, int p, int N);    //increment i with p for periodic boundary






