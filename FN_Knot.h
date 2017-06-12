//  FN_knot_code.h
//
//
//  Created by Carl Whitfield on 17/05/2016.
//
//  Last modified 3/11/16
#include "FN_Constants.h"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <time.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
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
    std::vector<Curvepoint> Curve; // the actual data of the curve
    int N;   // Number of points in the curve
    bool InSync;   // Says whether or not the geometric data (curvature torision, velocities) are in sync with the position data 
    double length;   //total lengthh of line
    double xavgpos;  // average position of the knot
    double yavgpos;
    double zavgpos;
};
// Main functions
void Update( vector<double>&ucvx, vector<double>&ucvy, vector<double>&ucvz, vector<double>& ucvmag,vector<double>&u,vector<knotcurve>& knotcurves,double t, gsl_multimin_fminimizer* minimizerstate, const griddata& griddata);
void PrintCurve( double t, vector<knotcurve>& knotcurves,const griddata& griddata);
void InitialiseCurve(string CurveFilename,Curve& Curve);
void cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *product);
// little inline guys
inline int sign(int i);
inline  int pt( int i,  int j,  int k,const griddata& griddata);       //convert i,j,k to single index
inline int circularmod(int i, int N);    // mod i by N in a cirucler fashion, ie wrapping around both in the +ve and -ve directions
inline int incp(int i, int p, int N);    //increment i with p for periodic boundary
inline int incw(int i, int p, int N);    //increment with reflecting boundary between -1 and 0 and N-1 and N
inline int gridinc(int i, int p, int N, int direction );    //increment with reflecting boundary between -1 and 0 and N-1 and N




