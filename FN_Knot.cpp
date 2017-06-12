/* Fitzhugh-Nagumo reaction diffusion simulation with arbitrary vortex lines
   OPENMP VERSION
   Created by Carl Whitfield
   Last modified 03/01/17

   Operational order of the code:
   1) The code takes as input an stl file (defined in knot_filename) which defines an orientable surface with a boundary.
   2) This surface is scaled to fill a box of size xmax x ymax x zmax.
   3) A nunmerical integral is performed to calculate a phase field (phi_calc) on the 3D grid which winds around the boundary of the surface.
   4) This is then used to initialise the Fitzhugh-Nagumo set of partial differential equations such that on initialisation (uv_initialise):
   u = 2cos(phi) - 0.4        v = sin(phi) - 0.4
   The pde's used are
   dudt = (u - u^3/3 - v)/epsilon + Del^2 u
   dvdt = epsilon*(u + beta - gam v)
   5) The update method is Runge-Kutta fourth order (update_uv) unless RK4 is set to 0, otherwise Euler forward method is used.
   6) A parametric curve for the knot is found at each unit T


   The parameters epsilon, beta and gam are set to give rise to scroll waves (see Sutcliffe, Winfree, PRE 2003 and Maucher Sutcliffe PRL 2016) which eminate from a closed curve, on initialisation this curve corresponds to the boundary of the surface in the stl file.

   See below for various options to start the code from previous outputted data.*/
#include "FN_Knot.h"    //contains user defined variables for the simulation, and the parameters used 
#include <omp.h>
#include <math.h>
#include <string.h>

int main (void)
{
    Curve Curve; 

    cout << "Reading input file...\n";
    InitialiseCurve(CurveFilename,Curve);

    // initialising timers
    time_t then = time(NULL);
    time_t rawtime;
    time (&rawtime);
    struct tm * timeinfo;
    /*
       while(CurrentTime <= TTime)
       {
       Update();
       iterationcounter++;
       CurrentTime  = starttime + ((double)(iterationcounter) * dtime);

    // at this point, let people know how things are going
    cout << "T = " << CurrentTime << endl;
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    cout << "current time \t" << asctime(timeinfo) << "\n";
    }
    */
    return 0;
}

/*************************Functions for knot initialisation*****************************/


/*
   void Update(curve& Curve)
   {
   NP = knotcurves[c].knotcurve.size();  //store number of points in knot curve
   double dxds, dyds, dzds, dxdm, dydm, dzdm,bx, by, bz;
   totlength = 0;
   double T[3][3];
   double N[2][3];
   double B[3];
   double deltas[3]; double ds;
   double curvature[2];double torsion;
   for(s=0; s<NP; s++)    //fwd diff (defined on connecting line) (cell data in paraview)
   {
   for(i=0;i<3;i++)
   {

   dx = (knotcurves[c].knotcurve[incp(s,i+1,NP)].xcoord - knotcurves[c].knotcurve[incp(s,i,NP)].xcoord);   //central diff as a is defined on the points
   dy = (knotcurves[c].knotcurve[incp(s,i+1,NP)].ycoord - knotcurves[c].knotcurve[incp(s,i,NP)].ycoord);
   dz = (knotcurves[c].knotcurve[incp(s,i+1,NP)].zcoord - knotcurves[c].knotcurve[incp(s,i,NP)].zcoord);
   deltas[i] = sqrt(dx*dx+dy*dy+dz*dz);
   T[i][0] = dx/(deltas[i]);
   T[i][1] = dy/(deltas[i]);
   T[i][2] = dz/(deltas[i]);
   if(i==0)
   {
   knotcurves[c].knotcurve[s].length = deltas[0];
   dxds = T[0][0];
   dyds = T[0][1];
   dzds = T[0][2];
   ds = deltas[0];
   }
   }
   for(i=0;i<2;i++)
   {
   N[i][0] = (T[i+1][0]-T[i][0])/deltas[i];
   N[i][1] = (T[i+1][1]-T[i][1])/deltas[i];
   N[i][2] = (T[i+1][2]-T[i][2])/deltas[i];
   curvature[i] = sqrt(N[i][0]*N[i][0]+N[i][1]*N[i][1]+N[i][2]*N[i][2]);
   N[i][0] /=curvature[i];
   N[i][1] /=curvature[i];
   N[i][2] /=curvature[i];
   }
//compute the binormal with a cross product 
B[0] = T[0][1]*N[0][2] - N[0][1]*T[0][2] ;
B[1] = T[0][2]*N[0][0] - N[0][2]*T[0][0] ;
B[2] = T[0][0]*N[0][1] - N[0][0]*T[0][1] ;
// this is dn/ds +kt
double vx= (N[1][0]-N[0][0])/deltas[0] +curvature[0]*T[0][0];
double vy= (N[1][1]-N[0][1])/deltas[0]+ curvature[0]*T[0][1];
double vz= (N[1][2]-N[0][2])/deltas[0] + curvature[0]*T[0][2];
// lets get the torsion by computing the x component of the binomral, and looking at the scale factor between it and dn/ds+kt
double torsion = sqrt(vx*vx +vy*vy +vz*vz);
// we have lost the sign of the torsion here. to get it , compare the signs of dn/dx+kt and the binromal
double sign = 1;
(vx*B[0] + vy*B[1] + vz*B[2] >0)? sign = 1: sign = -1;
torsion *=sign;

knotcurves[c].knotcurve[s].curvature = curvature[0];
knotcurves[c].knotcurve[s].torsion = torsion ;
knotcurves[c].knotcurve[s].tx = T[0][0] ;
knotcurves[c].knotcurve[s].ty = T[0][1] ;
knotcurves[c].knotcurve[s].tz = T[0][2] ;
knotcurves[c].knotcurve[s].nx = N[0][0] ;
knotcurves[c].knotcurve[s].ny = N[0][1] ;
knotcurves[c].knotcurve[s].nz = N[0][2] ;
knotcurves[c].knotcurve[s].bx = B[0] ;
knotcurves[c].knotcurve[s].by = B[1] ;
knotcurves[c].knotcurve[s].bz = B[2] ;

knotcurves[c].length += knotcurves[c].knotcurve[s].length;
knotcurves[c].xavgpos += knotcurves[c].knotcurve[s].xcoord/NP;
knotcurves[c].yavgpos += knotcurves[c].knotcurve[s].ycoord/NP;
knotcurves[c].zavgpos += knotcurves[c].knotcurve[s].zcoord/NP;
}

}
*/

void print_knot( double t, vector<knotcurve>& knotcurves,const griddata& griddata)
{
    for( int c=0; c < (knotcurves.size()) ; c++)
    {

        /***Write values to file*******/
        stringstream ss;
        ss << "globaldata" << "_" << c <<  ".txt";
        ofstream wrout (ss.str().c_str(), std::ofstream::app);
        wrout << t << '\t' << knotcurves[c].writhe << '\t' << knotcurves[c].twist << '\t' << knotcurves[c].length << '\n';
        wrout.close();

        ss.str("");
        ss.clear();       

        ss << "knotplot" << c << "_" << t <<  ".vtk";
        ofstream knotout (ss.str().c_str());

        int i;
        int n = knotcurves[c].knotcurve.size();

        knotout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET UNSTRUCTURED_GRID\n";
        knotout << "POINTS " << n << " float\n";

        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].xcoord << ' ' << knotcurves[c].knotcurve[i].ycoord << ' ' << knotcurves[c].knotcurve[i].zcoord << '\n';
        }

        knotout << "\n\nCELLS " << n << ' ' << 3*n << '\n';

        for(i=0; i<n; i++)
        {
            knotout << 2 << ' ' << i << ' ' << incp(i,1,n) << '\n';
        }

        knotout << "\n\nCELL_TYPES " << n << '\n';

        for(i=0; i<n; i++)
        {
            knotout << "3\n";
        }

        knotout << "\n\nPOINT_DATA " << n << "\n\n";

        knotout << "\nSCALARS Curvature float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].curvature << '\n'; }

        knotout << "\nSCALARS Torsion float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].torsion << '\n';
        }

        knotout << "\nVECTORS A float\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].ax << ' ' << knotcurves[c].knotcurve[i].ay << ' ' << knotcurves[c].knotcurve[i].az << '\n';
        }

        knotout << "\nVECTORS V float\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].vx << ' ' << knotcurves[c].knotcurve[i].vy << ' ' << knotcurves[c].knotcurve[i].vz << '\n';
        }
        knotout << "\nVECTORS t float\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].tx << ' ' << knotcurves[c].knotcurve[i].ty << ' ' << knotcurves[c].knotcurve[i].tz << '\n';
        }
        knotout << "\nVECTORS n float\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].nx << ' ' << knotcurves[c].knotcurve[i].ny << ' ' << knotcurves[c].knotcurve[i].nz << '\n';
        }
        knotout << "\nVECTORS b float\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].bx << ' ' << knotcurves[c].knotcurve[i].by << ' ' << knotcurves[c].knotcurve[i].bz << '\n';
        }
        knotout << "\nVECTORS vdotn float\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].vdotnx << ' ' << knotcurves[c].knotcurve[i].vdotny << ' ' << knotcurves[c].knotcurve[i].vdotnz << '\n';
        }
        knotout << "\nVECTORS vdotb float\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].vdotbx << ' ' << knotcurves[c].knotcurve[i].vdotby << ' ' << knotcurves[c].knotcurve[i].vdotbz << '\n';
        }
        knotout << "\n\nCELL_DATA " << n << "\n\n";
        knotout << "\nSCALARS Writhe float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].writhe << '\n';
        }

        knotout << "\nSCALARS Twist float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].twist << '\n';
        }

        knotout << "\nSCALARS Length float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].length << '\n';
        }
        knotout.close();
    }
}

void cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *product)
{
    double p1 = gsl_vector_get(u, 1)*gsl_vector_get(v, 2)
        - gsl_vector_get(u, 2)*gsl_vector_get(v, 1);

    double p2 = gsl_vector_get(u, 2)*gsl_vector_get(v, 0)
        - gsl_vector_get(u, 0)*gsl_vector_get(v, 2);

    double p3 = gsl_vector_get(u, 0)*gsl_vector_get(v, 1)
        - gsl_vector_get(u, 1)*gsl_vector_get(v, 0);

    gsl_vector_set(product, 0, p1);
    gsl_vector_set(product, 1, p2);
    gsl_vector_set(product, 2, p3);
}
inline int circularmod(int i, int N)    // mod i by N in a cirucler fashion, ie wrapping around both in the +ve and -ve directions
{
    if(i<0) return N - ((-i)%N);
    else return i%N;
}
// inlined functions for incrementing things respecting boundaries
inline int incp(int i, int p, int n)    //increment i with p for periodic boundary
{
    if(i+p<0) return (N+i+p);
    else return ((i+p)%N);
}
inline int sign(int i)
{
    if(i==0) return 0;
    else return i/abs(i);
}

void InitialiseCurve(string CurveFilename,Curve& Curve)
{
    double xt,yt,zt;    //temporary variables

    ifstream CurveInputStream;   //knot file(s)
    string filename,buff;
    stringstream ss;

    ss << CurveFilename;
    filename = ss.str();
    CurveInputStream.open(filename.c_str());

    while(CurveInputStream.good())   //read in points for knot
    {
        if(getline(CurveInputStream,buff))
        {
            ss.clear();
            ss.str("");
            ss << buff;
            ss >> xt >> yt >> zt;
        }
        else break;
        // construct a point and put it on the curve
        CurvePoint Point;
        Point.xcoord = xt;
        Point.ycoord = yt;
        Point.zcoord = zt;
        Curve.Curve.push_back(Point);
    }
    Curve.N = Curve.Curve.size();
    Curve.InSync = 0;

    CurveInputStream.close();

    return;
}
