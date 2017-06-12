// Curve Dynamics 
// Created by Jack Binysh 12/06/2017

#include "CurveDynamics.h"    
#include <math.h>
#include <string.h>

int main (void)
{
    Curve Curve; 

    cout << "Reading input file...\n";
    InitialiseCurve(CurveFilename,Curve);

    int iterationcounter = (int) (starttime/dtime);
    int TotalIterations = (int) (TTime/dtime);
    int CurvePrintInterval = (int) (CurvePrintTime/dtime);
    while(iterationcounter <= TotalIterations)
    {
        UpdateState(Curve);
        if(iterationcounter % CurvePrintInterval == 0 )
        {
            double CurrentTime  = starttime + ((double)(iterationcounter) * dtime);
            PrintCurve(CurrentTime,Curve);
            cout << "T = " << CurrentTime << endl;
            time_t rawtime;
            time (&rawtime);
            struct tm * timeinfo = localtime (&rawtime);
            cout << "current time \t" << asctime(timeinfo) << "\n";
        }
        UpdatePositionFromVelocity(Curve);
        iterationcounter++;
    }

    return 0;
}

void UpdateState(Curve& Curve)
{
    UpdateGeometryFromPosition(Curve);
    UpdateVelocityFromGeometry(Curve);
}

void UpdateGeometryFromPosition(Curve& Curve)
{
    int NP = Curve.Curve.size();  //store number of points in knot Curve
    double dxds, dyds, dzds,bx, by, bz;
    double T[3][3];
    double N[2][3];
    double B[3];
    double deltas[3]; double ds;
    double curvature[2];double torsion;
    for(int s=0; s<NP; s++)    //fwd diff (defined on connecting line) (cell data in paraview)
    {
        for(int i=0;i<3;i++)
        {

            double dx = (Curve.Curve[incp(s,i+1,NP)].xcoord - Curve.Curve[incp(s,i,NP)].xcoord);   //central diff as a is defined on the points
            double dy = (Curve.Curve[incp(s,i+1,NP)].ycoord - Curve.Curve[incp(s,i,NP)].ycoord);
            double dz = (Curve.Curve[incp(s,i+1,NP)].zcoord - Curve.Curve[incp(s,i,NP)].zcoord);
            deltas[i] = sqrt(dx*dx+dy*dy+dz*dz);
            T[i][0] = dx/(deltas[i]);
            T[i][1] = dy/(deltas[i]);
            T[i][2] = dz/(deltas[i]);
            if(i==0)
            {
                Curve.Curve[s].length = deltas[0];
                dxds = T[0][0];
                dyds = T[0][1];
                dzds = T[0][2];
                ds = deltas[0];
            }
        }
        for(int i=0;i<2;i++)
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

        Curve.Curve[s].curvature = curvature[0];
        Curve.Curve[s].torsion = torsion ;
        Curve.Curve[s].tx = T[0][0] ;
        Curve.Curve[s].ty = T[0][1] ;
        Curve.Curve[s].tz = T[0][2] ;
        Curve.Curve[s].nx = N[0][0] ;
        Curve.Curve[s].ny = N[0][1] ;
        Curve.Curve[s].nz = N[0][2] ;
        Curve.Curve[s].bx = B[0] ;
        Curve.Curve[s].by = B[1] ;
        Curve.Curve[s].bz = B[2] ;

        Curve.length += Curve.Curve[s].length;
        Curve.xavgpos += Curve.Curve[s].xcoord/NP;
        Curve.yavgpos += Curve.Curve[s].ycoord/NP;
        Curve.zavgpos += Curve.Curve[s].zcoord/NP;
    }

}

void UpdateVelocityFromGeometry(Curve& Curve)
{
    int NP = Curve.Curve.size();  //store number of points in knot Curve
    for(int s=0; s<NP; s++)    //fwd diff (defined on connecting line) (cell data in paraview)
    {

        Curve.Curve[s].vx = Curve.Curve[s].curvature * Curve.Curve[s].bx;
        Curve.Curve[s].vy = Curve.Curve[s].curvature * Curve.Curve[s].by;
        Curve.Curve[s].vz = Curve.Curve[s].curvature * Curve.Curve[s].bz;
        Curve.Curve[s].vdott = 0;
        Curve.Curve[s].vdotn = 0;
        Curve.Curve[s].vdotb = Curve.Curve[s].curvature;

    }
    // now the Curve object has updated geometry and velocity
    Curve.InSync = 1;
}

void UpdatePositionFromVelocity(Curve& Curve)
{
    int NP = Curve.Curve.size();  //store number of points in knot Curve
    for(int s=0; s<NP; s++)    //fwd diff (defined on connecting line) (cell data in paraview)
    {

        Curve.Curve[s].xcoord = Curve.Curve[s].xcoord + Curve.Curve[s].vx*dtime; 
        Curve.Curve[s].ycoord = Curve.Curve[s].ycoord + Curve.Curve[s].vy*dtime; 
        Curve.Curve[s].zcoord = Curve.Curve[s].zcoord + Curve.Curve[s].vz*dtime; 
    }
    // now the Curve object is out of sync
    Curve.InSync = 0;
}

void PrintCurve(double t, const Curve& Curve)
{
    stringstream ss;
    ss << "globaldata" <<  ".txt";
    ofstream globalout (ss.str().c_str(), std::ofstream::app);
    int printtime = (int)(floor(t+0.5));
    globalout << printtime << '\t' << Curve.xavgpos << '\t' << Curve.yavgpos << '\t' << Curve.zavgpos << '\t' << Curve.length << '\n';
    globalout.close();

    ss.str("");
    ss.clear();       

    ss << "Curve"<< "_" << printtime <<  ".vtk";
    ofstream CurveOutputStream (ss.str().c_str());

    int i;
    int n = Curve.Curve.size();

    CurveOutputStream << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    CurveOutputStream << "POINTS " << n << " float\n";

    for(i=0; i<n; i++)
    {
        CurveOutputStream << Curve.Curve[i].xcoord << ' ' << Curve.Curve[i].ycoord << ' ' << Curve.Curve[i].zcoord << '\n';
    }

    CurveOutputStream << "\n\nCELLS " << n << ' ' << 3*n << '\n';

    for(i=0; i<n; i++)
    {
        CurveOutputStream << 2 << ' ' << i << ' ' << incp(i,1,n) << '\n';
    }

    CurveOutputStream << "\n\nCELL_TYPES " << n << '\n';

    for(i=0; i<n; i++)
    {
        CurveOutputStream << "3\n";
    }



    CurveOutputStream << "\n\nPOINT_DATA " << n << "\n\n";




    CurveOutputStream << "\nSCALARS Curvature float\nLOOKUP_TABLE default\n";
    for(i=0; i<n; i++)
    {
        CurveOutputStream << Curve.Curve[i].curvature << '\n'; }

    CurveOutputStream << "\nSCALARS Torsion float\nLOOKUP_TABLE default\n";
    for(i=0; i<n; i++)
    {
        CurveOutputStream << Curve.Curve[i].torsion << '\n';
    }
    CurveOutputStream << "\nSCALARS vdotn float\nLOOKUP_TABLE default\n";
    for(i=0; i<n; i++)
    {
        CurveOutputStream << Curve.Curve[i].vdotn << '\n';
    }
    CurveOutputStream << "\nSCALARS vdott float\nLOOKUP_TABLE default\n";
    for(i=0; i<n; i++)
    {
        CurveOutputStream << Curve.Curve[i].vdott << '\n';
    }
    CurveOutputStream << "\nSCALARS vdotb float\nLOOKUP_TABLE default\n";
    for(i=0; i<n; i++)
    {
        CurveOutputStream << Curve.Curve[i].vdotb << '\n';
    }
    CurveOutputStream << "\nSCALARS Length float\nLOOKUP_TABLE default\n";
    for(i=0; i<n; i++)
    {
        CurveOutputStream << Curve.Curve[i].length << '\n';
    }





    CurveOutputStream << "\nVECTORS t float\n";
    for(i=0; i<n; i++)
    {
        CurveOutputStream << Curve.Curve[i].tx << ' ' << Curve.Curve[i].ty << ' ' << Curve.Curve[i].tz << '\n';
    }
    CurveOutputStream << "\nVECTORS n float\n";
    for(i=0; i<n; i++)
    {
        CurveOutputStream << Curve.Curve[i].nx << ' ' << Curve.Curve[i].ny << ' ' << Curve.Curve[i].nz << '\n';
    }
    CurveOutputStream << "\nVECTORS b float\n";
    for(i=0; i<n; i++)
    {
        CurveOutputStream << Curve.Curve[i].bx << ' ' << Curve.Curve[i].by << ' ' << Curve.Curve[i].bz << '\n';
    }
    CurveOutputStream << "\nVECTORS V float\n";
    for(i=0; i<n; i++)
    {
        CurveOutputStream << Curve.Curve[i].vx << ' ' << Curve.Curve[i].vy << ' ' << Curve.Curve[i].vz << '\n';
    }



    CurveOutputStream.close();
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
    Curve.InSync = 0;

    CurveInputStream.close();

    return;
}
// inlined functions for incrementing things respecting boundaries
inline int incp(int i, int p, int N)    //increment i with p for periodic boundary
{
    if(i+p<0) return (N+i+p);
    else return ((i+p)%N);
}
