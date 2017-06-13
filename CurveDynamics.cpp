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
    // reset the global accumulated quantities
    Curve.length = 0;
    int NP = Curve.Curve.size();  //store number of points in knot Curve
    for(int s=0; s<NP; s++)    //fwd diff (defined on connecting line) (cell data in paraview)
    {
        // forward difference on the tangents
        double dx = (Curve.Curve[incp(s,1,NP)].xcoord - Curve.Curve[incp(s,0,NP)].xcoord);   //central diff as a is defined on the points
        double dy = (Curve.Curve[incp(s,1,NP)].ycoord - Curve.Curve[incp(s,0,NP)].ycoord);
        double dz = (Curve.Curve[incp(s,1,NP)].zcoord - Curve.Curve[incp(s,0,NP)].zcoord);
        double deltas = sqrt(dx*dx+dy*dy+dz*dz);
        Curve.Curve[s].tx = dx/(deltas);
        Curve.Curve[s].ty = dy/(deltas);
        Curve.Curve[s].tz = dz/(deltas);
        Curve.Curve[s].length = deltas;
        Curve.length +=deltas;
    }
    for(int s=0; s<NP; s++)    //fwd diff (defined on connecting line) (cell data in paraview)
    {
        double nx = 2.0*(Curve.Curve[s].tx-Curve.Curve[incp(s,-1,NP)].tx)/(Curve.Curve[s].length+Curve.Curve[incp(s,-1,NP)].length);
        double ny = 2.0*(Curve.Curve[s].ty-Curve.Curve[incp(s,-1,NP)].ty)/(Curve.Curve[s].length+Curve.Curve[incp(s,-1,NP)].length);
        double nz = 2.0*(Curve.Curve[s].tz-Curve.Curve[incp(s,-1,NP)].tz)/(Curve.Curve[s].length+Curve.Curve[incp(s,-1,NP)].length);
        double curvature = sqrt(nx*nx+ny*ny+nz*nz);
        nx /=curvature;
        ny /=curvature;
        nz /=curvature;
        double tx = Curve.Curve[s].tx ;
        double ty =  Curve.Curve[s].ty ;
        double tz = Curve.Curve[s].tz ;
        double bx = ty*nz - tz*ny;
        double by = tz*nx - tx*nz;
        double bz = tx*ny - ty*nx;
        Curve.Curve[s].nx = nx ;
        Curve.Curve[s].ny = ny ;
        Curve.Curve[s].nz = nz ;
        Curve.Curve[s].bx = bx ;
        Curve.Curve[s].by = by ;
        Curve.Curve[s].bz = bz ;
        Curve.Curve[s].curvature = curvature ;
    }
    // torsions with a central difference
    for(int s=0; s<NP; s++)   
    {
        double bx = Curve.Curve[s].bx;
        double by =  Curve.Curve[s].by;
        double bz = Curve.Curve[s].bz;

        double dnxds = 2.0*(Curve.Curve[incp(s,1,NP)].nx-Curve.Curve[incp(s,-1,NP)].nx)/(Curve.Curve[incp(s,1,NP)].length+Curve.Curve[incp(s,-1,NP)].length);
        double dnyds = 2.0*(Curve.Curve[incp(s,1,NP)].ny-Curve.Curve[incp(s,-1,NP)].ny)/(Curve.Curve[incp(s,1,NP)].length+Curve.Curve[incp(s,-1,NP)].length);
        double dnzds = 2.0*(Curve.Curve[incp(s,1,NP)].nz-Curve.Curve[incp(s,-1,NP)].nz)/(Curve.Curve[incp(s,1,NP)].length+Curve.Curve[incp(s,-1,NP)].length);

        double torsion = bx*dnxds+by*dnyds+bz*dnzds;
        Curve.Curve[s].torsion = torsion ;
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

#ifdef SMOOTH
    // at this point, we apply the 5 point smoothing of Longuett-Higgins and Cokelet
    static std::vector<double> x(NP); 
    static std::vector<double> y(NP); 
    static std::vector<double> z(NP); 
    for(int s=0; s<NP; s++)    //fwd diff (defined on connecting line) (cell data in paraview)
    {
        x[s] = (1/(double)(16))*(-Curve.Curve[incp(s,-2,NP)].xcoord+4*Curve.Curve[incp(s,-1,NP)].xcoord+10*Curve.Curve[incp(s,0,NP)].xcoord+4*Curve.Curve[incp(s,1,NP)].xcoord-Curve.Curve[incp(s,2,NP)].xcoord);
        y[s] = (1/(double)(16))*(-Curve.Curve[incp(s,-2,NP)].ycoord+4*Curve.Curve[incp(s,-1,NP)].ycoord+10*Curve.Curve[incp(s,0,NP)].ycoord+4*Curve.Curve[incp(s,1,NP)].ycoord-Curve.Curve[incp(s,2,NP)].ycoord);
        z[s] = (1/(double)(16))*(-Curve.Curve[incp(s,-2,NP)].zcoord+4*Curve.Curve[incp(s,-1,NP)].zcoord+10*Curve.Curve[incp(s,0,NP)].zcoord+4*Curve.Curve[incp(s,1,NP)].zcoord-Curve.Curve[incp(s,2,NP)].zcoord);
    }
    
    static unsigned int counter = 0;
    if(counter%(int)(SmoothingTimescale/dtime)==0)
    {
        for(int s=0; s<NP; s++)    //fwd diff (defined on connecting line) (cell data in paraview)
        {
            Curve.Curve[s].xcoord = x[s]; 
            Curve.Curve[s].ycoord = y[s]; 
            Curve.Curve[s].zcoord = z[s]; 
        }
        counter = 0;
    }
    counter ++ ;
#endif // SMOOTH
    // now the Curve object is out of sync
    Curve.InSync = 0;
}

void PrintCurve(double t, const Curve& Curve)
{
    stringstream ss;
    ss << "globaldata" <<  ".txt";
    ofstream globalout (ss.str().c_str(), std::ofstream::app);
    int printtime = (int)(floor(t+0.5));
    globalout << printtime << '\t'  << Curve.length << '\n';
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
