#include <string>
/* CHANGE THESE OPTIONS */

const std::string CurveFilename = "test.txt";      //if FROM_SURFACE_FILE assumed input filename format of "XXXXX.stl"

// OPTION - how long should it run, when do you want data printed, what time value should it start at 
const double TTime = 100;       //total time of simulation (simulation units)
const double starttime = 0;       //total time of simulation (simulation units)
const double CurvePrintTime = 1; // print out the knot , without the velocity
// OPTION - what grid values do you want/ timestep
// timestep
const double dtime = 0.000001;         //size of each time step
