/**************************************************************************************/
// Common Utilities
/**************************************************************************************/

#define PI 3.14159265358979

// Calculate Delta Phi
double dphi(double phi1,double phi2)
{
    double a=phi1-phi2;
    while (a<-PI) a+=2*PI;
    while (a>PI) a-=2*PI;
    
    if (a<-PI/2) a=2*PI+a;
    return a;
}
