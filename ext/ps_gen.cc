//#include "head.h"

/* function generating the phasespace of e^-e^+->tau^+tau^- for vegas
 * with input x being 2 "random" numbers from the interval [0,1].
 * returns the respective jacobian and the 4 momenta for the outgoing particles. */
void PhSp_TaTa(double E12, double x[2], double intVars[2], double pOut[2][4], double *jacobian)
{
    int i;
    double sprime = E12*E12, // shat at this scale
           gamma = sqrt((E12-2*MT)*(E12+2*MT)), // gamma term
           absp = gamma/2,
           //E12*gamma/2,// absolute value of 3-momentum of final state taus
//            flux = 2*sprime*twoPi2, //flux factor Byckling 2.3
           phi = 2*M_PI*x[0], sin_phi = sin(phi), cos_phi = cos(phi),
//           theta = M_PI*x[2], sin_theta = sin(theta), cos_theta = cos(theta),
           cos_theta = 2*x[1] - 1, sin_theta = sqrt((1-cos_theta)*(1+cos_theta));

           
    // write the integration variables into an array to be used elsewhere
    intVars[0] = phi; intVars[1] = cos_theta; 
    
    // This consists of the Fluxfactor (III.2.3) and the R_2 (IV.1.8) from Byckling.
    *jacobian = gamma/16/sprime/E12/M_PI;
    
    // calculate the outgoing momenta for the process
    // tau^+ momentum
    pOut[0][0] = E12/2;
    pOut[0][1] = absp*sin_phi*sin_theta;
    pOut[0][2] = absp*cos_phi*sin_theta;
    pOut[0][3] = absp*cos_theta;
    
    // tau^- momentum
    pOut[1][0] = pOut[0][0]; 
    for (i=1;i<4;i++) pOut[1][i] = -pOut[0][i];
    
    // calculate the 2-particle phase space measure dPS2(...)
    
    // version with dcos(theta)
//     *jacobian = absp/flux
    // old vers *jacobian = absp/8/M_PI/E12/E12/E12; // Peskin-Schroeder 5.12 * 2*Pi * 2
}

// integrate the bare volume of a process e-e+->tau-tau+ using the generator
void PS2(double x[2], double *weight, double f[1]) 
{
//     double pIn[2][4] = { { Etot/2, 0.0, 0.0, Etot/2 }, { Etot/2, 0.0, 0.0, -Etot/2 } };
    double vars[2], jacobian, pOut[2][4];
           
//  just the volume 
    PhSp_TaTa(Etot,x,vars,pOut,&jacobian);
    
	f[0] = jacobian;

	return;
}

/* function calculating the phasespace measure of e^-e^+->tau^+tau^-photon for vegas
 * with input x being 5 "random" numbers from the interval [0,1].
 * includes the respective jacobians and gives the 4 momenta for the outgoing particles. */
void PhSp_TaTaPh(double E12, double x[5], double intVars[5], double pOut[3][4], double *jacobian)
{
    int i;
    double sprime = E12*E12, // shat at this energy
           flux = 2*sprime*twoPi5, //flux factor Byckling 2.3
           gamma = (E12-2*MT)*(E12+2*MT),
           s1 = x[0]*gamma + 4*MT2, //composite particle mass.
           phi1 = 2*M_PI*x[1], sin_phi1 = sin(phi1), cos_phi1 = cos(phi1),
//           theta1 = M_PI*x[2], sin_theta1 = sin(theta1), cos_theta1 = cos(theta1),
           cos_theta1 = 2*x[2] - 1, sin_theta1 = sqrt((1-cos_theta1)*(1+cos_theta1)),
           phi2 = 2*M_PI*x[3], sin_phi2 = sin(phi2), cos_phi2 = cos(phi2),
//           theta2 = M_PI*x[4], sin_theta2 = sin(theta2), cos_theta2 = cos(theta2),
           cos_theta2 = 2*x[4] - 1, sin_theta2 = sqrt((1-cos_theta2)*(1+cos_theta2)),
           E34 = (sprime+s1)/2/E12, // Energy of the p34 "particle" in cms
           absk = (sprime*sprime-s1*s1)/2/E12/(sprime+s1), // photon's 3 mom.
           e34 = sqrt(s1), // Energy of p34 in its rest frame
           absptil = sqrt(s1/4 - MT2); // abs. of 3-momentum in p34's rest frame

    // write the integration variables into an array to be used elsewhere
    intVars[0] = s1; intVars[1] = phi1; intVars[2] = cos_theta1; 
    intVars[3] = phi2; intVars[4] = cos_theta2; 
    
    // This consists of the Fluxfactor (III.2.3), R_3 and R_2 (IV.1.8) from Byckling.
    *jacobian = gamma * twoPi2 /** 4*/ // from mapping to [0..1]
                *(sprime-s1)/8/sprime // kallen contributions
                /**2*/*absptil/e34/*/8*/ // = sqrt(kallen(s1,MT2,MT2))/8/s1
                /flux; //divide by flux factor
           
    // calculate the outgoing momenta for the process
    // composite particle momentum in cms
    double p34[4] = { 
        E34, 
        absk*sin_phi1*sin_theta1, 
        absk*cos_phi1*sin_theta1, 
        absk*cos_theta1
    };
    // the boost of p3 and p4 to cms consists of p34/2 plus/minus this second part.
    double ptilde[4] = { 
        absk*absptil*cos_theta2/e34,
        
        absptil*(cos_phi1*sin_phi2*sin_theta2
        +sin_phi1*sin_theta1*cos_theta2*E34/e34
        +sin_phi1*cos_phi2*cos_theta1*sin_theta2),
        
        absptil*(-sin_phi1*sin_phi2*sin_theta2
        +cos_phi1*sin_theta1*cos_theta2*E34/e34
        +cos_phi1*cos_phi2*cos_theta1*sin_theta2),
        
        absptil*(cos_theta1*cos_theta2*E34/e34
        -cos_phi2*sin_theta1*sin_theta2)
    };
    
    for (i=0; i<4; i++)
    {
        pOut[0][i] = p34[i]/2 + ptilde[i]; // tau^+ momentum
        pOut[1][i] = p34[i]/2 - ptilde[i]; // tau^- momentum
        pOut[2][i] = -p34[i]; // photon momentum
    }
    pOut[2][0] = absk; // photon mass needs to vanish.
    
//     double CPS3 = gamma/128/M_PI/M_PI/M_PI/E12/E12; // constant for 3 particle phase space
    
    // calculate the 3-particle phase space measure dPS3(s1,...)
//     return sin_theta1*sin_theta2*sqrt(1-4*MT*MT/s1)*(E12*E12-s1)*
//     (E12*E12-4*MT*MT)/512/M_PI/E12/E12; version with dtheta
//     *jacobian = sqrt(1-4*MT*MT/s1)*(E12*E12-s1)*CPS3; // version with dcos(theta)
              
}

// Integrate the bare volume of a process e-e+->tau-tau+gamma using the generator
void PS3(double x[5], double *weight, double f[1]) 
{
//     double pIn[2][4] = { { Etot/2, 0.0, 0.0, Etot/2 }, { Etot/2, 0.0, 0.0, -Etot/2 } };
    double vars[5], jacobian, pOut[3][4];
           
//  just the volume 
    PhSp_TaTaPh(Etot,x,vars,pOut,&jacobian);
    
	f[0] = jacobian;

	return;
}
