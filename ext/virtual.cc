// Calculates the scalar onepoint function
void A0(std::vector<complex> &res, double m0)
{
    std::vector<double> m = {m0};
    std::vector<double> p = {};
    ql::QCDLoop<complex, double, double > qcdloop;
    
    qcdloop.integral(res,muR2,m,p);
    
    return;
}

//Calculates the scalar twopoint function
void B0(std::vector<complex> &res, double p1, double m0, double m1)
{
    std::vector<double> m = {m0,m1};
    std::vector<double> p = {p1};
    ql::QCDLoop<complex, double, double> qcdloop;
    
    qcdloop.integral(res,muR2,m,p);
    
    return;
}

//Calculates the scalar threepoint function
void C0(std::vector<complex> &res, double p1, double p2, double p3, 
           double m0, double m1, double m2)
{
    std::vector<double> m = {m0,m1,m2};
    std::vector<double> p = {p1,p2,p3};
    ql::QCDLoop<complex, double, double > qcdloop;
    
    qcdloop.integral(res,muR2,m,p);
    
    return;
}

//Calculates the scalar fourpoint function
void D0(std::vector<complex> &res, double p1, double p2, double p3, double p4,
        double p1p2, double p1p3, double m0, double m1, double m2, double m3)
{
    std::vector<double> m = {m0,m1,m2,m3};
    std::vector<double> p = {p1,p2,p3,p4,p1p2,p1p3};
    ql::QCDLoop<complex, double, double > qcdloop;
    
    qcdloop.integral(res,muR2,m,p);
    
    return;
}

//calculate the integral B0(p,M,M)
void B01(std::vector<complex> &res, double s, double m02, double m12)
{
    double //s = p*p, m2 = m0*m0
           beta = sqrt(1.0-4.0*m02/s),
           lambdap = (1.0+beta)/2.0,
           lambdam = (1.0-beta)/2.0;
           
    res[0] = 2 + log(muR2/m02) - log(lambdap/lambdam) - I*M_PI;
    res[1] = 0.0;
    res[2] = 0.0;
    res[3] = 1.0;  // set the divergence to be UV.
    
    return;
}

//calculate the integral B0(0,M,0)
void B02(std::vector<complex> &res, double s, double m02, double m12)
{
    //calculate the integral B0(0,MMU,MMU)
    
    res[0] = 1 + log(muR2/m02);
    res[1] = 0.0;
    res[2] = 0.0;
    res[3] = 1.0;  // set the divergence to be UV.
    
    return;
}

// calculates the terms that are phase space independent
void initVirtCor()
{   
    int i;
    
    //calculate the integral A0(MT)
    A0(A0MT,MT2);
    
    A0MT[3] = A0MT[1]; // set the divergence to be UV.
    A0MT[1] = 0.0;
    
    //calculate the integral A0(MMU)
    A0(A0MMU,MMU2);
    
    A0MMU[3] = A0MMU[1]; // set the divergence to be UV.
    A0MMU[1] = 0.0;
    
    //calculate the integral A0(ME)
    A0(A0ME,ME2);
    
    A0ME[3] = A0ME[1]; // set the divergence to be UV.
    A0ME[1] = 0.0;
    
    //set B0(0,0,0) = Duv - Dir.
    B00[0] = 0.0;
    B00[1] = -1.0;
    B00[2] = 0.0;
    B00[3] = 1.0;
    
    //calculate the integral B0(p1+p2,0,0)
    B0(B0s,shat,0.0,0.0);
    
    B0s[3] = B0s[1]; // set the divergence to be UV.
    B0s[1] = 0.0;
    
    //calculate the integral B0(MT,0,MT)
    B0(B0MT,MT2,0.0,MT2);
    
    B0MT[3] = B0MT[1]; // set the divergence to be UV.
    B0MT[1] = 0.0;
    
    //calculate the integral B0(0,MT,MT)
    B0(B00MT,0.0,MT2,MT2);
    
    B00MT[3] = B00MT[1]; // set the divergence to be UV.
    B00MT[1] = 0.0;
    
    //calculate the integral B0(0,MMU,MMU)
    B0(B00MMU,0.0,MMU2,MMU2);
    
    B00MMU[3] = B00MMU[1]; // set the divergence to be UV.
    B00MMU[1] = 0.0;
    
    //calculate the integral B0(0,ME,ME)
    B0(B00ME,0.0,ME2,ME2);
    
    B00ME[3] = B00ME[1]; // set the divergence to be UV.
    B00ME[1] = 0.0;
    
    //calculate the integral B0(0,MT,0)
    B02(B0MT0,0.0,MT2,0.0);
    
//     B0MT0[3] = B0MT0[1]; // set the divergence to be UV.
//     B0MT0[1] = 0.0;
    
    //calculate the integral B0(p3+p4,MT,MT)
    B01(B0sMT,shat,MT2,MT2);
    
//     B0sMT[3] = B0sMT[1]; // set the divergence to be UV.
//     B0sMT[1] = 0.0;
    
    //calculate the integral B0(p3+p4,MMU,MMU)
    B01(B0sMMU,shat,MMU2,MMU2);
    
//     B0sMMU[3] = B0sMMU[1]; // set the divergence to be UV.
//     B0sMMU[1] = 0.0; // B0 does usually not have an IR divergence
    
    //calculate the integral B0(p3+p4,ME,ME)
    B01(B0sME,shat,ME2,ME2);
    
    
    // for the fermion bubble this potentially needs to be calculated with different mfs!
    
    //calculate the integral C0(p1,p1+p2,0,0,0) for electron triangle and boxes
    C0(C0s,0.0,0.0,shat,0.0,0.0,0.0);
    
    C0s[3] = 0.0; // C0's do not have UV divergences.
    
    //calculate the integral C0(p3,p3+p4,MT,0,MT) for tau triangle
    C0(C0sMT,shat,MT2,MT2,MT2,MT2,0.0);
    
    C0sMT[3] = 0.0; // C0's do not have UV divergences.
    
    //calculate the integral C0(-p3/p4,-p1-p2,0,MT,0) for the boxes
    C0(C0MTs,MT2,MT2,shat,0.0,MT2,0.0);
    
    C0MTs[3] = 0.0; // C0's do not have UV divergences.
    
    //calculate the counter terms:
    
    for(i=0;i<4;i++) {
        //Electron Bubble:
        ElBCT[i] = B00[i];
        
        //massive Electron Bubble:
        mElBCT[i] = B00ME[i];
    
        //Muon Bubble:
        MuBCT[i] = B00MMU[i];
    
        //Tauon Bubble:
        TauBCT[i] = B00MT[i];
    }
    //ElBCT[1] = 0.0;// for testing so the thing is: i think i need to reintroduce electron mass for the bubble.
    
    //Electron Triangle: this needs some tweaking. what about the derivative.
    ElCT[0] = B00[0];// - 1.0;// + 4.0 + 2.0*logMuMT; //check if this is correct.
    ElCT[1] = B00[1];//2.0;
    ElCT[2] = B00[2];
    ElCT[3] = B00[3];
    
    //Tauon Triangle:
    TauCT[0] = B0MT0[0] - 1.0 + 4.0 + 2.0*logMuMT;
    TauCT[1] = 2.0;
    TauCT[2] = 0.0;
    TauCT[3] = B0MT0[3];
    
    
    return;
}


// calculates the polarized tree level matrix element from currents.
void currBornME(complex res[4][4],double p1p2,double p1p3, double &pref,
                complex j1[4][4], complex j2[4][4])
{
    int i,j;
    double p1p4 = p1p2 - p1p3; 
    
    pref = 1.0/(p1p2*(MT2*p1p2 - 2*p1p3*p1p4));
    
    // calculate the polarized matrix element
    for(i=0;i<4;i++) {
        res[0][i] = 0.0;
        
        // j1(1,1) part only nonzero for this configuration
        res[1][i] = e2*pref*j1[3][1]*(j2[1][i]*MT*p1p2 + j2[2][i]*p1p4 - j2[3][i]*p1p3);
//         printf("res[1][%d] =  % .10f %+.10fi\n",i,real(res[1][i]),imag(res[1][i]));
        
        // j1(1,-1) part only  nonzero for this configuration
        res[2][i] = e2*pref*j1[2][2]*(j2[0][i]*MT*p1p2 + j2[3][i]*p1p4 - j2[2][i]*p1p3);
//         printf("j11m1[2] =  % .10f %+.10fi\n",i,real(j11m1[2]),imag(j11m1[2]));
//         printf("res[2[%d] =  % .10f %+.10fi\n",i,real(res[2][i]),imag(res[2][i]));
        
        res[3][i] = 0.0;
    }
    
    return;
}


// calculates the most common non-factorizing rest matrix element from currents.
void currAntiBornME(complex res[4][4],double p1p2,double p1p3, double &pref,
                complex j1[4][4], complex j2[4][4])
{
    int i,j;
    double p1p4 = p1p2 - p1p3; 
    
//     pref = e2/(p1p2*(MT2*p1p2 - 2*p1p3*p1p4)); // its the same as born.
    
    // calculate the polarized matrix element
    for(i=0;i<4;i++) {
        res[0][i] = 0.0;
        
        // j1(1,1) part only nonzero for this configuration
        res[1][i] = e2*pref*j1[3][1]*(j2[0][i]*MT*p1p2 + j2[3][i]*p1p4 - j2[2][i]*p1p3);
//         printf("res[1][%d] =  % .10f %+.10fi\n",i,real(res[1][i]),imag(res[1][i]));
        
        // j1(1,-1) part only  nonzero for this configuration
        res[2][i] = e2*pref*j1[2][2]*(j2[1][i]*MT*p1p2 + j2[2][i]*p1p4 - j2[3][i]*p1p3);
//         printf("j11m1[2] =  % .10f %+.10fi\n",i,real(j11m1[2]),imag(j11m1[2]));
//         printf("res[2[%d] =  % .10f %+.10fi\n",i,real(res[2][i]),imag(res[2][i]));
        
        res[3][i] = 0.0;
    }
    
    return;
}



// calculates the rest that is neither born nor antiborn.
void currRestME(complex res[4][4], complex j1[4][4], complex j2[4][4])
{
    int i,j;
    
    for(i=0;i<4;i++)
        for(j=0;j<4;j++) res[i][j] = 0.0;
        
    // j1(1,1) part only nonzero for these configurations
    res[1][0] = j1[3][1]*(j2[0][0]+j2[1][0]);
//         printf("j1[3][1] =  % .16f %+.16fi\n",real(j1[3][1]),imag(j1[3][1]));
//         printf("j2[0][0] =  % .10f %+.10fi\n",real(j2[0][0]),imag(j2[0][0]));
//         printf("j2[1][0] =  % .10f %+.10fi\n",real(j2[1][0]),imag(j2[1][0]));
//         printf("res[1][0] =  % .10f %+.10fi\n",real(res[1][0]),imag(res[1][0]));
    res[1][3] = j1[3][1]*(j2[0][3]+j2[1][3]);
//         printf("sum      =  % .16f %+.16fi\n",real(j2[0][3]+j2[1][3]),
//                imag(j2[0][3]+j2[1][3]));
//         printf("j2[1][3] =  % .10f %+.10fi\n",real(j2[1][3]),imag(j2[1][3]));
//         printf("res[1][3] =  % .10f %+.10fi\n",real(res[1][3]),imag(res[1][3]));
        
    // j1(1,-1) part only  nonzero for this configuration
    res[2][0] = j1[2][2]*(j2[0][0]+j2[1][0]);
//         printf("res[2][0] =  % .10f %+.10fi\n",real(res[2][0]),imag(res[2][0]));
    res[2][3] = j1[2][2]*(j2[0][3]+j2[1][3]);
//         printf("res[2][3] =  % .10f %+.10fi\n",real(res[2][3]),imag(res[2][3]));
    
    return;
}

void getBubPref(complex bPref[4],double p1p2)
{
    int i,poles;
    double exp[4] = {1.0,0.0,0.0,0.0}; // finite epsilon expansion term.
    
    // So far only calculates for leptons let's see for the rest
    
    // decide whether poles are calculated
    if(cEps) poles = 4;
    else poles = 1;
    
    // electron part. IR divergent. Check if this gets cancelled
    for(i=0;i<poles;i++) {
        bPref[i] += e2/9/twoPi2*(
    // massless electron bubble
//                     + exp[i] // term from expansion in eps
//                     - B0s[i] * 3.0
//                     + ElBCT[i] * 3.0
    // massive electron bubble
                    + (exp[i] * (p1p2 - 3.0*ME2) // term from expansion in eps
                    + A0ME[i] * 3.0
                    - B0sME[i] * 3.0*(p1p2 + ME2)
                    + mElBCT[i] * 3.0*p1p2 // electron counterterm
                    )/p1p2
    // muon part
                    + (exp[i] * (p1p2 - 3.0*MMU2) // term from expansion in eps
                    + A0MMU[i] * 3.0
                    - B0sMMU[i] * 3.0*(p1p2 + MMU2)
                    + MuBCT[i] * 3.0*p1p2 // mu counterterm
                    )/p1p2
    
    // tauon part
                    + (exp[i] * (p1p2 - 3.0*MT2) // term from expansion in eps
                    + A0MT[i] * 3.0
                    - B0sMT[i] * 3.0*(p1p2 + MT2)
                    + TauBCT[i] * 3.0*p1p2// tau bubble counterterm
                    )/p1p2
                    );
        
//         printf("nPref[%d] =  % .16e %+.16ei\n",i,real(nPref[i]),imag(nPref[i]));
    }
    
    return;
}

// calculates the prefactors of the electron and tau triangle diagrams.
void getTriPref(complex bPref[4],complex rePref[4],double p1p2)
{
    int i,poles;
    double exp[4] = {1.0,0.0,0.0,0.0}; // Expansion term.
    
    // decide whether poles are calculated
    if(cEps) poles = 4;
    else poles = 1;
    
    for(i=0;i<poles;i++) {
        
        bPref[i] -= e2/4.0/twoPi2*(
            // Electron triangle part.
                    + exp[i] // term from expanding B0s prefactor in eps
                    - B00[i] * 4.0
                    + B0s[i] * 3.0
                    + ElCT[i]// counter term for electron triangle
                    + C0s[i] * 4.0*p1p2
                    
            // Tauon triangle part.
                    - B0MT[i] * 4.0
                    + exp[i] // term from expanding B0sMT prefactor in eps
                    + B0sMT[i] * 3.0
                    + TauCT[i] // counter term for tau triangle
                    - C0sMT[i] * 4.0*(MT2-p1p2)
                    );
        
//         printf("bPref[%d] =  % .16e %+.16ei\n",i,real(bPref[i]),imag(bPref[i];
        
        // Rest part only for the tau triangle.
        rePref[i] -= e4*MT/(2*MT2-p1p2)/p1p2/4.0/twoPi2*(
                     - A0MT[i] / MT2
                     + B0MT[i] * 2.0
                     - exp[i] // term from expanding B0sMT prefactor in eps
                     - B0sMT[i]);
        
//         printf("nPref[%d] =  % .16e %+.16ei\n",i,real(nPref[i]),imag(nPref[i]));
    }
    
    return;
}


// calculates the prefactors of the box and crossbox diagrams.
void getBoxPref(complex bPref[4], complex abPref[4],complex rePref[4],
                double p1p2, double p1p3, double p1p4, double pref)
{
    int i,poles;
    double p2p3 = p1p4, p2p4 = p1p3;
    
    std::vector<complex> B0box(4), B0cbox(4), C0box(4), C0cbox(4), D0box(4), D0cbox(4);
    
    // First calculate the phase space dependent integrals.
    //calculate the box' B0(p2-p4,0,MT)
    B0(B0box,MT2-2.0*p2p4,0.0,MT2);
    B0box[3] = B0box[1];
    B0box[1] = 0.0;
//     
    //calculate the box' C0(p2-p4,-p1,0,MT,0) = C0(p2,p2-p4,0,0,MT) integral
    C0(C0box,0.0,MT2,MT2-2.0*p2p4,0.0,0.0,MT2);
    C0box[3] = 0.0;
    
    //calculate the box' D0(p2,p2-p4,-p1,0,0,MT,0)=D0(p2,-p4,-p3,p1,p2-p4,p3+p4,{mi})
    D0(D0box,0.0,MT2,MT2,0.0,MT2-2.0*p2p4,2.0*p1p2,0.0,0.0,MT2,0.0);
    D0box[3] = 0.0;
    
//     printf("O(eps^0)  =  % .16e %+.16ei\n",real(D0box[0]),imag(D0box[0]));
//     printf("O(eps^-1) =  % .16e %+.16ei\n",real(D0box[1]),imag(D0box[1]));
//     printf("O(eps^-2) =  % .16e %+.16ei\n",real(D0box[2]),imag(D0box[2]));
    
    //calculate the cbox' B0(p2-p3,0,MT) integral
    B0(B0cbox,MT2-2.0*p2p3,0.0,MT2);
    B0cbox[3] = B0cbox[1];
    B0cbox[1] = 0.0;
//     
    //calculate the cbox' C0(p2-p3,-p1,0,MT,0) = C0(p2,p2-p3,0,0,MT).
    C0(C0cbox,0.0,MT2,MT2-2.0*p2p3,0.0,0.0,MT2);
    C0cbox[3] = 0.0;
    
    //calculate the cbox' D0(p2,p2-p3,-p1,0,0,MT,0)=D0(p2,-p3,-p4,p1,p2-p3,p3+p4,{mi})
    D0(D0cbox,0.0,MT2,MT2,0.0,MT2-2.0*p2p3,2.0*p1p2,0.0,0.0,MT2,0.0);
    D0cbox[3] = 0.0;
    
//     printf("O(eps^0)  =  % .16e %+.16ei\n",real(D0cbox[0]),imag(D0cbox[0]));
//     printf("O(eps^-1) =  % .16e %+.16ei\n",real(D0cbox[1]),imag(D0cbox[1]));
//     printf("O(eps^-2) =  % .16e %+.16ei\n",real(D0cbox[2]),imag(D0cbox[2]));
    
    // decide whether poles are calculated
    if(cEps) poles = 4;
    else poles = 1;
        
    for(i=0;i<poles;i++) {
         // on this level C0(0,0,s,0,0,0) cancels
        bPref[i] -= e2/4.0/twoPi2*(
            // Box part.
                    (
                    - B0MT[i] * 2.0*MT2/p1p3
                    + B0s[i] * 2.0
                    + B0box[i] * 2.0*(MT2-p1p3)/p1p3 
                    )*p1p4/(p1p4 - p1p3)
                    
                    + C0box[i] * 4.0*p1p3
                    + C0s[i] * 2.0*p1p2
                    + C0MTs[i] * 2.0*p1p2*(MT2 - p1p3 + p1p4)/(p1p4-p1p3)
                    
                    + D0box[i] * 4.0*p1p2*p1p3
            // Crossbox part.
                    +(
                    - B0MT[i] * 2.0*MT2/p1p4
                    + B0s[i] * 2.0
                    + B0cbox[i] * 2.0*(MT2-p1p4)/p1p4
                    )*p1p3/(p1p4 - p1p3)
                    
                    - C0cbox[i] * 4.0*p1p4
                    - C0s[i] * 2.0*p1p2
                    + C0MTs[i] * 2.0*p1p2*(MT2 + p1p3 - p1p4)/(p1p4-p1p3)
                    
                    - D0cbox[i] * 4.0*p1p2*p1p4
                    );
        
//         printf("bPref[%d] =  % .16e %+.16ei\n",i,real(bPref[i]),imag(bPref[i]));
        abPref[i] -= e2/4.0/twoPi2*(
            // Box part.
                     (
                     - B0MT[i] * 2.0*MT2 
                     + B0s[i] * 2.0*p1p3
                     + B0box[i] * 2.0*(MT2-p1p3)
                     )/(p1p4 - p1p3)
                     
                     + pref*p1p2*(
                     - C0box[i] * 8.0*p1p3*p1p3*p1p3
                     + C0s[i] * 4.0*p1p2*p1p3*p1p3
                     + C0MTs[i] * 2.0*p1p2*(MT2*MT2*p1p2 + 2.0*(p1p4-p1p3)*p1p3*p1p3
                     + MT2*(6.0*p1p3*p1p3-4.0*p1p2*p1p3))/(p1p4-p1p3)
                     
                     + D0box[i] * 8.0*p1p2*p1p3*p1p3*p1p3
                     )
            // Crossbox part.
                     +(
                     - B0MT[i] * 2.0*MT2 
                     + B0s[i] * 2.0*p1p4
                     + B0cbox[i] * 2.0*(MT2-p1p4)
                     )/(p1p4 - p1p3)
                     
                     + pref*p1p2*(
                     + C0cbox[i] * 8.0*p1p4*p1p4*p1p4
                     - C0s[i] * 4.0*p1p2*p1p4*p1p4
                     + C0MTs[i] * 2.0*p1p2*(MT2*MT2*p1p2 - 2.0*(p1p4-p1p3)*p1p4*p1p4
                     + 2.0*MT2*(p1p2*p1p2-4.0*p1p2*p1p3+3.0*p1p3*p1p3))/(p1p4-p1p3)
                     
                     - D0cbox[i] * 8.0*p1p2*p1p4*p1p4*p1p4
                     )
                    );

        rePref[i] -= e4/4.0/twoPi2*MT/(p1p4 - p1p3)*(
            // Box part.
                     + B0MT[i] * (2.0*MT2 + p1p3 - p1p4)/(2.0*MT2-p1p2)/p1p3
                     - B0s[i] * 2.0/(2.0*MT2-p1p2) 
                     - B0box[i] * 1.0/p1p3 
                     
//                      + C0box[i] * 0.0
//                      + C0s[i] * 0.0
                     - C0MTs[i] * 4.0*MT2/(2.0*MT2-p1p2)
                     
//                      + D0box[i] * 0.0
            // Crossbox part.
                     + B0MT[i] * (2.0*MT2 - p1p3 + p1p4)/(2.0*MT2-p1p2)/p1p4
                     - B0s[i] * 2.0/(2.0*MT2-p1p2) 
                     - B0cbox[i] * 1.0/p1p4
                     
//                      + C0cbox[i] * 0.0
//                      + C0s[i] * 0.0
                     - C0MTs[i] * 4.0*MT2/(2.0*MT2-p1p2)
                     
//                      + D0cbox[i] * 0.0
                    );
        
//         printf("nPref[%d] =  % .16e %+.16ei\n",i,real(nPref[i]),imag(nPref[i]));
    }
    
    return;
}

void virtualMatrix(complex res[4], double pIn[2][4], double pOut[2][4])
{
    int i,j,k,poles;
    
    double p1p2 = metric(pIn[0],pIn[1]),
           p1p3 = metric(pIn[0],pOut[0]),
           p1p4 = p1p2 - p1p3, p2p3 = p1p4, p2p4 = p1p3,
           pref;
    
    
    complex meBorn[4][4],meAntiBorn[4][4],meRest[4][4],summand[4][4][4],sum[4],
            j1[4][4],j2[4][4];
    
    getCurr(pIn,pOut,j1,j2);
    
    currBornME(meBorn,p1p2,p1p3,pref,j1,j2);
    
    currAntiBornME(meAntiBorn,p1p2,p1p3,pref,j1,j2);
    
    currRestME(meRest,j1,j2);
    
//     printf("pref =  % .16e\n",pref);
    
    complex bPref[4] = {0.0,0.0,0.0,0.0},// prefactor of the part factorizing in the born matrix
            abPref[4] = {0.0,0.0,0.0,0.0},// prefactor of the anti born part
            rePref[4] = {0.0,0.0,0.0,0.0};// prefactor of the rest part
            
    
    getBubPref(bPref,p1p2);
    
    getTriPref(bPref,rePref,p1p2);
    
    getBoxPref(bPref,abPref,rePref,p1p2,p1p3,p1p4,pref);
    
    // decide whether poles are calculated
    if(cEps) poles = 4;
    else poles = 1;
           
    for(i=0;i<poles;i++) {
        for(j=0;j<4;j++) {
            summand[i][0][j] = 0.0;
        
            // j1(1,1) part only nonzero for this configuration
            summand[i][1][j] = bPref[i]*meBorn[1][j]
                             + abPref[i]*meAntiBorn[1][j]
                             + rePref[i]*meRest[1][j];
        
            // j1(1,-1) part only  nonzero for this configuration
            summand[i][2][j] = bPref[i]*meBorn[2][j]
                             + abPref[i]*meAntiBorn[2][j]
                             + rePref[i]*meRest[2][j];
        
            summand[i][3][j] = 0.0;
            
            for(k=0;k<4;k++) sum[i] += summand[i][k][j]*conj(meBorn[k][j]);
        }
        res[i] = sum[i]/4.0; // average over polarizations;
    }
    
    
    return;
}
