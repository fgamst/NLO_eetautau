//tree level cross section
void bornIntegrand(double x[2], double *weight, double f[1]) 
{
    double pIn[2][4] = { { Etot/2, 0.0, 0.0, Etot/2 }, { Etot/2, 0.0, 0.0, -Etot/2 } },
           pOut[2][4], jac, vars[2], M2;
    
    f[0] = 0.0;
    
    // generate tau pair final state
    PhSp_TaTa(Etot,x,vars,pOut,&jac);
    
//     M2 = ME_Born_ee_tata(pIn,pOut); // my pure QED version
//     M2 = ME_Born_ee_AZ_tata(pIn,pOut); // my version with photon and Z
    M2 = ME_Born_ee_AZ_tata_check(pIn,pOut); // helas version
    
    binM2(M2,pIn,pOut,*weight,jac);
    
    f[0] += (M2 * jac);
    
    
	return;
}

//frame for calling intDipoles in vegas
void born(unsigned long nCall, int nIt)
{
    int i,j;
    double estim[1];   /* estimators for integrals                     */
    double std_dev[1]; /* standard deviations                          */
    double chi2a[1];   /* chi^2/n                                      */
    double reg[6];   /* integration domain                           */
    int init=0;
    wmax = 0.0;
    countN = 0;
    countIt = 0;
    
    
    // Initializing the integrated Dipoles before integration
    initIntDipoles();
    
    //  initializing the inegration range: 0...1
    for (i=0; i<2; i++) {
        reg[i] = 0.0;
        reg[i+2] = 1.0;
    }
    
    printf("\nPerforming tree level cross section integral\n");
    
    vegas(reg, 2, bornIntegrand, init,nCall,nIt,
          NPRN_INPUT | NPRN_RESULT,
          1, 0, 1, estim, std_dev,chi2a);
    
    printf("total LO cross sec. Estimator   : % .16e\n", estim[0]);
    printf("Error                           : % .16e\n", std_dev[0]);
    printf("chi2a                           : % .16e \n", chi2a[0]);
    
    // divide by the number of iterations.
//     for(i=0; i<3; i++) {
//         for(j=0; j<NBINS+1; j++)
//         fbin[0][j]
//     }
    
    // dump the result at the end of the bins.
    fbin[0][NBINS+1] = estim[0];
    fbin[1][NBINS+1] = std_dev[0];
    fbin[2][NBINS+1] = chi2a[0];
    
    return;
}


// integrate the 3 particle final state
void realCorIntegrand(double x[5], double *weight, double f[1]) 
{
	int i;
    double pIn[2][4] = { { Etot/2, 0.0, 0.0, Etot/2 }, { Etot/2, 0.0, 0.0, -Etot/2 } };
    double vars[5], M3, Dips, jacobian, 
           pOut[3][4];

           
//  generate the phase space 
    PhSp_TaTaPh(Etot,x,vars,pOut,&jacobian);
    
    //version to just check PS3 volume
//     M3 = 1.0;//ME_ee_tataph(pIn,pOut);
//     binM3(M3, pIn, pOut, *weight, jacobian);
//     f[0] = jacobian*M3;
//     return;
    
    
//     countN++;
//     if (countN>NSMP) {
//         countN = 1;
//         countIt++;
//         wmax = 0.0;
//     }
    
    f[0] = 0.0;
    M3 = 0.0;
    Dips = 0.0;
//         double pT = sqrt(pOut[2][1]*pOut[2][1]+pOut[2][2]*pOut[2][2]);
//     if (pT<acut) return;

//  check the cuts
    if (pOut[2][0]*pOut[2][0]/shat < scut // soft cut condition
        || metric(pIn[0],pOut[2])/shat < ccut // collinear cut conditions
        || metric(pIn[1],pOut[2])/shat < ccut) {
        return;
    }
     else {
        M3 = ME_ee_tataph(pIn,pOut);
        
        binM3(M3, pIn, pOut, *weight, jacobian);
        
        Dips = getDipSumOld(pIn,pOut,jacobian,*weight);
        // For testing.
//         Dips = getDipIF(pIn,pOut,pInTilde,pOutTilde,1,2,1);
//         Dips = getDipFI(pIn,pOut,pInTilde,pOutTilde,1,2,1);
//         Dips = getDipFF(pIn,pOut,pInTilde,pOutTilde,0,2,1);
    }
    
    f[0] = jacobian*(M3-Dips);
//     f[0] = 1/(metric(pIn[0],pOut[2]));
    if (f[0] != f[0]
        ||M3 != M3
        ||Dips != Dips) {
        printf("\nError encountered!\n");
        printf("x[1] = % .10e\n", x[0]);
        printf("x[2] = % .10e\n", x[1]);
        printf("x[3] = % .10e\n", x[2]);
        printf("x[4] = % .10e\n", x[3]);
        printf("x[5] = % .10e\n", x[4]);
        printf("jac  = % .10e\n", jacobian);
        printf("Dips = % .10e\n", Dips);
        printf("Dips = % .10e\n", M3);
        f[0] = 0.0;
    }
    
//     if (countIt<1) return;
//     if( fabs(f[0] * *weight) > wmax ){
//         wmax = fabs(f[0] * *weight);
//         printf("wmax     : % .10e\n",wmax);
//         printf("weight   : % .10e\n",*weight);
//         printf("jacobian : % .10e\n",jacobian);
//         printf("M3       : % .10e\n",M3);
//         printf("Dips     : % .10e\n",Dips);
//         printf("\n");
//         printf("p1      : % .10e % .10e % .10e % .10e\n",
//                    pIn[0][0],pIn[0][1],pIn[0][2],pIn[0][3]);
//         printf("p2      : % .10e % .10e % .10e % .10e\n",
//                    pIn[1][0],pIn[1][1],pIn[1][2],pIn[1][3]);
//         printf("p3      : % .10e % .10e % .10e % .10e\n",
//                    pOut[0][0],pOut[0][1],pOut[0][2],pOut[0][3]);
//         printf("p4      : % .10e % .10e % .10e % .10e\n",
//                    pOut[1][0],pOut[1][1],pOut[1][2],pOut[1][3]);
//         printf("p5      : % .10e % .10e % .10e % .10e\n",
//                    pOut[2][0],pOut[2][1],pOut[2][2],pOut[2][3]);
//       printf("\n");
//     }
    
	return;
}

//frame for calling realCorrection in vegas
void realCor(unsigned long nCall, int nIt)
{
    int i;
    double estim[1];   /* estimators for integrals                     */
    double std_dev[1]; /* standard deviations                          */
    double chi2a[1];   /* chi^2/n                                      */
    double reg[10];   /* integration domain                           */
    int init=0;
    
    wmax = 0.0;
    countN = 0;
    countIt = 0;
    
//     double muj = MUT, muj2 = muj*muj, muk = MUT, muk2 = muk*muk,
//            yp = 1 - 2*muk*(1-muk)/(1-muj2-muk2),
//            acut_FFp = acut_FF*yp;
//            
//     printf("alpha_FFp = %.5f\n", acut_FFp);
    
    // Initializing the integrated Dipoles before integration
    initIntDipoles();
    
    //  initializing the inegration range: 0...1
    for (i=0; i<5; i++) {
        reg[i] = 0.0;
        reg[i+5] = 1.0;
    }
    
    printf("\nPerforming total real correction integral\n");
    vegas(reg, 5, realCorIntegrand, init, nCall, nIt,
          NPRN_INPUT | NPRN_RESULT, 
          1, 0, 1, estim, std_dev, chi2a);
    
    printf("\n");
    printf("total real Cor. Estimator       : % .16e\n", estim[0]);
    printf("Error                           : % .16e\n", std_dev[0]);
    printf("chi2a                           : % .16e\n", chi2a[0]);
    
    // dump the result at the end of the bins.
    fbin[0][NBINS+1] = estim[0];
    fbin[1][NBINS+1] = std_dev[0];
    fbin[2][NBINS+1] = chi2a[0];
    
    return;
}


//integrated dipoles to be called by vegas.
void intDipIntegrand(double x[3], double *weight, double f[3]) 
{
    int i;
    double pIn[2][4] = { { Etot/2, 0.0, 0.0, Etot/2 }, { Etot/2, 0.0, 0.0, -Etot/2 } },
           pOut[2][4], pInTilde[2][4], pOutTilde[2][4], Dips[3];
//     double xint = x[2];//*(1-4*MT2/shat)+4*MT2/shat; // for the phase space convolution.
    // This will also change the jacobian!!
    
    for(i=0;i<3;i++) f[i] = 0.0;
    
    getIntDipSum(x, pIn, pOut, pInTilde, pOutTilde, *weight,Dips);
    
    for(i=0;i<3;i++) f[i] += Dips[i];
    
	return;
}

//frame for calling intDipoles in vegas
void intDip(unsigned long nCall, int nIt)
{
    int i;
    double estim[3];   /* estimators for integrals                     */
    double std_dev[3]; /* standard deviations                          */
    double chi2a[3];   /* chi^2/n                                      */
    double reg[6];   /* integration domain                           */
    int init=0;
    wmax = 0.0;
    countN = 0;
    countIt = 0;
    
    printf("\nmuR2 = %f\n",muR2);
    
    // Initializing the integrated Dipoles before integration
    initIntDipoles();
    
    printf("\nmuR2 = %f\n",muR2);
    
    //  initializing the inegration range: 0...1
    for (i=0; i<3; i++) {
        reg[i] = 0.0;
        reg[i+3] = 1.0;
    }
    
    printf("\nPerforming total integrated Dipole integral\n");
    
    vegas(reg, 3, intDipIntegrand, init,nCall,nIt,
          NPRN_INPUT | NPRN_RESULT,
          (1+cEps*2), 0, 1, estim, std_dev,chi2a);
    
    printf("total int. Dipole Estimator     : % .16e\n", estim[0]);
    printf("Error                           : % .16e\n", std_dev[0]);
    printf("chi2a                           : % .16e\n\n", chi2a[0]);
    
    if(cEps) {
        printf("O(eps^-1)                       : % .16e\n", estim[1]);
        printf("Error                           : % .16e\n", std_dev[1]);
        printf("chi2a                           : % .16e\n\n", chi2a[1]);
    
        printf("O(eps^-2)                       : % .16e\n", estim[2]);
        printf("Error                           : % .16e\n", std_dev[2]);
        printf("chi2a                           : % .16e\n\n", chi2a[2]);
    }
    
    // dump the result at the end of the bins.
    fbin[0][NBINS+1] = estim[0];
    fbin[1][NBINS+1] = std_dev[0];
    fbin[2][NBINS+1] = chi2a[0];
    
    printf("\nmuR2 = %f\n",muR2);
    
    return;
}

void virtCorIntegrand(double x[2], double *weight, double f[4]) 
{
	int i,j,poles;
    double pIn[2][4] = { { Etot/2, 0.0, 0.0, Etot/2 }, { Etot/2, 0.0, 0.0, -Etot/2 } },
           pOut[2][4],vars[2],jac,M2;
    
    complex res[4];
           
    PhSp_TaTa(Etot,x,vars,pOut,&jac);
    
    virtualMatrix(res,pIn,pOut);
    
    // decide whether poles are calculated
    if(cEps) poles = 4;
    else poles = 1;
    
    M2 = 2*real(res[0]);
    
    binM2(M2, pIn, pOut, *weight, jac);
    
    f[0] = M2*jac;
    
    for(i=1;i<poles;i++) f[i] = 2*real(res[i])*jac;
    
	return;
}

void virtCor(unsigned long nCall, int nIt)
{
    int i,poles;
    double estim[4];   /* estimators for integrals                     */
    double std_dev[4]; /* standard deviations                          */
    double chi2a[4];   /* chi^2/n                                      */
    double reg[4];   /* integration domain                           */
    int init=0;
    wmax = 0.0;
    countN = 0;
    countIt = 0;
    //  initializing the inegration range: 0...1
    for (i=0; i<2; i++) {
        reg[i] = 0.0;
        reg[i+2] = 1.0;
    }
    
    // Initializing the integrated Dipoles before integration
    initVirtCor();
    
    printf("\nPerforming virtual Correction integral\n");
    
    // decide whether poles are calculated
    if(cEps) poles = 4;
    else poles = 1;
    
    vegas(reg, 2, virtCorIntegrand, init,nCall,nIt,
          NPRN_INPUT | NPRN_RESULT,
          poles, 0, 1, estim, std_dev,chi2a);
    
    printf("total virtual Cor. Estimator    : % .16e\n", estim[0]);
    printf("Error                           : % .16e\n", std_dev[0]);
    printf("chi2a                           : % .16e\n\n", chi2a[0]);
//   
    if(cEps) {
    printf("IR O(eps^-1)                    : % .16e\n", estim[1]);
    printf("Error                           : % .16e\n", std_dev[1]);
    printf("chi2a                           : % .16e\n\n", chi2a[1]);
    
    printf("IR O(eps^-2)                    : % .16e\n", estim[2]);
    printf("Error                           : % .16e\n", std_dev[2]);
    printf("chi2a                           : % .16e\n\n", chi2a[2]);
    
    printf("UV O(eps^-1)                    : % .16e\n", estim[3]);
    printf("Error                           : % .16e\n", std_dev[3]);
    printf("chi2a                           : % .16e\n\n", chi2a[3]);
    }
    
    // dump the result at the end of the bins.
    fbin[0][NBINS+1] = estim[0];
    fbin[1][NBINS+1] = std_dev[0];
    fbin[2][NBINS+1] = chi2a[0];
    
    return;
}

void analyticPolesIntegrand(double x[2], double *weight, double f[1]) 
{
	int i,j;
    double pIn[2][4] = { { Etot/2, 0.0, 0.0, Etot/2 }, { Etot/2, 0.0, 0.0, -Etot/2 } },
           pOut[2][4],vars[2],jac,vijk,rho,pref,M2;
    
//     complex res[3];
           
    PhSp_TaTa(Etot,x,vars,pOut,&jac);
    
    vijk = sqrt(1-4*MUT*MUT)/(1-2*MUT*MUT);
    
    rho = sqrt((1-vijk)/(1+vijk));
    
    M2 = ME_Born_ee_tata(pIn,pOut);
    
    pref = (M2 * jac)*(alpha)/M_PI;
    
    f[0] = (0.0
//             + 1.5 + log(mu2/shat) // electron triangle part.
//             + 1.0 + log(rho)/vijk // tau triangle part.
            - 2.0*log(muR2/2/metric(pIn[0],pOut[0])) // log(mu2/2/metric(p1,p3)) box part.
            + log(muR2/MT2) // box part.
//             + 2.0*log(muR2/2/metric(pIn[0],pOut[1])) // log(mu2/2/metric(p1,p4)) crossbox part.
//             - log(muR2/MT2) // crossbox part.
           )*pref;
    
    f[1] = (0.0
//             + 1.0 // electron triangle part.
            - 1.0 // box part.
//             + 1.0 // crossbox part.
           )*pref;
    
	return;
}

void analyticPoles(unsigned long nCall, int nIt)
{
    int i;
    double estim[3];   /* estimators for integrals                     */
    double std_dev[3]; /* standard deviations                          */
    double chi2a[3];   /* chi^2/n                                      */
    double reg[4];   /* integration domain                           */
    int init=0;
    wmax = 0.0;
    countN = 0;
    countIt = 0;
    //  initializing the inegration range: 0...1
    for (i=0; i<2; i++) {
        reg[i] = 0.0;
        reg[i+2] = 1.0;
    }
    
    printf("\nPerforming integration of analytic Poles\n");
    
    vegas(reg, 2, analyticPolesIntegrand, init,nCall,nIt,
          NPRN_INPUT | NPRN_RESULT,
          2, 0, 1, estim, std_dev,chi2a);
    
//     printf("total virtual Cor. Estimator    : % .16e\n", estim[0]);
//     printf("Error                           : % .16e\n", std_dev[0]);
//     printf("chi2a                           : % .16e\n\n", chi2a[0]);
//     
    printf("O(eps^-1)                       : % .16e\n", estim[0]);
    printf("Error                           : % .16e\n", std_dev[0]);
    printf("chi2a                           : % .16e\n\n", chi2a[0]);
    
    printf("O(eps^-2)                       : % .16e\n", estim[1]);
    printf("Error                           : % .16e\n", std_dev[1]);
    printf("chi2a                           : % .16e\n\n", chi2a[1]);
    
    // dump the result at the end of the bins.
//     fbin[0][NBINS+1] = estim[0];
//     fbin[1][NBINS+1] = std_dev[0];
//     fbin[2][NBINS+1] = chi2a[0];
    
    return;
}
