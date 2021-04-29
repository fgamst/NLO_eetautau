/* integrated dipoles for Catani-Seymour factorization formula. */

// initialize some stuff so as to not calculate it all the time
void initIntDipoles()
{
    
    double /*mui = 0.0,*/ muj = MUT, muj2 = muj*muj, muk = MUT, muk2 = muk*muk,
            // mass distribution for Tauon pair production
           vijkT = sqrt(kallen(1,muj2,muk2))/(1-muj2-muk2), 
            // see the notebook int_dipoles.nb
           rho = sqrt((1-vijkT)/(1+vijkT)),//for all of these see massive Catani-Seymour 5.30.
           rhoj = sqrt((1-vijkT+2*muj2/(1-muj2-muk2))/(1+vijkT+2*muj2/(1-muj2-muk2))),
           rhok = sqrt((1-vijkT+2*muk2/(1-muj2-muk2))/(1+vijkT+2*muj2/(1-muj2-muk2)));
           
    // eikonal part of the integral.  Catani-Seymour 5.34
    // finite part.
    Ieik[0] = (-log(rho)*log(1-(muj+muk)*(muj+muk))
           - log(rhoj)*log(rhoj)/2 - log(rhok)*log(rhok)/2 + M_PI*M_PI/6
           + 2*Li2(-rho) - 2*Li2(1-rho) 
           - Li2(1-rhoj*rhoj)/2 - Li2(1-rhok*rhok)/2
           + logMuShat*log(rho)/2.0)/vijkT; // from expansion of the prefactor.
    // 1/eps part.
    Ieik[1] = log(rho)/vijkT/2.0;
    
    // collinear part of the integral. Catani-Seymour 5.35
    // finite part.
    IgQcol[0] = 3 - log(muj) -2*log((1-muk)*(1-muk) - muk2) + log(1-muk) - muk/(1-muk) 
             - 2*muj2/(1-muj2-muk2)*log(muj/(1-muk)) - 2*muk*(1-2*muk)/(1-muj2-muk2)
             + logMuShat; // from expansion of the prefactor.
    // 1/eps part.
    IgQcol[1] = 1.0;
             
    // calculate the upper limit of the yijk integral
    yp = 1 - 2*muk*(1-muk)/(1-muj2-muk2);
    acut_FFp = acut_FF*yp;
    
    // alpha contribuation to integrated dipoles Czakon A.20
    Del_IgQcol = 3*(1+acut_FFp)/2 + 1/(1-muk) - 2*(2-2*muj2-muk)/(1-muj2-muk2)
             + (1-acut_FFp)*muj2/2/(muj2+acut_FFp*(1-muj2-muk2))
             -2*log(acut_FFp*(1-muj2-muk2)/((1-muk)*(1-muk)-muj2)) 
             + (1+muj2-muk2)/2/(1-muj2-muk2)
             *log((muj2+acut_FFp*(1-muj2-muk2))/(1-muk)/(1-muk));
             
    double x = 0.0, // Czakon A.2
           a = 2*muk/(1-muj2-muk2), // Czakon A.10
           b = 2*(1-muk)/(1-muj2-muk2), // Czakon A.11
           c = 2*muk*(1-muk)/(1-muj2-muk2), // Czakon A.12
           d = (1-muj2-muk2)/2, // Czakon A.13
           vjk = vijkT, // Czakon A.15 they are the same
           xp = ((1-muk)*(1-muk)-muj2)/(1-muj2-muk2) + vjk, // Czakon A.14
           xm = ((1-muk)*(1-muk)-muj2)/(1-muj2-muk2) - vjk; // Czakon A.14
           
           
    if(acut_FFp < 0.704) // above this threshold x becomes complex.
        x = yp - acut_FFp + sqrt((yp-acut_FFp)*(1/yp - acut_FFp
               + 4*muj2*muk2/(muj2-(1-muk)*(1-muk))/(1-muj2-muk2))); // Czakon A.8
//     printf("\nx          = % .16e\n",x);
    // alpha contribuation to integrated dipoles Czakon A.20
    Del_Ieik = (-Li2((a+x)/(a+xp))+Li2((a)/(a+xp))// check
                  +Li2((xp-x)/(xp-b))-Li2((xp)/(xp-b))// check
                  
                  +Li2((c+x)/(c+xp))-Li2((c)/(c+xp))// check
                  +Li2((xm-x)/(xm+a))-Li2((xm)/(xm+a))// check
                  
                  -Li2((b-x)/(b-xm))+Li2((b)/(b-xm))// check
                  -Li2((xm-x)/(xm+c))+Li2((xm)/(xm+c))// check
                  
                  +Li2((b-x)/(b+a))-Li2((b)/(b+a))// check
                  -Li2((c+x)/(c-a))+Li2((c)/(c-a))// check
    
                  +log(c+x)*log((a-c)*(xp-x)/(a+x)/(c+xp))
                  -log(c)*log((a-c)*xp/a/(c+xp))// check
                  
                  +log(b-x)*log((a+x)*(xm-b)/(a+b)/(xm-x))
                  -log(b)*log(a*(xm-b)/(a+b)/xm)// check
                  
                  -log((a+x)*(b-xp))*log(xp-x)+log(a*(b-xp))*log(xp)// check
                  
                  +log(d)*log((a+x)*xp*xm/a/(xp-x)/(xm-x))
                  +log((xm-x)/xm)*log((c+xm)/(a+xm))
                  +log((a+x)/a)*log(a*(a+x)*(a+xp)*(a+xp))/2)/vjk;// check
    
//     printf("\n");
//     printf("vjk        = % .10e\n",vjk);
//     printf("-Li2((a+x)/(a+xp)) = % .10e\n",-Li2((a+x)/(a+xp)));
//     printf("+Li2((a)/(a+xp))   = % .10e\n",+Li2((a)/(a+xp)));
//     printf("+Li2((xp-x)/(xp-b))= % .10e\n",+Li2((xp-x)/(xp-b)));
//     printf("-Li2((xp)/(xp-b))  = % .10e\n",-Li2((xp)/(xp-b)));
//     printf("muj        = % .16e\n",muj);
//     printf("muk        = % .16e\n",muk);
//     printf("x          = % .16e\n",x);
//     printf("x+         = % .16e\n",xp);
//     printf("x-         = % .16e\n",xm);
//     printf("a          = % .16e\n",a);
//     printf("b          = % .16e\n",b);
//     printf("c          = % .16e\n",c);
//     printf("d          = % .16e\n",d);
//     printf("Del_Ieik   = % .16e\n",Del_Ieik);
//     printf("Del_IgQcol = % .16e\n",Del_IgQcol);
}


/* integrated Dipoles */

/* initial emitter and initial observer. See massless Catani-Seymour Sec 5.5 Eq. 5.155
   for the respective phasespace see Eq. 5.149 and following.
       
   result part is split between 0: DiracDelta part, 1: PlusDistro part, 2: Regular part.
   */
void getIntDipII(double x, double result[3][4]) 
{   
    // Sort the results into Dirac, PlusDist and regular part
    // For the results see mathematica notebook "int_dipoles.nb"
    double plusTerm;
    
    
    result[0][0] = -M_PI*M_PI/(6.0); // DiracDelta part   
//     if(colQED) result[0][0] -= 1.5*logMeShat; // Dittmaier 2000 collinear part
    
    // poles for cross checking
    if(cEps) {
        // 1/eps^2
        result[2][0] = 1.0; //part from massless CS 5.32
        
        // 1/eps with expansion part from 1/eps^2
        result[1][0] = 1.5 + logMuShat; // part from massless CS 5.32
    }
    
    // resulting expansion part from poles
    result[0][0] += 1.5*logMuShat + logMuShat*logMuShat/2.0;
//     result[0][0] += logMuShat*logMuShat/2.0;// part where the collinear singularity is replaced
    
    // avoid the singularity at x=1
    if(x < 1.0-ecut) {
        plusTerm = 4.0*log(1.0-x)/(1.0-x);
        
        result[0][3] = plusTerm; // x=1-part of PlusDistro
        if(colQED) result[0][3] -= logMeShat*2/(1.0-x); // Dittmaier 2000 collinear part
    
        if(shat*x < 4*(MT2)) return; // dont calculate parts with Mx below the threshold
        
        result[0][2] = plusTerm;  // x-part of PlusDistro
        if(colQED) result[0][2] -= logMeShat*(1+x*x)/(1.0-x); // Dittmaier 2000 collinear part
    
        // regular part
        result[0][1] = (1.0-x)+(1.0+x)*(log(x)-2.0*log(1.0-x))-2.0*log(x)/(1.0-x); 
        
        // 1/eps part
        //result[1][1] = -P^qq(x); // nonperturbative collinear part from massles CS 5.83.
        
        // resulting expansion part
        //result[0][1] = -P^qq(x)*shat;

        // Add the alpha parameter part. A
        if ((1.0-x) >= acut_II) 
            result[0][1] += (1.0+x*x)/(1.0-x)*log((acut_II)/(1.0-x));
    }
    
    
    return;
}

/* initial emitter and final observer. See massive Catani-Seymour Sec 5.3 Eq. 5.88
 * and for alpha dependencies alpha_paper appendix A eq A9
   for the respective phasespace see Eq. 5.78 and following.
   
   indices:
   a emitter
   j observer

   */
void getIntDipIF(double x, double sja, double sTja, double mu2, double muT2,
//                  double pIn[2][4],double pOut[2][4],
//                  double pInTilde[2][4],double pOutTilde[2][4],
                 double result[3][4])  
{   
    double //mu2 = MT2/2/metric(pIn[a],pOut[j]),
//            muT2 = MT2/2/metric(pIn[a],pOutTilde[j]),
           muTT2 = muT2/x,
           zp = (1-x)/(1-x+muT2),// is it muTT or muT??
           logmu = log(1+mu2),
           logs = log(sja), // actually this is log(2sja)
           logsT = log(sTja),
           logmus = logMu - logs, // log(mu^2/2sja)
           logMesja = logMe - logs, // log(me^2/2sja)
           logMesTja = logMe - logsT, // log(me^2/2sja)
           //logmusT = logMu - logsT, // log(mu^2/2sT)
           plusPart = 0.0; 
           
           
           
    // Sort the results into Dirac, regular and PlusDist part
    // For these result see mathematica notebook "int_dipoles.nb" and A9 in alpha paper
    result[0][0] = //M_PI*M_PI/(6.0) not sure this is correct
                   +2*Li2(-mu2)+2*log(mu2)*logmu-logmu*logmu/2; // DiracDelta part
//     if(colQED) result[0][0] -= 1.5*logMesja; // Dittmaier 2000 collinear part
    
    // poles for cross checking
    if(cEps) {
        // 1/eps^2
        result[2][0] = 1.0; // massless part from massive CS 5.90
        
        // 1/eps
        result[1][0] = logmu // massive part from massive CS 5.88
                + 1.5 // massless part from massive CS 5.90
                + logmus; // expansion of 1/eps^2 part
                
        //result[1][1] = -Pqq(x); // nonperturbative collinear part.
        //result[0][1] = -Pqq(x) * logmus; // nonperturbative collinear part from expansion.
    }
    //resulting expansion parts from poles:
    result[0][0] += (logmu + 1.5 + logmus/2.0)*logmus;
    
    if(x < 1.0-ecut) { // avoid the singularity at x=1
        plusPart = 4*log(1-x)/(1-x);
        result[0][3] = plusPart - 2*logmu/(1-x); // x=1-part of PlusDistro
        if(colQED) result[0][3] -= logMesja*2/(1.0-x); // Dittmaier 2000 collinear part
    }
    
    if(muT2 < 0.0) return; // dont calculate parts with Mx below the threshold
    
    if(x < 1.0-ecut) {
        result[0][2] = plusPart - 2*log(1+muTT2)/(1-x); // x-part of PlusDistro
        if(colQED) result[0][2] -= logMesja*(1+x*x)/(1.0-x); // Dittmaier 2000 collinear part
        
    
        result[0][1] = -2/(1-x)*log((2-x+muT2)/(1+muTT2)) - (1+x)*log((1-x)*(1-x)/(1-x+muT2))
                + (1 - x); // regular part

        // Add the alpha parameter part. See alpha paper A9
        if (zp >= acut_IF) 
            result[0][1] += -(2/(1.0-x)*log((zp*(1.0-x+acut_IF))/(acut_IF*(1.0-x+zp)))
                                 -(1.0+x)*log(zp/acut_IF));
    }
    else {
//         result[0][1] = 0.0; // regular part
        
        if (zp >= acut_IF) result[0][1] += 2/zp + 2*log(zp/acut_IF) - 2/acut_IF;
    }
    
    return;
}


/* final emitter and initial observer. See massless Catani-Seymour Sec 5.2.3 Eq. 5.58-60
 * as well as alpha paper A13/14
   for the respective phasespace see Eq. 5.53 and following.
   
   indices:
   a emitter
   j observer
       
   
   */
void getIntDipFI(double x, double sja, double sTja, double mu2, double muT2,
//                  double pIn[2][4],double pOut[2][4],
//                  double pInTilde[2][4],double pOutTilde[2][4],
                 double result[3][4]) 
{   
    double //muT2 = MT2/2/metric(pIn[a],pOutTilde[j]),
           muTT2 = muT2/x,//MT2/2/metric(pInTilde[a],pOutTilde[j]),
//            mu2 = MT2/2/metric(pIn[a],pOut[j]);
           lnmu = log(mu2),
           ln1mu = log(1+mu2),
           logmu = ln1mu - lnmu - 1.0,//log((1+mu2)/mu2)-1.0
           logs = log(sja),
           //logsT = log(sTja),
           logmus = logMu - logs; // log(mu^2/sja)
           //logmusT = logMu - logsT; // log(mu^2/sTja)
    // Sort the results into Dirac, PlusDist and regular part
    // For these result see mathematica notebook "int_dipoles.nb" as well as alpha paper A13/14.
    
//     printf("logmus  =  % .16e\n",logmus);
//     printf("logmusT =  % .16e\n",logmusT);
    
    // DiracDelta part
    result[0][0] = -2*Li2(-mu2) - M_PI*M_PI/3 + 2 
                +lnmu*lnmu/2 + ln1mu*ln1mu/2 -2*lnmu*ln1mu + lnmu;
    
    result[0][0] += 2*log(acut_FI)*logmu; // alpha contribution
    
    
           
    // poles for cross checking
    if(cEps) {
        // 1/eps^2 is 0 since the final state particles mass does not vanish
        // 1/eps
        result[1][0] = -logmu;//1.0 + log(mu2/(1+mu2))
    }
    
    result[0][0] += -logmus*logmu; // resulting expansion part from poles
    
    // Add the alpha parameter part. See alpha paper A13/14?
    if (x >= 1 - acut_FI) {
        if(x < 1.0-ecut) result[0][3] = 2/(1-x)*logmu; // x=1-part of PlusDistro
        
        if (muT2 < 0.0) return;
        
        // x-part of PlusDistro
        if(x < 1.0-ecut) {
            result[0][2] = 2/(1-x)*(log((1+muTT2)/muTT2)-1);
            
            // regular part this guy becomes nan if x goes too close to 1
            result[0][1] = (1-x)/2/(1-x+muT2)/(1-x+muT2) 
            + 2/(1-x) * log((2-x+muT2)*muTT2/(1+muTT2)/(1-x+muT2));
//             printf("term 1    =  % .16e\n",(1-x)/2/(1-x+muT2)/(1-x+muT2));
//             printf("term 2    =  % .16e\n",
//                         2/(1-x) * log((2-x+muT2)*muTT2/(1+muTT2)/(1-x+muT2)));
        }
        // limit from mathematica
        else result[0][1] = -2/(mu2+mu2*mu2);
        
    }
    
    return;
}


/* final emitter and final observer. 
 * See massive Catani-Seymour Sec 5.1 Eqs. 5.23, 5.26, 5.34, 5.35
   for the respective phasespace see Eq. 5.10.

 */
void getIntDipFF(double result[3][4]) 
// result part is split between 
// 0: DiracDelta part, 1: Regular part, 2: PlusDistro part .
{   
    // Catani-Seymour 5.23 and alpha FF paper A.17.
    // delta part.
    
    result[0][0] = /*CF**/2*Ieik[0]+IgQcol[0];
    
    result[0][0] += /*CF**/2*Del_Ieik+Del_IgQcol;
    
    // poles for cross checking
    if(cEps) {
        // 1/eps^2 is 0!
        // 1/eps
        result[1][0] = /*CF**/2*Ieik[1]+IgQcol[1];
    }
    
    //resulting expansion parts from poles:
//     result[0][0] += (2*Ieik[1]+IgQcol[1])*logMuShat;
    return;
}

void getIntDipSum(double x[3], double pIn[2][4], double pOut[2][4],
                    double pInTilde[2][4], double pOutTilde[2][4],
                    double wgt, double res[3])
{
    int i,j,k,l, thresh = 0;
    double vars[2], M1, Mx = 0.0, jac1, jacx = 0.0;
    
    
    PhSp_TaTa(Etot, x, vars, pOut, &jac1); // get the phase space for the untransformed kinematics.
    
    M1 = ME_Born_ee_tata(pIn, pOut); // get the unscaled squared matrix element.
    
//     printf("p3    : % .10e % .10e % .10e % .10e\n",
//                 pOut[0][0],pOut[0][1],pOut[0][2],pOut[0][3]);
//     printf("p4    : % .10e % .10e % .10e % .10e\n",
//                 pOut[1][0],pOut[1][1],pOut[1][2],pOut[1][3]);
    
    double costh = pOut[bin][3]/abs3(pOut[bin]), costhp[2] = {costh, costh}, 
           pT = sqrt(pOut[bin][1]*pOut[bin][1]+pOut[bin][2]*pOut[bin][2]), pTp = pT,
           y = 0.5*log((pOut[bin][0]+pOut[bin][3])/(pOut[bin][0]-pOut[bin][3])), yp = y;
            
    double scalar[2][2], scalarT[2][2], // scalar products of regular and rescaled momenta
           mu2[2][2], muT2[2][2]; // reduced masses of regular and rescaled particles
    
    //     for e-e+->ta+ta- particle scattering via momentum conservation we have
    scalar[0][0] = metric(pIn[0],pOut[0]); // p1.p3
    scalar[0][1] = shat/2 - scalar[0][0]; // p1.p4 = p1.p2 - p1.p3
    scalar[1][0] = scalar[0][1]; // p2.p3 = p1.p4
    scalar[1][1] = scalar[0][0]; // p2.p4 = p1.p3;
           
    //calculate the products and reduced masses.
    for(i=0;i<2;i++) {
        for(j=0;j<2;j++) {
//             scalar[i][j] = metric(pIn[i],pOut[j]);
            mu2[i][j] = MT2/2/scalar[i][j];
            muT2[i][j] = -1.0;
        }
    }
        
//     check if the threshold condition is met.
    double xint = x[2], sprime = xint*shat; // rescaled phase space energy.
    
    if (sprime >= 4*(MT2)) {//return; 
        thresh = 1;
        // boost to the rescaled center of mass
        double E12 = sqrt(sprime), // calculate the new CMS energy
               gamma = ((1+xint)*Etot)/2/E12,
               beta = ((1-xint)*Etot)/2/E12/gamma, 
               //calculate the boost parameters
               p3Til[2];
               
        
//         printf("gamma  = % .6f\n",gamma);
//         printf("beta   = % .6f\n",beta);

        pInTilde[0][0] = E12/2;
        pInTilde[0][1] = 0.0;
        pInTilde[0][2] = 0.0;
        pInTilde[0][3] = pInTilde[0][0];
        
        for (i=0;i<3;i++) pInTilde[1][i] = pInTilde[0][i];
        pInTilde[1][3] = -pInTilde[0][3];
        
        PhSp_TaTa(E12, x, vars, pOutTilde, &jacx); // generate the rescaled phase space.
//         printf("pT3   : % .10e % .10e % .10e % .10e\n",
//                 pOutTilde[0][0],pOutTilde[0][1],pOutTilde[0][2],pOutTilde[0][3]);
//         printf("pT4   : % .10e % .10e % .10e % .10e\n",
//                 pOutTilde[1][0],pOutTilde[1][1],pOutTilde[1][2],pOutTilde[1][3]);
        p3Til[0] = gamma*(pOutTilde[bin][3] - beta*pOutTilde[bin][0]); // boost to rescaled system.
        p3Til[1] = gamma*(pOutTilde[bin][3] + beta*pOutTilde[bin][0]); // boost to rescaled system.
//         printf("p3Til1 = % .6f\n",p3Til[0]);
//         printf("p3Til2 = % .6f\n",p3Til[1]);
        for(i=0;i<2;i++) 
            costhp[i] = p3Til[i]/sqrt(pOutTilde[bin][1]*pOutTilde[bin][1]
                                     +pOutTilde[bin][2]*pOutTilde[bin][2]
                                     +p3Til[i]*p3Til[i]);
            
        pTp = sqrt(pOutTilde[bin][1]*pOutTilde[bin][1]+pOutTilde[bin][2]*pOutTilde[bin][2]);
        yp = 0.5*log((pOutTilde[bin][0]+pOutTilde[bin][3])/(pOutTilde[bin][0]-pOutTilde[bin][3]));
    
        Mx = ME_Born_ee_tata(pInTilde, pOutTilde); // calculate the rescaled matrix element.
//         muT2 = MT2/2/metric(pIn[i],pOutTilde[j]);
        
        for(i=0;i<2;i++) {
//             for(j=0;j<2;j++) {
//                 scalarT[i][j] = metric(pIn[i],pOutTilde[j]);
//                 muT2[i][j] = MT2/2/scalarT[i][j];
//             }
            scalarT[i][0] = metric(pIn[i],pOutTilde[0]);
            scalarT[1-i][1] = scalarT[i][0];
            muT2[i][0] = MT2/2/scalarT[i][0];
            muT2[1-i][1] = muT2[i][0];
        }
    }
//     printf("mu2    : % .10e % .10e % .10e % .10e\n",
//            mu2[0][0],mu2[0][1],mu2[1][0],mu2[1][1]);
//     printf("muT2   : % .10e % .10e % .10e % .10e\n",
//            muT2[0][0],muT2[0][1],muT2[1][0],muT2[1][1]);
    
    double D[4][4][3][4];
    
    for(i=0;i<4;i++) {
        for(j=0;j<4;j++) {
            for(k=0;k<3;k++) {
                for(l=0;l<4;l++) D[i][j][k][l] = 0.0;
            }
        }
    }
    
//     Call the different integrated Dipoles.
    getIntDipII(xint,D[0][1]);
//     printf("iD12   : % .10e % .10e % .10e % .10e\n",
//            D[0][1][0],D[0][1][1],D[0][1][2],D[0][1][3]);
    
    getIntDipIF(xint,2*scalar[0][0],2*scalarT[0][0],mu2[0][0],muT2[0][0],D[0][2]);
//     printf("iD13   : % .10e % .10e % .10e % .10e\n",
//            D[0][2][0],D[0][2][1],D[0][2][2],D[0][2][3]);
    
    getIntDipIF(xint,2*scalar[0][1],2*scalarT[0][1],mu2[0][1],muT2[0][1],D[0][3]);
//     printf("iD14   : % .10e % .10e % .10e % .10e\n",
//            D[0][3][0],D[0][3][1],D[0][3][2],D[0][3][3]);
    
//     getIntDipII(xint,D[1][0]);
//     printf("iD21   : % .10e % .10e % .10e % .10e\n",
//            D[1][0][0],D[1][0][1],D[1][0][2],D[1][0][3]);
    
//     //iD23 = iD14 so delete this.
//     getIntDipIF(xint,mu2[1][0],muT2[1][0],D[1][2]);
//     printf("iD23   : % .10e % .10e % .10e % .10e\n",
//            D[1][2][0],D[1][2][1],D[1][2][2],D[1][2][3]);
    
//     //iD24 = iD13 so delete this.
//     getIntDipIF(xint,mu2[1][1],muT2[1][1],D[1][3]);
//     printf("iD24   : % .10e % .10e % .10e % .10e\n",
//            D[1][3][0],D[1][3][1],D[1][3][2],D[1][3][3]);
    
    getIntDipFI(xint,2*scalar[0][0],2*scalarT[0][0],mu2[0][0],muT2[0][0],D[2][0]);
//     printf("iD31   : % .10e % .10e % .10e % .10e\n",
//            D[2][0][0],D[2][0][1],D[2][0][2],D[2][0][3]);
    
//     //iD41 = iD32 so delete this.
//     getIntDipFI(xint,mu2[1][0],muT2[1][0],D[2][1]);
//     printf("iD32   : % .10e % .10e % .10e % .10e\n",
//            D[2][1][0],D[2][1][1],D[2][1][2],D[2][1][3]);
    
    getIntDipFF(D[2][3]);
//     printf("iD34   : % .10e % .10e % .10e % .10e\n",
//            D[2][3][0],D[2][3][1],D[2][3][2],D[2][3][3]);
    
    getIntDipFI(xint,2*scalar[0][1],2*scalarT[0][1],mu2[0][1],muT2[0][1],D[3][0]);
//     printf("iD41   : % .10e % .10e % .10e % .10e\n",
//            D[3][0][0],D[3][0][1],D[3][0][2],D[3][0][3]);
    
//     //iD42 = iD31 so delete this.
//     getIntDipFI(xint,mu2[1][1],muT2[1][1],D[3][1]);
//     printf("iD42   : % .10e % .10e % .10e % .10e\n",
//            D[3][1][0],D[3][1][1],D[3][1][2],D[3][1][3]);
    
//     //iD43 = iD34 so delete this.
//     getIntDipFF(D[3][2]);
//     printf("iD43   : % .10e % .10e % .10e % .10e\n",
//            D[3][2][0],D[3][2][1],D[3][2][2],D[3][2][3]);
        
    // sum them up
//  calculate the charge correlated Dipoles.
//  IIRC it's incoming fermion or outgoing antifermion Q = +1 
//  incoming antifermion or outgoing fermion Q = -1. Multiplication yields the sign.
    double sum[3][4];
    for(i=0;i<(1+cEps*2);i++) { // counter for the poles part
        for(j=0;j<4;j++)
        // sum up the contributions for the rescaled momentum p1
        sum[i][j] = /*+D[0][0][i]=0.0*/
                    -D[0][1][i][j]
                    +D[0][2][i][j]
                    -D[0][3][i][j]
                    +D[2][0][i][j]
                    -D[3][0][i][j];
                    
//         sum up the contributions for the rescaled momentum p2 which is the same as for p1.
//         sum[1][i] = -/*D[1][0][i]=*/D[0][1][i]
// //                     D[1][1][i] = 0.
//                     -/*D[1][2][i] = */D[0][3][i]
//                     +/*D[1][3][i] = */D[0][2][i]
//                     -/*D[2][1][i] = */D[3][0][i]
//                     +/*D[3][1][i] = */D[2][0][i];
                    
//         sum[1][i] = -/*D[1][0][i]=*/D[0][1][i]
//                     /*D[1][1][i] = 0.*/
//                     -D[1][2][i]
//                     +D[1][3][i]
//                     -D[2][1][i]
//                     +D[3][1][i];
    }
    
    
//     printf("sum1   : % .10e % .10e % .10e % .10e\n",
//            sum[0],sum[1],sum[2],sum[3]);
//     printf("sum2   : % .10e % .10e % .10e % .10e\n",
//            sum[1][0],sum[1][1],sum[1][2],sum[1][3]);
//     printf("iD34   : % .10e % .10e % .10e % .10e\n",
//            2*D[2][3][0],2*D[2][3][1],2*D[2][3][2],2*D[2][3][3]);

    // delta and x=1 of plus part - bin this guy directly
    double delta[3],// regular and x plus part - bin this guy for p1 and p2 separately
           xpart[3];
    
    for(i=0;i<(1+2*cEps);i++) { // counter for the poles part
        delta[i] = alpha/(2*M_PI)*(2*sum[i][0]-2*sum[i][3]-2*D[2][3][i][0])*M1*jac1;
        xpart[i] = alpha/(2*M_PI)*(sum[i][1]+sum[i][2])*Mx*jacx;
        res[i] = -(delta[i]+2*xpart[i]);
    }

    binIntDip(delta[0], xpart[0], costh, costhp, pT, pTp, y, yp, wgt);
           
//     return -(delta+2*xpart);
}



