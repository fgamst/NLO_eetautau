/* unintegrated dipoles for Catani-Seymour factorization formula. */


/* unintegrated Dipoles */

void remap3to2II(double pIn[2][4], double pOut[3][4], double pInTilde[2][4],
                 double pOutTilde[2][4], int a, int i, int b, double xiab)
{
    int j, mu;
    double K[4], KT[4], KplKT[4]; // final state "total momentum"    
    
    for (mu=0;mu<4;mu++) 
    // remapping initial state and K momenta (Catani-Seymour 5.137,138,140).
    {
        pInTilde[a][mu] = pIn[a][mu]*xiab; // rescale the emitter momentum
        pInTilde[b][mu] = pIn[b][mu]; // the observer momentum stays the same
        K[mu] = pIn[a][mu] + pIn[b][mu] - pOut[i][mu]; // K^mu
        KT[mu] = pInTilde[a][mu] + pInTilde[b][mu]; // Ktilde^mu
        KplKT[mu] = K[mu] + KT[mu];
    }
    
    // remapping of final state (Catani Seymour 5.139)
    for (j=0;j<2;j++) for (mu=0;mu<4;mu++) 
        pOutTilde[j][mu] 
        = pOut[j][mu] - 2*metric(pOut[j],KplKT)/metric(KplKT,KplKT)*KplKT[mu] 
        + 2*metric(pOut[j],K)/metric(K,K)*KT[mu];
}

/* initial emitter and initial observer. See massless Catani-Seymour Sec 5.5

   indices:
   
   emitter  a
   photon   i
   observer b   */
double getDipII(double pIn[2][4], double pOut[3][4],
                double pInTilde[2][4],double pOutTilde[2][4],
                int a, int i, int b, int* res)
{
    int i1,i2;
    double papb = shat/2,//metric(pIn[a],pIn[b]), // scalar product of emitter and observer. For e^-e^+->tau^+tau^- always shat/2.
           papi = metric(pIn[a],pOut[i]); // scalar product of emitter and soft particle
    
    if(papi/papb > acut_II) {
        *res = 0;
        return 0.0;
    }// actual cut from alpha paper
//     if(papi/papb <= acut_II) return 0.0; // inverse cut for testing

    *res = 1;
    
    double pbpi = metric(pIn[b],pOut[i]);
           
    double xiab = 1.0 - (papi+pbpi)/papb;// rescaling factor.
           //(metric(pIn[a],pIn[b])-metric(pOut[i],pIn[a])-metric(pOut[i],pIn[b]))/metric(pIn[a],pIn[b]) // original version
           //K[4], KT[4], KplKT[4], // final state "total momentum"
    
    
    //remap the momenta from 3 to 2 particle phasespace
    remap3to2II(pIn, pOut, pInTilde, pOutTilde, a, i, b, xiab);
    
    // calculate the splitting function for fermion->fermion+boson (Catani-Seymour 5.145)
    double V = 8*M_PI*alpha*CF*(2/(1-xiab)-(1+xiab)/*-eps*(1-xiab)*/);//*muTo2Eps
    
    return -V*ME_Born_ee_tata(pInTilde,pOutTilde)/(2*metric(pIn[a],pOut[i]))/xiab;
}

//remap the 2 particle final state momenta. See 5.73
void remap3to2IF(double pIn[2][4], double pOut[3][4], double pInTilde[2][4],
                 double pOutTilde[2][4], int a, int i, int j, double xija)
{
    int mu; // Lorentz index
    
    for(mu=0;mu<4;mu++) {
        pInTilde[a][mu] = pIn[a][mu]*xija;
        pOutTilde[j][mu] = pOut[i][mu] + pOut[j][mu] - (1-xija)*pIn[a][mu];
        
        pInTilde[(a+1)%2][mu] = pIn[(a+1)%2][mu];
        pOutTilde[(j+1)%2][mu] = pOut[(j+1)%2][mu];
    }
}

/* initial emitter final observer. See massive Catani-Seymour Sec 5.3
   
   indices:
   
   a emitter
   i photon
   j observer
    TODO 
*/
double getDipIF(double pIn[2][4], double pOut[3][4],
                double pInTilde[2][4],double pOutTilde[2][4],
                int a, int i, int j, int* res)
{
    double papi = metric(pIn[a],pOut[i]),
           papj = metric(pIn[a],pOut[j]),
           ziT = papi/(papi+papj); // massive Catani-Seymour 5.74
           
    if (ziT > acut_IF) { // actual cutoff
        *res = 0;
        return 0.0;
    }
//     if (ziT <= acut_IF) return 0.0; // inverse cutoff for testing
    
    *res = 1;
    
    double pipj = metric(pOut[i],pOut[j]),
           xija = (papi+papj-pipj)/(papi+papj),
           zjT = papj/(papi+papj); // for xija, ziT, zjT see massive Catani-Seymour 5.74
    
    // remap the 2 particle final state momenta. See 5.73
    remap3to2IF(pIn, pOut, pInTilde,pOutTilde, a, i, j, xija);
//    for(i1=0;i1<4;i1++) {
//        pInTilde[a][i1] = pIn[a][i1]*xija;
//        pOutTilde[j][i1] = pOut[i][i1] + pOut[j][i1] - (1-xija)*pIn[a][i1];
        
//        pInTilde[(a+1)%2][i1] = pIn[(a+1)%2][i1];
//        pOutTilde[(j+1)%2][i1] = pOut[(j+1)%2][i1];
//    }
    
    // splitting function for f->f+b (massive Catani-Seymour 5.81)
    double V = 8*M_PI*alpha*CF*(2/(2-xija-zjT)-1-xija/*-eps*(1-xaij)*/)/**muTo2Eps*/;
    
    // return dipole from massive Catani-Seymour 5.71
    return -V*ME_Born_ee_tata(pInTilde,pOutTilde)/(2*papi)/xija;
}

/* final emitter initial observer. See massive Catani-Seymour Sec 5.2
   
   indices:
   
   i photon
   j emitter
   a observer
    TODO 
*/
double getDipFI(double pIn[2][4], double pOut[3][4],
                double pInTilde[2][4],double pOutTilde[2][4],
                int j, int i, int a, int* res)
{
    int i1;
    double papi = metric(pIn[a],pOut[i]),
           papj = metric(pIn[a],pOut[j]),
           pipj = metric(pOut[i],pOut[j]),
           xija = (papi+papj-pipj)/(papi+papj);
           
    if (xija < 1-acut_FI) {  // actual cutoff
        *res = 0;
        return 0.0;
    }
//     if (xija >= 1-acut_FI) return 0.0; // inverse cutoff for testing
    
    *res = 1;
    
    double //ziT = papi/(papi+papj),
           zjT = papj/(papi+papj);
           // for xaij, ziT, zjT see massive Catani-Seymour 5.42
    
    // remap the 2 particle final state momenta. See 5.43
    for(i1=0;i1<4;i1++) {
        pInTilde[a][i1] = pIn[a][i1]*xija;
        pOutTilde[j][i1] = pOut[i][i1] + pOut[j][i1] - (1-xija)*pIn[a][i1];
        
        
        pInTilde[(a+1)%2][i1] = pIn[(a+1)%2][i1];
        pOutTilde[(j+1)%2][i1] = pOut[(j+1)%2][i1];
    }
    
    // splitting function for f->f+b (massive Catani-Seymour 5.50)
    double V = 8*M_PI*alpha*CF*(2/(2-xija-zjT)-1-zjT-MT2/pipj/*-eps*(1-zjT)*/);
    
    // return dipole from massive Catani-Seymour 5.40
    return -V*ME_Born_ee_tata(pInTilde,pOutTilde)/(2*pipj)/xija;
}

/* final emitter final observer. See massive Catani-Seymour Sec 5.1
   
   indices:
   
   i emitter
   j photon
   k observer
   D^a_{ij}
 *works with ziT instead of zjT - compare massless with massive Catani-Seymour */
double getDipFF(double pIn[2][4], double pOut[3][4],
                double pInTilde[2][4],double pOutTilde[2][4],
                int i, int j, int k, int* res)
{
    int i1;
    double //muTo2Eps = 1.0, // mu^2eps regulator
           Q[4],/*Qc[4],*/ pij[4],pij2, // final state "total momentum", composite particle momentum
           //mij,mui,muj,muk,muij, // composite mass and rescaled masses
           pipj = metric(pOut[i],pOut[j]),
           pipk = metric(pOut[i],pOut[k]),
           pjpk = metric(pOut[j],pOut[k]),
           yijk = pipj/(pipj+pipk+pjpk);
           //metric(pOut[i],pOut[j])
           //     /(metric(pOut[i],pOut[j])+metric(pOut[i],pOut[k])+metric(pOut[j],pOut[k]));

    
    if (yijk > acut_FFp) { // actual cut
        *res = 0;
        return 0.0;
    }
//     if (yijk <= acut_FFp) return 0.0; // inverse cut for testing
    
    *res = 1;
    
    double ziT = pipk/(pipk+pjpk),
//            zjT = 1.0 - ziT,
           // for yijk, ziT, zjT see Catani-Seymour 5.12
           vijk,vijkT;//,vijkc,vijkc2,vijkTc; // relative velocity, rescaled relative velocity.
    
    // calculate the momenta.
    // ATTENTION 
    // in general Q is not {Etot,0,0,0}, only in case of 3 final state particles
    for(i1=0;i1<4;i1++) {
        Q[i1] = 0.0;
//         Qc[i1] = pOut[0][i1] + pOut[1][i1] + pOut[2][i1];
        pij[i1] = pOut[i][i1] + pOut[j][i1];
    }
    Q[0] = Etot;
    
    pij2 = metric (pij,pij);
    
    // calculate the rescaled momenta. See 5.9
    // ATTENTION 
    // only for Q = {Etot,0,0,0}, mi=mk=mij=MT, mj=0 and 3 particle final state.
    for(i1=0;i1<4;i1++) {
        pOutTilde[k][i1] = sqrt(shat*shat - 4*shat*MT2)
        /sqrt(shat*shat + pij2*pij2 + MT2*MT2 - 2*shat*pij2 - 2*shat*MT2 - 2*pij2*MT2)
        *(pOut[k][i1] - Q[i1]*pOut[k][0]*Etot/shat)
        +Q[i1]/2;
        
        pOutTilde[i][i1] = Q[i1] - pOutTilde[k][i1];
        pInTilde[0][i1] = pIn[0][i1];
        pInTilde[1][i1] = pIn[1][i1];
    }
    
//     vijkc = sqrt(1-(metric(pij,pij)*metric(pOut[k],pOut[k]))
//              /(metric(pij,pOut[k])*metric(pOut[k],pij)));
//     vijkTc = sqrt(1-(metric(pOutTilde[i],pOutTilde[i])*metric(pOutTilde[k],pOutTilde[k]))
//             /(metric(pOutTilde[i],pOutTilde[k])*metric(pOutTilde[k],pOutTilde[i])));
    
    // calculate the masses. Massive Catani-Seymour 5.5
    // ATTENTION 
    // again these only hold only for a 3 final state tau^+ tau^- photon
//     mij = MT;
//     mui = MT/Etot;
//     muj = 0;
//     muij = mui;
//     muk = mui;
    // maybe implement 5.8 for v then, since it seems more advantageous.
    double pijpk = metric(pij,pOut[k]);
    // calculate the relative velocities. 5.14, 5.8
    vijk = sqrt(1-pij2*MT2/pijpk/pijpk); // massive Catani-Seymour 5.5
    vijkT = sqrt(1-4*MUT*MUT)/(1-2*MUT*MUT); // see the notebook int_dipoles.nb
//     sqrt(1.0 + muij*muij*muij*muij + muk*muk*muk*muk - 2*muij*muij - 2*muk*muk 
//             - 2*muij*muij*muk*muk)/(1.0 - muij*muij - muk*muk);
//     vijkc2 = sqrt((2*muk*muk + (1.0 - mui*mui /*- muj*muj*/ - muk*muk)*(1-yijk))
//            *(2*muk*muk + (1.0 - mui*mui /*- muj*muj*/ - muk*muk)*(1-yijk))-4*muk*muk)
//            /(1.0 - mui*mui /*- muj*muj*/ - muk*muk)/(1-yijk);
//     printf("vijkT   : % .10e\n",vijkT);
//     printf("vijkTc  : % .10e\n",vijkTc);
//     printf("vijk   : % .16e\n",vijk);
//     printf("vijkc  : % .16e\n",vijkc);
//     printf("vijkc2 : % .16e\n",vijkc2);
    
    // calculate the splitting function for fermion->fermion+boson (Catani-Seymour 5.16)
    // ATTENTION
    // actually inserting ziT gives the correct dipole compare massless 5.16 with massive
    // 5.6 Catani-Seymour. This might be an error in the paper!
    double V = 8*M_PI*/*muTo2Eps**/alpha*/*CF**/(2/(1-ziT*(1-yijk))
             -vijkT/vijk*(1+ziT+MT2/pipj/*-eps*(1-ziT)*/));
    
    return -V*ME_Born_ee_tata(pIn,pOutTilde)/(2*pipj);
}

double getDipSumOld(double pIn[2][4], double pOut[3][4],double jacobian,double weight) {
    int i,j,res;
    double pInTilde[2][4], pOutTilde[2][4], D[4][4], sum = 0.0;
                
//  calculate the charge correlated Dipoles.
//  Some details on why which sign goes where would be nice.
//  IIRC it's incoming fermion or outgoing antifermion Q = +1 
//  incoming antifermion or outgoing fermion Q = -1. Multiplication yields the sign.
    D[0][0] = 0.0;
    D[0][1] = 0.0-getDipII(pIn,pOut,pInTilde,pOutTilde,0,2,1,&res);
    if(res) binDip(D[0][1],pInTilde,pOutTilde,jacobian,weight);
    
    D[0][2] = 0.0+getDipIF(pIn,pOut,pInTilde,pOutTilde,0,2,0,&res);
    if(res) binDip(D[0][2],pInTilde,pOutTilde,jacobian,weight);
    
    D[0][3] = 0.0-getDipIF(pIn,pOut,pInTilde,pOutTilde,0,2,1,&res);
    if(res) binDip(D[0][3],pInTilde,pOutTilde,jacobian,weight);
    
    D[1][0] = 0.0-getDipII(pIn,pOut,pInTilde,pOutTilde,1,2,0,&res);
    if(res) binDip(D[1][0],pInTilde,pOutTilde,jacobian,weight);
    
    D[1][1] = 0.0;
    D[1][2] = 0.0-getDipIF(pIn,pOut,pInTilde,pOutTilde,1,2,0,&res);
    if(res) binDip(D[1][2],pInTilde,pOutTilde,jacobian,weight);
    
    D[1][3] = 0.0+getDipIF(pIn,pOut,pInTilde,pOutTilde,1,2,1,&res);
    if(res) binDip(D[1][3],pInTilde,pOutTilde,jacobian,weight);
    
    D[2][0] = 0.0+getDipFI(pIn,pOut,pInTilde,pOutTilde,0,2,0,&res);
    if(res) binDip(D[2][0],pInTilde,pOutTilde,jacobian,weight);
    
    D[2][1] = 0.0-getDipFI(pIn,pOut,pInTilde,pOutTilde,0,2,1,&res);
    if(res) binDip(D[2][1],pInTilde,pOutTilde,jacobian,weight);
    
    D[2][2] = 0.0;
    D[2][3] = 0.0-getDipFF(pIn,pOut,pInTilde,pOutTilde,0,2,1,&res);
    if(res) binDip(D[2][3],pInTilde,pOutTilde,jacobian,weight);
    
    D[3][0] = 0.0-getDipFI(pIn,pOut,pInTilde,pOutTilde,1,2,0,&res);
    if(res) binDip(D[3][0],pInTilde,pOutTilde,jacobian,weight);
    
    D[3][1] = 0.0+getDipFI(pIn,pOut,pInTilde,pOutTilde,1,2,1,&res);
    if(res) binDip(D[3][1],pInTilde,pOutTilde,jacobian,weight);
    
    D[3][2] = 0.0-getDipFF(pIn,pOut,pInTilde,pOutTilde,1,2,0,&res);
    if(res) binDip(D[3][2],pInTilde,pOutTilde,jacobian,weight);
    
    D[3][3] = 0.0;
    
//  sum them up
    for(i=0;i<4;i++) {
        for(j=0;j<4;j++) sum += D[i][j];
    }
    
    return sum;
}

/* sign from charge correlations
   index    particle    sign contribution 
   0        e^-         +1
   1        e^+         -1
   2        tau^+       +1
   3        tau^-       -1
   
   */


// function that decides from the indices which dipole to call for our process
// double getDip(double pIn[2][4], double pOut[3][4],
//               double pInTilde[2][4], double pOutTilde[2][4],
//               int em, int ph, int sp) {
//                 
// //  calculate the charge correlated Dipoles.
// //  Some details on why which sign goes where would be nice.
// //  IIRC it's incoming fermion or outgoing antifermion Q = +1 
// //  incoming antifermion or outgoing fermion Q = -1. Multiplying yields the sign.
//     if (em<2) {
//         if (sp<2) return 0.0;/*getDipII(pIn,pOut,pInTilde,pOutTilde,em,ph-2,sp);*/
//         else return 0.0;/*getDipIF(pIn,pOut,pInTilde,pOutTilde,em,ph-2,sp-2);*/
//     }
//     else {
//         em -= 2;
//         if (sp<2) return 0.0;/*getDipFI(pIn,pOut,pInTilde,pOutTilde,em-2,ph-2,sp);*/
//         else return getDipFF(pIn,pOut,pInTilde,pOutTilde,em,ph-2,sp-2);
//     }
// }

// double getDipSum(double pIn[2][4], double pOut[3][4],double jacobian,double weight) {
//     int i,j,k,Qc;
//     double pInTilde[2][4], pOutTilde[2][4], D[4][4], sum = 0.0,costh;
//                 
// //  calculate the charge correlated Dipoles.
// //  Some details on why which sign goes where would be nice.
// //  IIRC it's incoming fermion or outgoing antifermion Q = +1 
// //  incoming antifermion or outgoing fermion Q = -1. Multiplying yields the sign.
//     
//     for(i=0;i<4;i++) {
//         for(j=0;j<4;j++) {
//             pInTilde[i][j] = 0.0;
//             pOutTilde[i][j] = 0.0;
//         }
//     }

//     for(i=0;i<4;i++) {
//         for(j=0;j<4;j++) {
//             if (i!=j) {
//                 // calculate the charge correlations.
//                 // ATTENTION
//                 // This is only valid for e- e+ -> ta+ ta- (A) parametrisation
//                 if ((i+j) % 2 == 0) Qc = 1; 
//                 else Qc = -1;
//                 D[i][j] = Qc*getDip(pIn,pOut,pInTilde,pOutTilde,i,4,j);
//                 costh = dot3(pInTilde[0],pOutTilde[bin])
//                         /abs3(pInTilde[0])/abs3(pOutTilde[bin]);
//                 k = (int)((costh+1)/2*(NBINS));
//                 fbin[0][k] -= (jacobian * D[i][j] * weight);
//                 sum += D[i][j]; // sum them up
//             }
//             else D[i][j] = 0.0;
//         }
//     }
//     for(i=0;i<4;i++) {
//         for(j=0;j<4;j++) {
//             if (i!=j) {
//                 // calculate the charge correlations.
//                 // ATTENTION
//                 // This is only valid for e- e+ -> ta+ ta- (A) parametrisation
//                 if ((i+j) % 2 == 0) Qc = 1; 
//                 else Qc = -1;
//                 if (i<2) {
//                     if (j<2)
//                     D[i][j] = Qc*getDipII(pIn,pOut,pInTilde,pOutTilde,i,2,j);
//                     
//                     else 
//                     D[i][j] = Qc*getDipIF(pIn,pOut,pInTilde,pOutTilde,i,2,j-2);
//                 }
//                 else {
//                     if (j<2)
//                     D[i][j] = Qc*getDipFI(pIn,pOut,pInTilde,pOutTilde,i-2,2,j);
//                     
//                     else 
//                     D[i][j] = Qc*getDipFF(pIn,pOut,pInTilde,pOutTilde,i-2,2,j-2);
//                 }
//                 sum += D[i][j]; // sum them up
//             }
//         }
//     }
//    
//    return sum;
// }
