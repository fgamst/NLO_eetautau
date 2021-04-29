// bin the tree level matrix element.
void binM2(double M2, double pIn[2][4], double pOut[2][4], double wgt, double jac) {
    int i;
    double costh, pT, y,/* w = wgt/(NIT), w2 = w*w,*/ fx = (jac * M2 * wgt), fx2 = fx*fx;
    
    // Bin M3 according to M2 kinematic
    // Calculate the equivalent of cos(theta) in 2 particle phase space.
    costh = pOut[bin][3]/abs3(pOut[bin]);
    i = (int)((NBINS)*(costh+1.0)/2.0);
    itBin[countIt][0][i] += fx;
//     dbin[0][i] += fx2;
//     dbin[0][i] += w2;
//     printf("costh = % .6f\n",costh);
//     printf("i     = %d\n",i);
    
    // bin according to transverse momentum
    pT = sqrt(pOut[bin][1]*pOut[bin][1]+pOut[bin][2]*pOut[bin][2]);
    i = (int)(pT/absp3*(NBINS));
    itBin[countIt][1][i] += fx;
//     dbin[1][i] += fx2;
//     dbin[1][i] += w2;
//     printf("pT    = % .6f\n",pT);
//     printf("i     = %d\n",i);
    
    // bin according to rapidity
//     y = 0.5*log((pOut[bin][0]+pOut[bin][3])/(pOut[bin][0]-pOut[bin][3]));
//     i = (int)(y/ymax*(NBINS));
//     fbin[2][i] += fx;
//     dbin[2][i] += fx2;
    return;
}

// bin the real correction matrix element.
void binM3(double M3, double pIn[2][4], double pOut[3][4], double wgt, double jac) {
    int i;
    double/* idouble, */costh, pT, y,/* w = wgt/(NIT), w2 = w*w,*/ fx = (jac * M3 * wgt), fx2 = fx*fx;
    
    
    // Bin M3 according to M2 kinematic
    // Calculate the equivalent of cos(theta) in 2 particle phase space.
    costh = pOut[bin][3]/abs3(pOut[bin]);
    i = (int)((NBINS)*(costh+1.0)/2.0);
    itBin[countIt][0][i] += fx;
//     dbin[0][i] += fx2;
//     dbin[0][i] += w2;
//     printf("costh = % .6f\n",costh);
//     printf("i     = %d\n",i);
    
    // bin according to transverse momentum
    pT = sqrt(pOut[bin][1]*pOut[bin][1]+pOut[bin][2]*pOut[bin][2]);
    i = (int)(pT/absp3*(NBINS));
    itBin[countIt][1][i] += fx;
//     dbin[1][i] += fx2; 
//     dbin[1][i] += w2;
//     printf("pT    = % .6f\n",pT);
//     printf("i     = %d\n",i);
    
    // bin according to rapidity
//     y = 0.5*log((pOut[bin][0]+pOut[bin][3])/(pOut[bin][0]-pOut[bin][3]));
//     i = (int)(y/ymax*(NBINS));
//     fbin[2][i] += fx;
//     dbin[2][i] += fx2;
    return;
}

//bin the remapped tilde momenta
void binDip(double Dip, double pIn[2][4], double pOut[2][4],double wgt, double jac) {
    int i;
    double costh, pT, y,/* w = wgt/(NIT), w2 = w*w,*/ fx = (jac * Dip * wgt), fx2 = fx*fx;
//     if(abs(pOut[0][3]+pOut[1][3])>1e-14) {
//         double beta = (pOut[0][3]+pOut[1][3])/(pOut[0][0]+pOut[1][0]),
//                gamma = 1/sqrt((1-beta)*(1+beta));
// //         printf("beta  = % .10f\n",beta);
// //         printf("gamma = % .10f\n",gamma);
// //         // boost to the center of mass
//         double p3Til = (pOut[bin][3]-beta*pOut[bin][0])*gamma;
// //         printf("p3'   = % .10f\n",p3Til);
//         costh = p3Til/sqrt(pOut[bin][1]*pOut[bin][1]+pOut[bin][2]*pOut[bin][2]
//                                      +p3Til*p3Til);
// //         printf("costh = % .10f\n",costh);
//     }
//     else costh = dot3(pIn[0],pOut[bin])/abs3(pIn[0])/abs3(pOut[bin]);
    costh = dot3(pIn[0],pOut[bin])/abs3(pIn[0])/abs3(pOut[bin]);
    // bin according to costh
    i = (int)((costh+1)/2*(NBINS));
    itBin[countIt][0][i] -= fx;
//     dbin[0][i] += fx2;
//     dbin[0][i] += w2;
//     printf("costh = % .6f\n",costh);
//     printf("i     = %d\n",i);
    
    // bin according to transverse momentum
    pT = sqrt(pOut[bin][1]*pOut[bin][1]+pOut[bin][2]*pOut[bin][2]);
    i = (int)((pT)/absp3*(NBINS));
    itBin[countIt][1][i] -= fx;
//     dbin[1][i] += fx2;
//     dbin[1][i] += w2;
//     printf("pT    = % .6f\n",pT);
//     printf("i     = %d\n",i);
    
    // bin according to rapidity
//     y = 0.5*log((pOut[bin][0]+pOut[bin][3])/(pOut[bin][0]-pOut[bin][3]));
//     i = (int)(y/ymax*(NBINS));
//     fbin[2][i] += fx;
//     dbin[2][i] += fx2;
    return;
}

// bin the integrated dipoles 
void binIntDip(double delta, double xpart,
               double costh, double costhp[2], double pT, double pTp,
               double y, double yp,double wgt) {
    int i, j;
    double //w = wgt / (NIT), w2 = w*w,
           del = delta * wgt,
           del2 = del*del,
           xp = xpart * wgt,
           xp2 = xp*xp;
    //  bin the delta part just at cos(theta).
    i = (int)((costh+1)/2*(NBINS));
    itBin[countIt][0][i] -= del;
//     dbin[0][i] += del2;
//     dbin[0][i] += w2;
//     dbin[0][i] += fb
//     printf("costh  = % .6f\n",costh);
//     printf("i      = %d\n",i);
        
    //  bin the rescaled part at the respective cos(theta_prime).
    for(j=0;j<2;j++) {
        i = (int)((costhp[j]+1)/2*(NBINS));
        itBin[countIt][0][i] -= xp; 
//         dbin[0][i] += xp2;
//         dbin[0][i] += w2;
//     printf("costhp%d= % .6f\n",j+1,costhp[j]);
//     printf("i      = %d\n",i);
    }
    
    // bin delta part according to transverse momentum
    i = (int)(pT/absp3*(NBINS));
    itBin[countIt][1][i] -= del;
//     dbin[1][i] += del2;
//     dbin[1][i] +=w2;
//     printf("pT     = % .6f\n",pT);
//     printf("i      = %d\n",i);
    
    // bin x dependent part according to transverse momentum
    i = (int)(pTp/absp3*(NBINS));
    itBin[countIt][1][i] -= 2*xp;
//     dbin[1][i] += 4*xp2;
//     dbin[1][i] +=4*w2;
//     printf("pTp    = % .6f\n",pTp);
//     printf("i      = %d\n",i);
    
    // bin according to rapidity
//     i = (int)(y/ymax*(NBINS));
//     fbin[2][i] += (delta * wgt / (NIT));
//     dbin[2][i] += del2;
//     
//     // bin x dependent part according to rapidity
//     i = (int)(yp/ymax*(NBINS));
//     fbin[2][i] += ((2*xpart) * wgt / (NIT));
//     dbin[2][i] += 4*xp2;
    return;
}

// Turans empfehlung für die Fehler.
// res[0]=(this->*integrand)(x,&wgt); 
// 
//   if ( fillHistograms ) {
//     const double w=res[0]*wgt/integrator->getcalls()/it;
//     fillOwnHistograms( w, w*w );
//     logEvent(x, w, res[0]);
//   }
//   return(res[0]);
// am ende die wurzel ziehen.
// scrap that we are doing standard deviation
