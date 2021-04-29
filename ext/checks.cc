// Check the tree level squared matrix element N times.
void checkM2(int N) {
    int i, j;
    double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}}; 
    double pOut[2][4], vars[2], x[2], jac, M2, M2_check;
    
        // call the helas coupling constants of the standard model.
//     i = 0; // i = 0 is default configuration
//     helascoupsm(&i);
    
    
    for(i=0;i<N;i++) {
        for(j=0;j<2;j++) x[j]=((double)rand())/((double) RAND_MAX);
        PhSp_TaTa(Etot,x,vars,pOut,&jac);
        M2 = ME_Born_ee_tata(pIn, pOut);
        M2_check = ME_Born_ee_tata/*_check*/(pIn, pOut);
        printf("\n");
        for(j=0;j<2;j++) printf("p%d    : % .10e % .10e % .10e % .10e\n",
                            j+3,pOut[j][0],pOut[j][1],pOut[j][2],pOut[j][3]);
        printf("M2       =  % .16e \n", M2);
        printf("M2_check =  % .16e \n", M2_check);
//         double oCanc = -log10(abs((M2-M2_check)/M2));
//         printf("oCanc    =  % .16e\n",oCanc);
    }
    
    return;
}

// Checks the tree level squared matrix element including a photon and a Z mediated process.
void checkM2gZ(int N) {
    int i, j;
    double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}}; 
    double pOut[2][4], vars[2], x[2], jac, M2, M2_check;
    
        // call the helas coupling constants of the standard model.
//     i = 0; // i = 0 is default configuration
//     helascoupsm(&i);
    
    // if we make some changes to gp and gm we can get closer to what madgraph+helas get
//     gp += -0.000202;
//     gm += 0.00026;
    // i dunno whats happening here cw^2 and sw^2 are definitely the same
    
    printf("sin^2(W) =  % .16f \n", sw2);
    printf("cos^2(W) =  % .16f \n", cw2);
    printf("gp       =  % .16f \n", gp);
    printf("gm       =  % .16f \n", gm);
    printf("af       =  % .16f \n", af);
    printf("vf       =  % .16f \n", vf);
    
    for(i=0;i<N;i++) {
        for(j=0;j<2;j++) x[j]=((double)rand())/((double) RAND_MAX);
        PhSp_TaTa(Etot,x,vars,pOut,&jac);
//         M2 = ME_Born_ee_tata(pIn, pOut);
        M2 = ME_Born_ee_AZ_tata(pIn, pOut);
        M2_check = ME_Born_ee_AZ_tata_check(pIn, pOut);
        printf("\n");
        for(j=0;j<2;j++) printf("p%d    : % .10e % .10e % .10e % .10e\n",
                            j+3,pOut[j][0],pOut[j][1],pOut[j][2],pOut[j][3]);
        printf("M2       =  % .16e \n", M2);
        printf("M2_check =  % .16e \n", M2_check);
//         printf("ratio    =  % .16e\n",M2_check/M2);
        double oCanc = -log10(abs((M2-M2_check)/M2));
        printf("oCanc    =  % .2f\n",oCanc);
    }
    
    return;
}

// Makes N testruns to check three particle phase space. 
void checkPS3(int N) {
    int i, j;
    double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}}; 
    double pOut[3][4], vars[2], x[5], jac;
    
    
    
    for(i=0;i<N;i++) {
        for(j=0;j<5;j++) x[j]=((double)rand())/((double) RAND_MAX);
        PhSp_TaTaPh(Etot,x,vars,pOut,&jac);
        printf("\n");
        for(j=0;j<3;j++) {
            printf("p%d    : % .10e % .10e % .10e % .10e\n",
                            j+3,pOut[j][0],pOut[j][1],pOut[j][2],pOut[j][3]);
            printf("abs3(p%d)    =  % .16e \n",j+3,abs3(pOut[j])/absp3);
            if(abs3(pOut[j])/absp3<0.00001/*>0.9999999*/) {
                printf("ERROR: abs3(p%d) too large\n",j+3);
                return;
            }
        }
    }
    
    return;
}

// Checks the real correction squared matrix element N times.
void checkM3(int N) {
    int i, j;
    double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}}; 
    double pOut[3][4], vars[2], x[5], jac, M3, M3_check;
    
        // call the helas coupling constants of the standard model.
//     i = 0; // i = 0 is default configuration
//     helascoupsm(&i);
    
    
    for(i=0;i<N;i++) {
        for(j=0;j<5;j++) x[j]=((double)rand())/((double) RAND_MAX);
        PhSp_TaTaPh(Etot,x,vars,pOut,&jac);
        M3 = ME_ee_tataph(pIn, pOut);
        M3_check = ME_ee_tataph_check(pIn, pOut);
        printf("\n");
        for(j=0;j<3;j++) printf("p%d    : % .10e % .10e % .10e % .10e\n",
                            j+3,pOut[j][0],pOut[j][1],pOut[j][2],pOut[j][3]);
        printf("M3       =  % .16e \n", M3);
        printf("M3_check =  % .16e \n", M3_check);
//         printf("ratio    =  % .16e\n",M3_check/M3);
        double oCanc = -log10(abs((M3-M3_check)/M3));
        printf("oCanc    =  % .2f\n",oCanc);
    }
    
    return;
}

// check the cancellation of squared matrix element and dipoles in the singular limits
void checkDip(int N) {
    int i, j, col1, col2, NPart=3, Steps = N;
    double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}}; 
    double Masses[3] = {MT, MT, 0.0}; 
    double pOut[3][4], vars[2], x[5], M3, Dip, CMSEnergy = Etot, Jacobian, oCanc, weight = 1.0, SingDepth = 1e-16;

    
    acut_II = 1.0;
    acut_IF = 1.0;
    acut_FI = 1.0;
    acut_FF = 1.0;
    
    initIntDipoles();
    
//     //  p5||p1
//     col1=1 -1;
//     col2=5 -1;
    
//     //  p5||p2   
    col1=2 -1;
    col2=5 -1;
    
//     //  p5||p3
//     col1=1 -1;
//     col2=5 -1;
    
//     //  p5||p4   
//     col1=4 -1;
//     col2=5 -1;
    
//     //  p5 soft
//     col1=5 -1;
//     col2=5 -1;
    
    for( i = 1; i <= N; i++){    
//  generate singular phase space:  
        printf("Step %d\n", i);
        gensing_(&NPart,&CMSEnergy,Masses,pOut,&col1,&col2,&SingDepth,&Steps);

//  print output 
        printf("Mom1: % .10e % .10e % .10e % .10e \n",
           pOut[0][0],pOut[0][1],pOut[0][2],pOut[0][3]);
        printf("Mom2: % .10e % .10e % .10e % .10e \n",
           pOut[1][0],pOut[1][1],pOut[1][2],pOut[1][3]);
        printf("Mom3: % .10e % .10e % .10e % .10e \n",
           pOut[2][0],pOut[2][1],pOut[2][2],pOut[2][3]);
        printf("Sum : % .10e % .10e % .10e % .10e \n",
           pOut[0][0]+pOut[1][0]+pOut[2][0],
           pOut[0][1]+pOut[1][1]+pOut[2][1],
           pOut[0][2]+pOut[1][2]+pOut[2][2],
           pOut[0][3]+pOut[1][3]+pOut[2][3]);
        printf("\n\n");

        printf("p1.p5 : %e \n",lsp(pIn[0],pOut[2])/pIn[0][0]/pOut[2][0]);
        printf("p2.p5 : %e \n",lsp(pIn[1],pOut[2])/pIn[1][0]/pOut[2][0]);
        printf("p3.p5 : %e \n",lsp(pOut[0],pOut[2])/pOut[0][0]/pOut[2][0]);
        printf("p4.p5 : %e \n",lsp(pOut[1],pOut[2])/pOut[1][0]/pOut[2][0]);
        printf("E5/E  : %e \n",pOut[2][0]*pOut[2][0]/shat);
    
        M3 = ME_ee_tataph(pIn,pOut);
        Dip = getDipSumOld(pIn,pOut,Jacobian,weight);
        
        printf("M3      =  % .16e \n", M3);
        printf("Dipole  =  % .16e \n", Dip);
//         printf("ratio    =  % .16e\n",Dip/M3);
        oCanc = -log10(abs((M3-Dip)/M3));
        printf("oCanc   =  % .2f\n",oCanc);
    
//         getchar();
    
        printf("\n\n");
//         printf("\n\n\n");
    }
    return;
}

// Checks the real correction integrand for numeric stability in singular limits.
void checkRealCorIntegrand(int N) {
    int i, j;
    double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}}; 
    double pOut[3][4], vars[2], x[5], jac, M3, M3_check, *weight, Dips, result[1];
    
    *weight = 1.0;
        // call the helas coupling constants of the standard model.
//     i = 0; // i = 0 is default configuration
//     helascoupsm(&i);
    
    
    for(i=0;i<N;i++) {
        for(j=0;j<5;j++) x[j]=((double)rand())/((double) RAND_MAX);
        PhSp_TaTaPh(Etot,x,vars,pOut,&jac);
        printf("Point %d\n", i+1);
        M3 = ME_ee_tataph(pIn, pOut);
        M3_check = ME_ee_tataph_check(pIn, pOut);
        Dips = getDipSumOld(pIn,pOut,jac,*weight);
        realCorIntegrand(x, weight, result);
        printf("\n");
        for(j=0;j<3;j++) printf("p%d    : % .10e % .10e % .10e % .10e\n",
                            j+3,pOut[j][0],pOut[j][1],pOut[j][2],pOut[j][3]);
        printf("jacobian =  % .16e \n", jac);
        printf("M3       =  % .16e \n", M3);
        printf("Dips     =  % .16e \n", Dips);
//         printf("M3_check =  % .16e \n", M3_check);
        printf("check    =  % .16e \n", (M3-Dips)*jac);
        printf("result   =  % .16e \n", result[0]);
//         printf("ratio    =  % .16e\n",M3_check/M3);
//         double oCanc = -log10(abs((M3-M3_check)/M3));
//         printf("oCanc    =  % .2f\n",oCanc);
    }
    
    return;
}

// check M3 bins for functionality.
void checkM3Bin(int N) {
    int i, j;
    double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}}; 
    double pOut[3][4], x[5], vars[5], wgt = 1.0, jac, M3;
    
    for(i=0;i<N;i++) {
        for(j=0;j<5;j++) x[j]=((double)rand())/((double) RAND_MAX);
        PhSp_TaTaPh(Etot,x,vars,pOut,&jac);
        M3 = ME_ee_tataph(pIn, pOut);
        
        printf("\n");
        printf("p%d    : % .10e % .10e % .10e % .10e\n",
                            bin+3,pOut[bin][0],pOut[bin][1],pOut[bin][2],pOut[bin][3]);
        binM3(M3,pIn, pOut,wgt,jac);
    }
    
    return;
}

// check the dipole bins for functionality.
void checkDipoleBin(int N) {
    int i, j, mu, res;
    double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}}; 
    double pOut[3][4], pInTilde[2][4], pOutTilde[2][4],x[5], vars[5], wgt = 1.0, jac,
           Dips, beta, gamma, p3Til, costh;
    
    double /*mui = 0.0,*/ muj = MUT, muj2 = muj*muj, muk = MUT, muk2 = muk*muk,
                yp = 1 - 2*muk*(1-muk)/(1-muj2-muk2);// upper limit of the yijk integral
                
    acut_FFp = acut_FF*yp;
    
    for(i=0;i<N;i++) {
        for(j=0;j<5;j++) x[j]=((double)rand())/((double) RAND_MAX);
        PhSp_TaTaPh(Etot,x,vars,pOut,&jac);
        printf("\n");
        
        // final state momenta
        for(j=0;j<3;j++){
            printf("p%d   :",j+3);
            for(mu=0;mu<4;mu++)printf(" % .10e",pOut[j][mu]);
            printf("\n");
        }
        printf("\n");
//         Dips = getDipII(pIn,pOut,pInTilde,pOutTilde,0,2,1,&res);
//         Dips = getDipIF(pIn,pOut,pInTilde,pOutTilde,0,2,1,&res);
//         Dips = getDipFI(pIn,pOut,pInTilde,pOutTilde,0,2,1,&res);
        Dips = getDipFF(pIn,pOut,pInTilde,pOutTilde,0,2,1,&res);
        if(res)
        {
            // print remapped momenta
            for(j=0;j<2;j++){
                printf("pT%d   :",j+1);
                for(mu=0;mu<4;mu++)printf(" % .10e",pInTilde[j][mu]);
                printf("\n");
                printf("pT%d^2 : % .10e",j+1,metric(pInTilde[j],pInTilde[j]));
                printf("\n");
            }
            printf("Sum   :");
            for(mu=0;mu<4;mu++)printf(" % .10e",pInTilde[0][mu]+pInTilde[1][mu]);
            printf("\n\n");
        
            for(j=0;j<2;j++){
                printf("pT%d   :",j+3);
                for(mu=0;mu<4;mu++)printf(" % .10e",pOutTilde[j][mu]);
                printf("\n");
                printf("pT%d^2 : % .10e",j+3,metric(pOutTilde[j],pOutTilde[j]));
                printf("\n");
            }
            printf("Sum   :");
            for(mu=0;mu<4;mu++)printf(" % .10e",pOutTilde[0][mu]+pOutTilde[1][mu]);
            printf("\n\n");
        
            binDip(Dips,pInTilde,pOutTilde,wgt,jac);
        
            beta = (pInTilde[0][3]+pInTilde[1][3])/(pInTilde[0][0]+pInTilde[1][0]);
            gamma = 1/sqrt((1-beta)*(1+beta));
            printf("beta  = % .10f\n",beta);
            printf("gamma = % .10f\n",gamma);
            // boost to the center of mass
            p3Til = (pOutTilde[bin][3]-beta*pOutTilde[bin][0])*gamma;
            printf("p3'   = % .10f\n",p3Til);
            costh = p3Til/sqrt(pOutTilde[bin][1]*pOutTilde[bin][1]
                                     +pOutTilde[bin][2]*pOutTilde[bin][2]
                                     +p3Til*p3Til);
        
            printf("costh = % .10f\n",costh);
        }
        else printf("Outside alpha cut\n");
        
    }
    
    return;
}

// check the binning of the integrated Dipoles as well as the lorentz boost to the pure rescaled system.
void checkIntDipoleBin(int N) {
    int i, j;
    double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}}; 
    double pOut[3][4], pInTilde[2][4], pOutTilde[2][4],x[3], wgt = 1.0, Dips[3];
    
    initIntDipoles();
    
    for(i=0;i<N;i++) {
//         for(j=0;j<3;j++) x[j]=((double)rand())/((double) RAND_MAX);
        for(j=0;j<2;j++) x[j]=((double)rand())/((double) RAND_MAX);
        x[2]= 0.00000000000001*(i+99999999999991);
        printf("\n");
        
        if (x[2]*shat < 4*(MT2)) printf("xint     =  % .14f => below threshold\n",x[2]);
        else printf("xint     =  % .14f => above threshold\n",x[2]);
    
        getIntDipSum(x,pIn,pOut,pInTilde,pOutTilde,wgt,Dips);
        
//         printf("pT1   : % .10e % .10e % .10e % .10e\n",
//                 pInTilde[0][0],pInTilde[0][1],pInTilde[0][2],pInTilde[0][3]);
//         printf("pT2   : % .10e % .10e % .10e % .10e\n",
//                 pInTilde[1][0],pInTilde[1][1],pInTilde[1][2],pInTilde[1][3]);
//         printf("p3    : % .10e % .10e % .10e % .10e\n",
//                 pOut[0][0],pOut[0][1],pOut[0][2],pOut[0][3]);
//         printf("pT3   : % .10e % .10e % .10e % .10e\n",
//                 pOutTilde[0][0],pOutTilde[0][1],pOutTilde[0][2],pOutTilde[0][3]);
//         printf("p4    : % .10e % .10e % .10e % .10e\n",
//                 pOut[1][0],pOut[1][1],pOut[1][2],pOut[1][3]);
//         printf("pT4   : % .10e % .10e % .10e % .10e\n",
//                 pOutTilde[1][0],pOutTilde[1][1],pOutTilde[1][2],pOutTilde[1][3]);
        
//         printf("Dips     =  % .16e\n",Dips[0]);
    }
    
    return;
}

// a check to see whether the plus PlusDistro causes instabilities
void checkEta(int low, int high) {
    int i,j,k, mu;
    double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}}; 
    double pOut[3][4], pInTilde[2][4], pOutTilde[2][4], vars[5],x[100][5], Dip[3][4];
    double jac, sprime, E12, xint, Mx, M1, eta, sja, sTja, mu2, muT2;
    
    x[0][0]=((double)rand())/((double) RAND_MAX);
    x[0][1]=((double)rand())/((double) RAND_MAX);
    for(i=low;i<high;i++) {
        eta = pow(0.1,i);
        xint = 1.0 - eta;
        printf("\n");
        printf("xint      =  % .16e \n", xint);
        PhSp_TaTa(Etot,x[0],vars,pOut,&jac);
        M1 = ME_Born_ee_tata(pIn, pOut);
        printf("M1        =  % .16e \n", M1);
        sprime = xint*shat;
        E12 = sqrt(sprime);
//         printf("E12      =  % .16e\n",E12);
        pInTilde[0][0] = E12/2;
        pInTilde[0][1] = 0.0;
        pInTilde[0][2] = 0.0;
        pInTilde[0][3] = E12/2;
        for(mu=0;mu<3;mu++) pInTilde[1][mu] = pInTilde[0][mu];
        pInTilde[1][3] = -pInTilde[0][3];
        PhSp_TaTa(E12,x[0],vars,pOutTilde,&jac);
        Mx = ME_Born_ee_tata(pInTilde, pOutTilde); // calculate the rescaled matrix element.
        printf("Mx        =  % .16e \n", Mx);
//         printf("p3    : % .10e % .10e % .10e % .10e\n",
//                 pOut[0][0],pOut[0][1],pOut[0][2],pOut[0][3]);
//         printf("pT3   : % .10e % .10e % .10e % .10e\n",
//                 pOutTilde[0][0],pOutTilde[0][1],pOutTilde[0][2],pOutTilde[0][3]);
//         printf("p4    : % .10e % .10e % .10e % .10e\n",
//                 pOut[1][0],pOut[1][1],pOut[1][2],pOut[1][3]);
//         printf("pT4   : % .10e % .10e % .10e % .10e\n",
//                 pOutTilde[1][0],pOutTilde[1][1],pOutTilde[1][2],pOutTilde[1][3]);
        sja = 2*metric(pIn[0],pOut[0]);
        sTja = 2*metric(pIn[0],pOutTilde[0]);
        
        mu2 = MT2/sja;
        muT2 = MT2/sTja;
        for(j=0;j<3;j++)
            for(k=0;k<4;k++) Dip[j][k] = 0.0;
//         getIntDipII(xint,Dip);
//         getIntDipIF(xint,sja,sTja,mu2,muT2,Dip);
        getIntDipFI(xint,sja,sTja,mu2,muT2,Dip);
//         printf("Dip[2]   =  % .16e\n",Dip[2]);
//         printf("Dip[3]   =  % .16e\n",Dip[3]);
//         printf("Mx-M1    =  % .16e\n",Mx-M1);
        double xPart = Mx*Dip[0][2], subPart = M1*Dip[0][3],
               plusPart = xPart - subPart,
               oCanc = -log10(abs((subPart-xPart)/xPart));
        printf("subPart   =  % .16e\n",subPart);
        printf("xPart     =  % .16e\n",xPart);
        printf("plusPart  =  % .16e\n",plusPart);
        printf("deltaPart =  % .16e\n",M1*Dip[0][0]);
        printf("regPart   =  % .16e\n",Mx*Dip[0][1]);
        printf("result    =  % .16e\n",plusPart + M1*Dip[0][0] + Mx*Dip[0][1]);
        printf("oCanc     =  % .2f\n",oCanc);
    }
    
    return;
}

// Check the chi implementation.
void checkChi(int N) {    
    int i,j;
    double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}},
          pOut[2][4], x[2], vars[2], jac;
    complex chi[2];
    
    printf("\nI = % .10e %+.10ei\n",real(I),imag(I));
          
    printf("\nCalculating chi:\n");
    
    
    x[0]=((double)rand())/((double) RAND_MAX);
    for(i=0;i<N;i++) {
        x[1]= 0.0 + pow(0.1,i+1);
        PhSp_TaTa(Etot,x,vars,pOut,&jac);
        printf("\nPoint %d\n",i+1);
        printf("\nx[1] = %+.2e\n",x[1]);
        
        printf("p3       :");
        for(j=0;j<4;j++) printf(" % .10f,",pOut[0][j]);
        printf("\n");
        
//         printf("p4       :");
//         for(j=0;j<4;j++) printf(" % .10f,",pOut[1][j]);
//         printf("\n");
        
        getChi(1,chi, pOut[0]);       
        printf("chi+(p3))    :\n");
        for(j=0;j<2;j++) printf("  % .10e %+.10ei\n",real(chi[j]),imag(chi[j]));
        printf("\n");

        getChi(-1,chi, pOut[0]);       
        printf("chi-(p3))    :\n");
        for(j=0;j<2;j++) printf("  % .10e %+.10ei\n",real(chi[j]),imag(chi[j]));
        printf("\n");
    }
    return;
}

// Check the spinor implementation.
void checkSpinors(int N) {    
    int i,j;
    double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}},
          pOut[2][4], x[2], vars[2], jac;
    complex u1[4],v2[4],vb2[4],v3[4],u4[4],ub4[4],
                         pSlash[4][4][4];
    
    printf("\nI = % .10e %+.10ei\n",real(I),imag(I));
          
    printf("\nCalculating spinors:\n");
    
    
    for(i=0;i<N;i++) {
        x[0]=((double)rand())/((double) RAND_MAX);
        x[1]=((double)rand())/((double) RAND_MAX);
        PhSp_TaTa(Etot,x,vars,pOut,&jac);
        getU(-1,u1,pIn[0]);
        getV(-1,v2,pIn[1]);
        sbar(vb2,v2);
        getV(-1,v3,pOut[0]);
        getU(-1,u4,pOut[1]);
        sbar(ub4,u4);
        for(j=0;j<2;j++) {
            slash(pSlash[j],pIn[j]);
            slash(pSlash[j+2],pOut[j]);
        }
        printf("\nPoint %d\n",i+1);

// boring since they are always the same...        
        printf("u(p1)    :\n");
        for(j=0;j<4;j++) printf("  % .10e %+.10ei\n",real(u1[j]),imag(u1[j]));
        printf("\n");
        
        printf("vbar(p2) :\n");
        for(j=0;j<4;j++) printf("  % .10e %+.10ei\n",real(vb2[j]),imag(vb2[j]));
        printf("\n");
        
        
        printf("p3       :");
        for(j=0;j<4;j++) printf(" % .10f,",pOut[0][j]);
        printf("\n");
        
        printf("v(p3)    :\n");
        for(j=0;j<4;j++) printf(" % .10f %+.10fi\n",real(v3[j]),imag(v3[j]));
        printf("\n");
        
        printf("ubar(p4) :\n");
        for(j=0;j<4;j++) printf(" % .10f %+.10fi\n",real(ub4[j]),imag(ub4[j]));
        printf("\n");
    }
    return;
}


// Check the currents implementation.
void checkCurrents(int N) {
    int i,j,k;
    double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}},
           pOut[2][4], x[2], vars[2], jac, M2sq, M2c;
           
    complex j1[4][4],j2[4][4];
    
//     printf("\nI = % .10f %+.10fi\n",real(I),imag(I));
          
    printf("\nCalculating Currents:\n");
    
    
    for(i=0;i<N;i++) {
        x[0]=((double)rand())/((double) RAND_MAX);
        x[1]=((double)rand())/((double) RAND_MAX);
        PhSp_TaTa(Etot,x,vars,pOut,&jac);
        printf("\nPoint %d\n",i+1);
        
        printf("p3       :");
        for(j=0;j<4;j++) printf(" % .10f,",pOut[0][j]);
        printf("\n");
        
        getCurr(pIn,pOut,j1,j2);
        

// initital state currents are simple
        for(j=0;j<4;j++) {
            for(k=0;k<4;k++) {
                printf("j1[%d][%d]  = % .10f %+.10fi\n",j,k,real(j1[j][k]),imag(j1[j][k]));
            }
            printf("\n");
        }
        
// final state currents
        
        for(j=0;j<4;j++) {
            for(k=0;k<4;k++) {
                printf("j2[%d][%d]  = % .10f %+.10fi\n",j,k,real(j2[j][k]),imag(j2[j][k]));
            }
            printf("\n");
        }
        
    }
   return;
}


// this calculates the tree level unsquared matrix element from 
void checkCurrentBornME(int N)
{
    int i,j,k;
    
    double pIn[2][4]= {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}},
           pOut[2][4], p1p2 = metric(pIn[0],pIn[1]), p1p3,
           M2, aM2, M2c, x[2], vars[2], jac, pref;
    
    complex j1[4][4],j2[4][4],bornMatrix[4][4],bornMatrix2[4][4],sum(0.0,0.0);
    

    for(i=0;i<N;i++) {
        sum = 0.0;
        x[0]=((double)rand())/((double) RAND_MAX);
        x[1]=((double)rand())/((double) RAND_MAX);
        PhSp_TaTa(Etot,x,vars,pOut,&jac);
        p1p3 = metric(pIn[0],pOut[0]);
        printf("\n");
        
        getCurr(pIn,pOut,j1,j2);
        
        currBornME(bornMatrix,p1p2,p1p3,pref,j1,j2);
        
//         printf("sum       = % .10f %+.10fi\n",real(sum),imag(sum));
        // sum up the matrix element
        for(j=0;j<4;j++) {
            for(k=0;k<4;k++) sum += bornMatrix[j][k]*conj(bornMatrix[j][k]);
        }
        // average
        M2 = real(sum)/4.0;
        
        currAntiBornME(bornMatrix2,p1p2,p1p3,pref,j1,j2);
        
        sum = 0.0;
        
        // sum up the matrix element
        for(j=0;j<4;j++) {
            for(k=0;k<4;k++) {
                sum += bornMatrix2[j][k]*conj(bornMatrix[j][k]);
            }
        }
        // average
        aM2 = real(sum)/4.0;
//         printf("M2       = % .10f %+.10fi\n",real(M2),imag(M2));
//         M2sq = real(M2*conj(M2))/4;
        printf("M2       = % .16e\n",M2);
        printf("aM2      = % .16e\n",aM2);
        M2c = ME_Born_ee_tata(pIn,pOut);
        printf("M2c      = % .16e\n",M2c);
        printf("M2/M2c   = % .10f\n",M2/M2c);
        printf("1-M2/M2c = % .10e\n",1.0-M2/M2c);
    }
    
    return;
}


// check the qcdloop integrals
void checkQCDLoop(int N) {
    std::vector<complex> res1(3),res2(3),res3(3),res4(4);
    std::vector<double> m(2),m1(3),m2(3),m3(3)/*{0.0,0.0,MT2,0.0}*/;
    std::vector<double> p(1),p1(3),p2(3),p3(3)/*{0.0,0.0,MT2,MT2,shat,0.0}*/;
    ql::QCDLoop<complex, double, double > qcdloop;
//     double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}}; 
//     double pOut[3][4], pInTilde[2][4], pOutTilde[2][4], vars[5],x[100][5], Dip[4];
    
    int i;
    double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}},
           pOut[2][4], x[2], vars[2], jac, p1p2 = shat/2.0,p1p3,p1p4,p2p3,p2p4;
        
//     B0(res1,shat,MT2,MT2);
//         printf("O(eps^0)  =  % .16e %+.16ei\n",real(res1[0]),imag(res1[0]));
//         printf("O(eps^-1) =  % .16e %+.16ei\n",real(res1[1]),imag(res1[1]));
//         printf("O(eps^-2) =  % .16e %+.16ei\n",real(res1[2]),imag(res1[2]));
//     B01(res4,shat,MT2,MT2);
//         printf("O(eps^0)  =  % .16e %+.16ei\n",real(res4[0]),imag(res4[0]));
//         printf("O(eps^-1) =  % .16e %+.16ei\n",real(res4[1]),imag(res4[1]));
//         printf("O(eps^-2) =  % .16e %+.16ei\n",real(res4[2]),imag(res4[2]));
//         printf("O(eps^-1) =  % .16e %+.16ei\n",real(res4[3]),imag(res4[3]));
//     B0(res1,0.0,ME2,ME2);
//         printf("O(eps^0)  =  % .16e %+.16ei\n",real(res1[0]),imag(res1[0]));
//         printf("O(eps^-1) =  % .16e %+.16ei\n",real(res1[1]),imag(res1[1]));
//         printf("O(eps^-2) =  % .16e %+.16ei\n",real(res1[2]),imag(res1[2]));
//     B02(res4,shat,MT2,MT2);
//         printf("O(eps^0)  =  % .16e %+.16ei\n",real(res4[0]),imag(res4[0]));
//         printf("O(eps^-1) =  % .16e %+.16ei\n",real(res4[1]),imag(res4[1]));
//         printf("O(eps^-2) =  % .16e %+.16ei\n",real(res4[2]),imag(res4[2]));
//         printf("O(eps^-1) =  % .16e %+.16ei\n",real(res4[3]),imag(res4[3]));
//     B0(res1,0.0,MMU2,MMU2);
//         printf("O(eps^0)  =  % .16e %+.16ei\n",real(res1[0]),imag(res1[0]));
//         printf("O(eps^-1) =  % .16e %+.16ei\n",real(res1[1]),imag(res1[1]));
//         printf("O(eps^-2) =  % .16e %+.16ei\n",real(res1[2]),imag(res1[2]));
//     m[0] = MT2; m[1] = MT2; p[0] = 0.0;
//     qcdloop.integral(res1,MT2,m,p);
//     printf("O(eps^0)  =  % .16e %+.16ei\n",real(res1[0]),imag(res1[0]));
//     printf("O(eps^-1) =  % .16e %+.16ei\n",real(res1[1]),imag(res1[1]));
//     printf("O(eps^-2) =  % .16e %+.16ei\n",real(res1[2]),imag(res1[2]));
// //     p[0] = shat;
// //     m[0] = MT2; m[1] = MT2;
// //     qcdloop.integral(resultc,1.0,m,p);
// //     printf("O(eps^0)  =  % .16e %+.16ei\n",real(resultc[0]),imag(resultc[0]));
// //     printf("O(eps^-1) =  % .16e %+.16ei\n",real(resultc[1]),imag(resultc[1]));
// //     printf("O(eps^-2) =  % .16e %+.16ei\n",real(resultc[2]),imag(resultc[2]));
// //     printf("p[5] = % .10e\n",p[5]);
//     return;
//     for(i=0;i<N;i++) {
//         x[0]=((double)rand())/((double) RAND_MAX);
//         x[1]=((double)rand())/((double) RAND_MAX);
//         PhSp_TaTa(Etot,x,vars,pOut,&jac);
//         
//         p1p3 = metric(pIn[0],pOut[0]);
//         p1p4 = p1p2 - p1p3; 
//         p2p3 = p1p4; p2p4 = p1p3;
//         
//         
//         printf("\nPoint %d\n",i+1);
//         
//         m1[0]=0.0;
//         m1[1]=0.0;
//         m1[2]=MT2;
//         
// //         m1[0]=0.0;
// //         m1[1]=0.0;
// //         m1[2]=MT2;
// //         m1[3]=0.0;
//         
//         p1[0]=0.0;
//         p1[1]=MT2;
//         p1[2]=MT2-2*p2p3;;
//         
// //         p1[0]=0.0;
// //         p1[1]=MT2;
// //         p1[2]=MT2;
// //         p1[3]=0.0;
// //         p1[4]=MT2-2*p2p4;
// //         p1[5]=2*p1p2;
//         
//         qcdloop.integral(res1,1.0,m1,p1);
//         
//         printf("O(eps^0)  =  % .16e %+.16ei\n",real(res1[0]),imag(res1[0]));
//         printf("O(eps^-1) =  % .16e %+.16ei\n",real(res1[1]),imag(res1[1]));
//         printf("O(eps^-2) =  % .16e %+.16ei\n",real(res1[2]),imag(res1[2]));
//         
//         
//         m2[0]=0.0;
//         m2[1]=MT2;
//         m2[2]=0.0;
//         
// //         m2[0]=0.0;
// //         m2[1]=0.0;
// //         m2[2]=MT2;
// //         m2[3]=0.0;
//         
//         p2[0]=MT2-2*p2p3;
//         p2[1]=MT2;
//         p2[2]=0.0;
//         
// //         p2[0]=0.0;
// //         p2[1]=MT2;
// //         p2[2]=MT2;
// //         p2[3]=0.0;
// //         p2[4]=MT2-2*p2p3;
// //         p2[5]=2*p1p2;
//         
//         qcdloop.integral(res2,1.0,m2,p2);
//         
//         printf("O(eps^0)  =  % .16e %+.16ei\n",real(res2[0]),imag(res2[0]));
//         printf("O(eps^-1) =  % .16e %+.16ei\n",real(res2[1]),imag(res2[1]));
//         printf("O(eps^-2) =  % .16e %+.16ei\n",real(res2[2]),imag(res2[2]));
//    }
//         m3[0]=0.0;
//         m3[1]=MT2;
//         m3[2]=0.0;
//         
//         p3[0]=MT2;
//         p3[1]=MT2;
//         p3[2]=shat;
//         qcdloop.integral(res3,1.0,m3,p3);
//         printf("O(eps^0)  =  % .16e %+.16ei\n",real(res3[0]),imag(res3[0]));
//         printf("O(eps^-1) =  % .16e %+.16ei\n",real(res3[1]),imag(res3[1]));
//         printf("O(eps^-2) =  % .16e %+.16ei\n",real(res3[2]),imag(res3[2]));
// this guy (C0(-p4/p3, -p1 - p2, 0, MT, 0) is indeed finite.
   return;
}


// check the fermion bubble
void checkBubble(int N) {
    std::vector<double> pA={},pB={shat};
    std::vector<complex> resA0,resA[9],resB0,resB[9];
    ql::QCDLoop<complex, double, double > qcdloop;
    
    int i;
    double pIn[2][4] = {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}},
           pOut[2][4], x[2], vars[2], jac, p1p2 = shat/2.0,p1p3,p1p4,p2p3,p2p4,
           musq;
    
    std::vector<double> mA0={0.0},mB0={0.0,0.0},
                        mfA[9] = {{ME2},{MMU*MMU},{MT2},
                                {MU*MU},{MC*MC},{MTO*MTO},
                                {MD*MD},{MS*MS},{MB*MB}},
                        mfB[9] = {{ME2,ME2},{MMU*MMU,MMU*MMU},{MT2,MT2},
                                {MU*MU,MU*MU},{MC*MC,MC*MC},{MTO*MTO,MTO*MTO},
                                {MD*MD,MD*MD},{MS*MS,MS*MS},{MB*MB,MB*MB}};
    
    musq = MZ2;
                                
    qcdloop.integral(resA0,musq,mA0,pA);
    qcdloop.integral(resB0,musq,mB0,pB);
    for(i=0;i<9;i++) {
        qcdloop.integral(resA[i],musq,mfA[i],pA);
        
        printf("mf%d^2  =  % .10f\n",i,mfB[i][0]);
        printf("Result for A0(mf%d):\n",i);
        printf("O(eps^0)  =  % .16e %+.16ei\n",real(resA[i][0]),imag(resA[i][0]));
        printf("O(eps^-1) =  % .16e %+.16ei\n",real(resA[i][1]),imag(resA[i][1]));
        printf("O(eps^-2) =  % .16e %+.16ei\n\n",real(resA[i][2]),imag(resA[i][2]));
        
        qcdloop.integral(resB[i],musq,mfB[i],pB);
        
        printf("Result for B0(s,mf%d,mf%d):\n",i,i);
        printf("O(eps^0)  =  % .16e %+.16ei\n",real(resB[i][0]),imag(resB[i][0]));
        printf("O(eps^-1) =  % .16e %+.16ei\n",real(resB[i][1]),imag(resB[i][1]));
        printf("O(eps^-2) =  % .16e %+.16ei\n\n",real(resB[i][2]),imag(resB[i][2]));
    }
          
    return;
}

// this calculates the tree level unsquared matrix element from 
void checkVirtME(int N)
{
    int i,j,k;
    
    double pIn[2][4]= {{Etot/2, 0.0, 0.0, Etot/2},{Etot/2, 0.0, 0.0, -Etot/2}},
           pOut[2][4], Mv, M2c, x[2], vars[2], jac, weight, result[2];
    
    complex res[4];
    
    initVirtCor();
    
//     printf("\nC0s values:\n");
//     printf("O(eps^0)  =  % .16e %+.16ei\n",real(C0s[0]),imag(C0s[0]));
//     printf("O(eps^-1) =  % .16e %+.16ei\n",real(C0s[1]),imag(C0s[1]));
//     printf("O(eps^-2) =  % .16e %+.16ei\n",real(C0s[2]),imag(C0s[2]));
    

    for(i=0;i<N;i++) {
        x[0]=((double)rand())/((double) RAND_MAX);
        x[1]=((double)rand())/((double) RAND_MAX);
        PhSp_TaTa(Etot,x,vars,pOut,&jac);
        printf("\nPoint %d\n",i+1);
        
        for(j=0;j<2;j++) {
            printf("p%d       :",j+1);
            for(k=0;k<4;k++) printf(" % .10f,",pIn[j][k]);
            printf("\n");
        }
        for(j=0;j<2;j++) {
            printf("p%d       :",j+3);
            for(k=0;k<4;k++) printf(" % .10f,",pOut[j][k]);
            printf("\n");
        }
        printf("\n");
        
        virtualMatrix(res,pIn,pOut);
        
        printf("finite       =  % .16e %+.16ei\n",real(res[0]),imag(res[0]));
        printf("IR O(eps^-1) =  % .16e %+.16ei\n",real(res[1]),imag(res[1]));
        printf("IR O(eps^-2) =  % .16e %+.16ei\n",real(res[2]),imag(res[2]));
        printf("UV O(eps^-1) =  % .16e %+.16ei\n",real(res[3]),imag(res[3]));
        
//         analyticPolesIntegrand(x, &weight, result);
//         
//         printf("IR O(eps^-1) =  % .16e\n",result[0]);
//         printf("IR O(eps^-2) =  % .16e\n",result[1]);
    }
    
    return;
}

// // this calculates the analytically extracted poles from the virt. corrections
// void checkVirtIR(int N)
// {
//     
//     
//     
//     return;
// }
