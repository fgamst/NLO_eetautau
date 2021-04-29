#include "head.h"
// foreign code
#include "vegas/nvegas.c"
#include "PSgen/genps.cc"
// own code
#include "ext/io.cc"
#include "ext/lorentz.cc"
#include "ext/ps_gen.cc"
#include "ext/spinors.cc"
#include "ext/currents.cc"
#include "ext/matrixElements.cc"
#include "ext/virtual.cc"
#include "ext/binning.cc"
#include "dipoles/dipoles.cc"
#include "dipoles/intdipoles.cc"
#include "ext/integrands.cc"
#include "ext/checks.cc"

extern "C" {
    // prototype for the wrapped helas subroutine that sets the standard model couplings
    void helascoupsm(int* N);
}


double ApproximateIntegral(int N)
{
	int i;
	double sum = 0.0;
	
	for (i=0; i<NBINS; i++) sum += fbin[0][i]/(N);

	return sum;
}

void initialize() 
{
    int i, j, k;
    
        // call the helas coupling constants of the standard model.
    i = 0; // i = 0 is default configuration
    helascoupsm(&i);
    
    //  initializing the bins per iteration
    for (i=0;i<NIT;i++) {
        for (j=0;j<3;j++) {
            for(k=0;k<=(NBINS+1);k++) itBin[i][j][k] = 0.0;
//             dbin[i][j] = 0.0;
        }
    }
    
    absp3 = sqrt((Etot-2*MT)*(Etot+2*MT))/2;//+sqrt((Etot-2*MT)*(Etot+2*MT))/2/(NBINS);
    ymax = 0.5*log((Etot + absp3)/(Etot - absp3));
    
    // set the default renormalization sccale.
    muR = Etot;
    muR2 = muR*muR;
    muF2 = muR2;
    
    // set the cut parameter.
    acut_II = 0.01;
    acut_IF = 0.01;
    acut_FI = 0.01;
    acut_FF = 0.01;
//     acut_II = 0.1;
//     acut_IF = 0.1;
//     acut_FI = 0.1;
//     acut_FF = 0.1;
//     acut_II = 0.5;
//     acut_IF = 0.5;
//     acut_FI = 0.5;
//     acut_FF = 0.5;
//     acut_II = 1.0;
//     acut_IF = 1.0;
//     acut_FI = 1.0;
//     acut_FF = 1.0;
    
    //calculate some frequent logs.
    logMu = log(muR2);
    logMuMT = log(muR2/MT2);
    logMuMMU = log(muR2/MMU2);
    logMuME = log(muR2/ME2);
    logMuShat = log(muR2/shat);
    logShat = log(shat);
    logMeShat = log(ME2/shat);
    logMe = log(ME2);
}

int main(int argc, char **argv)
{
    int i,j = 1,k;
    double sum = 0.0, stdev;
    
//     printf("argc = %d\n",argc);
//     
//     for(i=0;i<argc;i++) printf("%s\n",argv[i]);
    
    
    std::string fname("");//("default_");
    
    initialize();
    
    if(argc>1) {
        if(argc>3) {
            fname.append(argv[1]);
            fname.append("_");
            fname.append (argv[2]);
            fname.append("_");
            fname.append(argv[3]);
//             fname.append = strcat(fname,"_");
//             fname.append = strcat(fname,argv[2]);
//             fname.append = strcat(fname,"_");
//             fname.append = strcat(fname,argv[3]);
            j = atoi(argv[3]);
            printf("j = %d\n",j);
            if(argv[2][0] == 'a') {
                if(argv[2][1] == 'I') {
                    if(argv[2][2] == 'I') {
                        acut_II = pow(0.1,j);
                        printf("alpha_II = %.5f\n", acut_II);
                    }
                    if(argv[2][2] == 'F') {
                        acut_IF = pow(0.1,j);
                        printf("alpha_IF = %.5f\n", acut_IF);
                    }
                }
                if(argv[2][1] == 'F') {
                    if(argv[2][2] == 'I') {
                        acut_FI = pow(0.1,j);
                        printf("alpha_FI = %.5f\n", acut_FI);
                    }
                    if(argv[2][2] == 'F') {
                        acut_FF = pow(0.1,j);
                        printf("alpha_FF = %.5f\n", acut_FF);
                    }
                }
            }
            if(argv[2][0] == 'm') {
                if(argv[2][1] == 'u') {
                    if(argv[4][0] == 'm') {
                        if(argv[4][1] == 'e') muR = j*(ME);
                        if(argv[4][1] == 'm') muR = j*(MMU);
                        if(argv[4][1] == 't') muR = j*(MT);
                    }
                    if(argv[4][0] == 'M') {
                        if(argv[4][1] == 'W') muR = j*(MW);
                        if(argv[4][1] == 'Z') muR = j*(MZ);
                    }
                    if(argv[4][0] == 'E') muR = j*(Etot);
                    
                    fname.append(argv[4]);
                    
                    muR2 = muR*muR;
                    muF2 = muR2;
                    printf("mu = %d%s = %.5f\n",j,argv[4], muR);
                    
                    //calculate some frequent logs.
                    logMu = log(muR2);
                    logMuMT = log(muR2/MT2);
                    logMuMMU = log(muR2/MMU2);
                    logMuME = log(muR2/ME2);
                    logMuShat = log(muR2/shat);
                }
            }
        }
        else {
            fname.append("default_");
            fname.append(argv[1]);
//             strcat(fname,argv[1]);
        }
        printf("\n");
        std::cout << "filename : " << fname << "\n";
//         printf("filename : %s\n",fname);
        printf("NSMP     = %d\n",NSMP);
        if(argv[1][0] == 'b') {
            born(NSMP, NIT); // calculate the leading order matrix element.
        }
        if(argv[1][0] == 'r') {
            realCor(NSMP,NIT); // integrate the real correction with vegas.
        }
        if(argv[1][0] == 'v') {
            virtCor(NSMP,NIT); // integrate the virtual matrix element with vegas.
        }
        if(argv[1][0] == 'i') {
            intDip(NSMP,NIT); // integrate the integrated dipoles with vegas.
        }
        if(argv[1][0] == 'a') {     
            analyticPoles(NSMP,NIT); // integrate the analytic IR-poles of the virtual diagrams.
        }
        // calculate the means of the bins and the standard deviation.
        for (i=0;i<3;i++) {
            for (j=0;j<=(NBINS);j++) {
                sum = 0.0;
                for(k=0;k<(NIT);k++) {
                    sum += itBin[k][i][j];
                }
                fbin[i][j] = sum/(NIT);
                sum = 0.0;
                for(k=0;k<(NIT);k++) {
                    sum += (itBin[k][i][j]-fbin[i][j])*(itBin[k][i][j]-fbin[i][j]);
                }
                dbin[i][j] = sqrt(sum/(NIT-1.0));
            }
        }
        OutputBins(fname); // output the bins to fname
    }
    else {
        printf("No valid integral specified\n");
    }
    
    
// Checks part
    
//     checkIntDipoleBin(10);
//     checkEta(1,5);
//    checkM2gZ(5);
//     checkM3(20);
//     checkM3Bin(20);
//     checkDip(20);
//     checkDipoleBin(10);
//     checkRealCorIntegrand(20);
//     checkPS3(100000);
//     checkChi(17);
//     checkSpinors(1);
//     checkCurrents(1);
//     checkCurrentBornME(20);
//     checkQCDLoop(10); 
//     checkBubble(10); 
//     checkVirtME(30);
    
            

  return(0);
}


