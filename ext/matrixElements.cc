extern "C" {
    // prototype for the wrapped helas e-e+->ta+ta-A process matrix element
    void m3sq(int* Npart, double mom[5][4], double* Result);
    void m2sqew(int* Npart, double mom[4][4], double* Result);
}

double ME_Born_ee_tata(double pIn[2][4], double pOut[2][4])
{
    double p1p2 = metric(pIn[0],pIn[1]),
           p1p3 = metric(pIn[0],pOut[0]),
           p1p4 = p1p2 - p1p3,//metric(pIn[0],pOut[1]), 
           p2p3 = p1p4,//metric(pIn[1],pOut[0]),
           p2p4 = p1p3;//metric(pIn[1],pOut[1]),
    
    return 2*e4*(p1p2*MT2 + p1p4*p2p3 + p1p3*p2p4)/(p1p2*p1p2);
}

// wrapped madgraph version 
// ATTENTION OBSOLETE since my version is correct!
// double ME_Born_ee_tata_check(double pIn[2][4], double pOut[2][4])
// {
//     int i,j,set = 4;
//     double Msq, mom[5][4];
//     
//     for(i=0;i<2;i++) {
//         for(j=0;j<4;j++) {
//             mom[i][j] = pIn[i][j];
//             mom[i+2][j] = pOut[i][j];
//         }
//     }
//     
//     m2sq(&set,mom,&Msq);
//     
//     return Msq;
// }

double ME_Born_ee_AZ_tata(double pIn[2][4], double pOut[2][4])
{
    double p1p2 = metric(pIn[0],pIn[1]),
           p1p3 = metric(pIn[0],pOut[0]),
           p1p4 = p1p2 - p1p3,//metric(pIn[0],pOut[1]), 
           p2p3 = p1p4,//metric(pIn[1],pOut[0]),
           p2p4 = p1p3,//metric(pIn[1],pOut[1]),
           gm2 = gm*gm, gp2 = gp*gp,
           propZ = 2*p1p2 - MZ2, // Z propagator
//            inv = p1p2*MT2 + p1p4*p2p3 + p1p3*p2p4, // same for all
           AAsq = 0.0, AZsq = 0.0, Aint = 0.0;
    
    // squared matrix element for QED process
    AAsq = 2*e4*(p1p2*MT2 + p1p4*p2p3 + p1p3*p2p4)/(p1p2*p1p2);
    
    // squared matrix element for weak process via a single Z
    AZsq = 4*e4*(p1p2*MT2*gm*gp*(gm2+gp2) + 2*p1p4*p2p3*gm2*gp2 + p1p3*p2p4*(gm2*gm2+gp2*gp2))
            /(propZ*propZ);
    
    // interference matrix element
    Aint = e4*(p1p2*MT2*(gm+gp)*(gm+gp) + 4*p1p4*p2p3*gm*gp + 2*p1p3*p2p4*(gm2+gp2))
            /(propZ*p1p2);
    
    
    return (AAsq + AZsq + 2*Aint);
}

// wrapped madgraph version 
// ATTENTION they dont exactly match which i think is related to different gm gp.
double ME_Born_ee_AZ_tata_check(double pIn[2][4], double pOut[2][4])
{
    int i,j,set = 4;
    double Msq, mom[5][4];
    
    for(i=0;i<2;i++) {
        for(j=0;j<4;j++) {
            mom[i][j] = pIn[i][j];
            mom[i+2][j] = pOut[i][j];
        }
    }
    
    m2sqew(&set,mom,&Msq);
    
    return Msq;
}

//////////////////////////////////////////
/*    real correction matrix elements   */
//////////////////////////////////////////


//My version of the real correction matrix element with set shat
double ME_ee_tataph(double pIn[2][4], double pOut[3][4])
{
    double sum = 0.0, p1p2 = shat/(2.0),// metric(pIn[0],pIn[1],
           p1p3 = metric(pIn[0],pOut[0]), p1p4 = metric(pIn[0],pOut[1]),
           p2p3 = metric(pIn[1],pOut[0]), p2p4 = metric(pIn[1],pOut[1]),
           p3p4 = metric(pOut[0],pOut[1]),
           p1k = metric(pIn[0],pOut[2]), p2k = metric(pIn[1],pOut[2]),
           p3k = metric(pOut[0],pOut[2]), p4k = metric(pOut[1],pOut[2]);
    
    // calculate the initial matrix element
    sum += 2*e6*(MT2*(p1k*p1k - 2*p1k*p1p2 + 2*p1p2*p1p2 - 2*p1p2*p2k + p2k*p2k)
           - p1k*p1p4*p2p3 + 2*p1p2*p1p4*p2p3 - p1p4*p2k*p2p3 + 2*p1k*p2p3*p2p4 + p1k*p1p4*p3k 
           - p1p2*p1p4*p3k - p1p2*p2p4*p3k + p2k*p2p4*p3k - p1p2*p2p3*p4k + p2k*p2p3*p4k
           + p1p3*(2*p1p4*p2k - p1k*p2p4 + 2*p1p2*p2p4 - p2k*p2p4 + p1k*p4k - p1p2*p4k))
           /(p1k*p2k*(p1k - p1p2 + p2k)*(p1k - p1p2 + p2k));
           
    // calculate the final matrix element.
    sum += -2*e6*(MT2*MT2*p1p2*(p3k*p3k + p4k*p4k)
           - p3k*p4k*(-2*p1p4*p2p4*p3k + p1p3*p2k*(p3k + p3p4) + p1k*p2p3*(p3k + p3p4)
           - 2*p1p3*p2p3*p4k + p1p4*p2k*(p3p4 + p4k) + p1k*p2p4*(p3p4 + p4k)
           + p1p4*p2p3*(p3k + 2*p3p4 + p4k) + p1p3*p2p4*(p3k + 2*p3p4 + p4k))
           + MT2*(p1p3*p2k*p3k*p3k + p1k*p2p3*p3k*p3k + p1p4*p2p3*p3k*p3k
           - 2*p1p2*p3k*p3p4*p4k + p1p4*p2k*p4k*p4k + p1p4*p2p3*p4k*p4k
           + p1k*p4k*(2*p2k*p3k + p2p4*p4k) + p1p3*p2p4*(p3k*p3k + p4k*p4k)))
           /(p1p2*p1p2*p3k*p3k*p4k*p4k);
           
    // calculate the interference matrix element.
    sum += -2*e6*(MT2*(-(p2k*(p1p4*(p2k - 2*p1p2) + p1k*(p1p3 + p1p4 + p2k - p2p3 - p1k))) 
            + p1k*(p1k - 2*p1p2 + p2k)*p2p4)*p3k + p3k*(-(p1k*p1k*p2p3*p2p4) 
            + p1p3*(p2k + 2*p2p4)*(p1p4*p2k - p1k*p2p4) + p1p4*p2k*(2*p1p4*p2p3 - p2p4*p3k) 
            + p1k*p1p4*(p2p4*(-2*p2p3 + p3k) + p2k*(p2p3 - p3p4)) + p1k*p2k*p2p4*p3p4) 
            + MT2*(p1p3*p2k*(p2k - 2*p1p2) + 2*p1k*p1p2*p2p3 - p1k*p1k*(p2k + p2p3) 
            + p1k*p2k*(p1p3 + p1p4 + p2k - p2p3 - p2p4))*p4k 
            + (-2*p1p3*p1p3*p2k*p2p4 + p1k*p1k*p2p3*p2p4 + p1p4*p2k*(p1p2 - p2p3)*p3k 
            + p1k*(p1p4*p2p3*(p2k + 2*p2p3) - p1p4*(2*p2k + p2p3)*p3k
            + (p1p2 - 2*p2k)*(p2p3 - p2p4)*p3k - p2k*p2p3*p3p4) + p1p3*(-(p1p4*p2k*(p2k + 2*p2p3))
            + p2k*(-p1p2 + p2p4)*p3k + p1k*p2p4*(2*p2p3 + p3k) 
            + p1k*p2k*(-p2p4 + 2*p3k + p3p4)))*p4k + p1p3*(-p1k + p2k)*p2p3*p4k*p4k)
            /(p1k*p1p2*p2k*(p1k - p1p2 + p2k)*p3k*p4k);
           
    return sum;//(initialM + interfM + finalM);
}


// Madgraph version for the real correction Matrix element 
// ATTENTION OBSOLETE since my version is correct!
double ME_ee_tataph_check(double pIn[2][4], double pOut[3][4])
{
    int i,j,NPart=5;
    double mom[5][4],MGRes;
     
    for (i=0;i<4;i++) {
        mom[0][i]=pIn[0][i];
        mom[1][i]=pIn[1][i];
        mom[2][i]=pOut[0][i];
        mom[3][i]=pOut[1][i];
        mom[4][i]=pOut[2][i];
    }        
    //  print output 
//     for(i=0;i<5;i++) {
//         printf("p%d  :",i+1);
//         for(j=0;j<4;j++) printf(" % .10e",mom[i][j]);
//         printf("\n");
//     }
//     printf("Sum : % .10e % .10e % .10e % .10e \n\n",
//             mom[2][0]+mom[3][0]+mom[4][0], mom[2][1]+mom[3][1]+mom[4][1],
//             mom[2][2]+mom[3][2]+mom[4][2], mom[2][3]+mom[3][3]+mom[4][3]);
        
    
    m3sq(&NPart,mom,&MGRes);
    
//     printf("MGRes    : % .10e \n\n",MGRes);
    
    return MGRes;
}
