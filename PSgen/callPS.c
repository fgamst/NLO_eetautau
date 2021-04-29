#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "genps.h"


inline double lsp(double *Mom1, double *Mom2){
   return Mom1[0]*Mom2[0] - Mom1[1]*Mom2[1] - Mom1[2]*Mom2[2] - Mom1[3]*Mom2[3];
};


int main(int argc, char *argv[])
{

  int NPart=3, i;
  double x[5];
  double CMSEnergy,Jacobian;
  CMSEnergy = 11.7;
  
  double pOut[3][4];
  double Masses[3];
  Masses[0] = 1.7;
  Masses[1] = 1.7;
  Masses[2] = 0.0;

  double pIn[2][4];
  pIn[0][0] =+CMSEnergy/2.0;
  pIn[0][1] = 0.0;
  pIn[0][2] = 0.0;
  pIn[0][3] =+CMSEnergy/2.0;

  pIn[1][0] =+CMSEnergy/2.0;
  pIn[1][1] = 0.0;
  pIn[1][2] = 0.0;
  pIn[1][3] =-CMSEnergy/2.0;

  
  
  
  
  for( i = 1; i <= 10; i++){
  
//     //generate regular phase space
//     x[0]=(rand() % RAND_MAX) / (double)RAND_MAX;
//     x[1]=(rand() % RAND_MAX) / (double)RAND_MAX;
//     x[2]=(rand() % RAND_MAX) / (double)RAND_MAX;
//     x[3]=(rand() % RAND_MAX) / (double)RAND_MAX;
//     x[4]=(rand() % RAND_MAX) / (double)RAND_MAX;       
//     genps_(&NPart, &CMSEnergy, x, Masses, pOut, &Jacobian );

    
    
//  generate singular phase space:  
//     //  p5||p1
//     int col1=1 -1;
//     int col2=5 -1;
    
//     //  p5||p2   
//     int col1=2 -1;
//     int col2=5 -1;
    
    //  p5 soft
    int col1=2 -1;
    int col2=5 -1;
    
    double SingDepth = 1e-8;
    int Steps = 10;
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

    printf("p1.p5 : %e \n",lsp(pIn[0],pOut[2])/CMSEnergy/CMSEnergy);
    printf("p2.p5 : %e \n",lsp(pIn[1],pOut[2])/CMSEnergy/CMSEnergy);
    printf("p3.p5 : %e \n",lsp(pOut[0],pOut[2])/CMSEnergy/CMSEnergy);
    printf("p4.p5 : %e \n",lsp(pOut[1],pOut[2])/CMSEnergy/CMSEnergy);
    printf("E5/E  : %e \n",pOut[2][0]/CMSEnergy);
    
    getchar();
    
    printf("\n\n\n");
  };

  
  
  
  return(0);
}

