void OutputBins(std::string name)
{
    int i,j,k;
    std::string fname("out/");
//     std::string fname(name);
    std::ofstream output;
    
    fname.append(name); 
//     std::cout << "filename : " << fname << "\n";
//     fname.append("_realCor_");
//     fname.append(std::to_string(a));
    
//     printf("j = %d \n",j);
//     std::cout << fname << "\n";
    
    output.open (fname);
//     output.open ("bla");
    
    output.precision(16);
    
 	for (i=0; i<=NBINS; i++) {
        output.precision(3);
        output << 2*((double)i)/NBINS-1 << ",";// cos(theta)
        output.precision(16);
        output << fbin[0][i] <<","<< dbin[0][i] <<","; // cos(theta) bin
        // dump pT bins
        output << ((double)i)*absp3/NBINS << "," << fbin[1][i] << "," << dbin[1][i] << ",";
        // dump rapidity bins
//         output << ((double)i)*ymax/NBINS << "," << fbin[2][i] << "," << dbin[2][i] << "\n";
        output << 0.0 << "," << 0.0 << "," << 0.0 << "\n";
        for(j=0;j<3;j++) {
            for(k=0;k<NIT;k++) itBin[k][j][i] = 0.0;
//             dbin[j][i] = 0.0;
        }
    }
    output << fbin[0][NBINS+1] << "," <<fbin[1][NBINS+1] << "," << fbin[2][NBINS+1];

    output.close();
    
	return;
    
}

// void OutputBinsC(char *fname, int a)
// {
// 	int i,j;
//     /* = "aFI_int_"*/;
// //     fname.append(std::to_string(a));
//     FILE *output;
//     
// 	output = fopen("histo/histo.dat", "w");
// //     output.open (fname);
// // 
// //  	for (i=0; i<=NBINS; i++) {
// //         output << 2*((double)i)/NBINS-1 << "," << fbin[0][i] <<",";
// //         output << ((double)i)*absp3/NBINS << "," << fbin[1][i] << "\n";
// //         for(j=0;j<3;j++) fbin[j][i] = 0.0;
// //     }
// //     output << fbin[0][NBINS+1] << "," <<fbin[1][NBINS+1] << "," << fbin[2][NBINS+1];
// // 
// //     output.close();
//     
// 	return;
//     
// }

//eventually i want to do this with C again so i can specify the number of digits.

// Solution 1 from stackoverflow:
// char fname[128];
// printf("Enter .txt file name\n");
// scanf("%123s",fname);
// strcat(fname,".txt");
// FILE *inputf;
// inputf=fopen(fname,"w");

// Solution 2 from stackoverflow:
// #include <stdio.h>
// 
// void read_name(char *);
// 
// int main(void)
//                {
//                  char name[BUFSIZ];
//                  char line[BUFSIZ];
//                  FILE *f;
//                  printf("Name ");
//                  read_name(name);
//                  if ( (f=fopen(name,"r"))==NULL)
//                  return -1;
//                  else
//                  return 0;
//                  fclose(f);
//                }     
// void read_name(char *s)
//                {
//                  int i;
//                  fgets(s,BUFSIZ,stdin);
//                  for (i=0; s[i]!='\n'; i++);
//                  s[i]='\0';
//                  return;
//                 }
// 


// void OutputPS3old()
// {
// 	int i;
// 	FILE *output;
//     
//     // set xbins for PS3
// //     for (i=0;i<=NBINS;i++) {
// //         xbin[0][i] = (shat-4*MT2)*i/NBINS + 4*MT2;
// //         xbin[1][i] = 2*M_PI*i/NBINS;
// //         xbin[2][i] = 2*i/NBINS-1;
// //         xbin[3][i] = 2*M_PI*i/NBINS;
// //         xbin[4][i] = 2*i/NBINS-1;
// //     }
// 	
// 	output = fopen("sampling.csv", "w");
//  	for (i=0; i<=NBINS; i++) {
//         fprintf(output, "% .10e % .10e ",
//                 (shat-4*MT2)*((double)i)/NBINS + 4*MT2, fbin[0][i]);
//         fprintf(output, "% .10e % .10e ", 2*M_PI*((double)i)/NBINS, fbin[1][i]);
//         fprintf(output, "% .10e % .10e ", 2*((double)i)/NBINS-1, fbin[2][i]);
//         fprintf(output, "% .10e % .10e ", 2*M_PI*((double)i)/NBINS, fbin[3][i]);
//         fprintf(output, "% .10e % .10e\n", 2*((double)i)/NBINS-1, fbin[4][i]);
//     }
// 	fclose(output);
// 
// 	return;
// }
/*
void OutputPS3()
{
	int i;
	FILE *output;
    
    // set xbins for PS3
//     for (i=0;i<=NBINS;i++) {
//         xbin[0][i] = (shat-4*MT2)*i/NBINS + 4*MT2;
//         xbin[1][i] = 2*M_PI*i/NBINS;
//         xbin[2][i] = 2*i/NBINS-1;
//         xbin[3][i] = 2*M_PI*i/NBINS;
//         xbin[4][i] = 2*i/NBINS-1;
//     }
	
	output = fopen("histo/histo.dat", "w");
 	for (i=0; i<=NBINS; i++) {
        fprintf(output, "% .10e % .10e \n", 2*((double)i)/NBINS-1, fbin[i]);
    }
	fclose(output);

	return;
}*/
