#include "integrand.h"
#include "eventn.h"
#include <string.h>

double integrandsum(double a, double b, double c)
{
int i=0,k;
int j=0;
std::string quack;
double xxx[10000][13] = {0.0};
int evsize=0;
std::ifstream file;
std::string line;
file.open("../Daten/events/events.dat");
const char *formatStr = "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f";
char * pch;
if (file.is_open())
	{
	while(std::getline(file,line))
		{
		quack = line;
		char *cstr = new char[line.length() + 1];
		strcpy(cstr, line.c_str());
		
		  pch = strtok (cstr,",");
		  while (pch != NULL)
		  {
		    printf ("%s\n",pch);
			xxx[i][j] = strtod (pch, NULL);
		    pch = strtok (NULL, ",");
			j++;
		  }
		j=0;
		delete [] cstr;
		
		i++;
		evsize++;
		}
	file.close();
	}	
else
	{
	std::cerr << "Fehler";
	getchar();
	}
 const double mt=173.2;
 const double mt2=mt*mt;
 const double muF=mt;
 const double sqrtsg=13000.;
 const double sg=sqrtsg*sqrtsg;
 const double um=3.89379323e8;

 const double smin=mt2;
 const double smax=sg;
 const double srange=smax-smin;

 const double x2max=1.;

 const unsigned int n=2;
 double jacobi;
 double kout[n][4];
 const double masses[2]={mt,0.};

 double xx[4];
 double s;
 double flux;
 double x2min;
 double x2range;
 double x2;
 double sqrts;
 double x1;
 double sjacob;
 double p2[4];
 double p1[4];
 double x[2];
 double ergebnis;
 double ergebnissum=0.;

 std::cout << quack << "\n " << xxx[0][2] << " " << xxx[0][3] << " " << xxx[0][4] << std::endl;

 for(k=0; k<evsize; k++)
	{
	xx[0]=xxx[k][1];
	xx[1]=xxx[k][2];
	xx[2]=xxx[k][3];
	xx[3]=xxx[k][4];

 	s=smin+srange*xx[2];

 	flux=2.*s;

 	x2min=s/sg;
 	x2range=x2max-x2min;
 	x2=x2min+x2range*xx[3];

 	sqrts=sqrt(s);
 	x1=s/sg/x2;
 	sjacob=1./sg/x2;
 	p2[0]=sqrts/2.;
	p2[1]=0.;
	p2[2]=0.;
	p2[3]=sqrts/2.;
 	p1[0]=sqrts/2.;
	p1[1]=0.;
	p1[2]=0.;
	p1[3]=-sqrts/2.;
	x[0]=xx[0];
	x[1]=xx[1];
 	
 	eventn(n, sqrts, x, masses, kout, jacobi); 
 	ergebnis=um/flux*srange*x2range*sjacob*sig(p1, p2, kout[0], kout[1], x1, x2, a, b, c)*jacobi/pow(2.*M_PI,3*n-4);
	ergebnissum+=ergebnis;
	}
 return(ergebnissum);
}

