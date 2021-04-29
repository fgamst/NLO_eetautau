/*  spinor representation borrowed from: MURAYAMA et al.: HELAS, 1992; Appendix A.1              */

// helicity eigenspinor chi - tested
void getChi(char lam, std::complex<double> chi[2], double p[4])
{
	double abs3p = abs3(p);

// 	if((abs3p + p[3])/abs3p < 1e-14) {
	if(abs3p == - p[3]) {
		chi[(int)((1-lam)/2)] = 0.0;
		chi[(int)((1+lam)/2)] = lam;

	}
	else {
		chi[(int)((1-lam)/2)] = (abs3p + p[3])/sqrt(2*abs3p*(abs3p+p[3]));
		chi[(int)((1+lam)/2)] = (lam*p[1] + I*p[2])/sqrt(2*abs3p*(abs3p+p[3]));
	}

	return;
}


// u spinor - tested
void getU(char lam, std::complex<double> u[4],double p[4])
{
	double abs3p = abs3(p);
	std::complex<double> chi[2];
	getChi(lam,chi,p);

	u[0] = sqrt(p[0]-lam*abs3p)*chi[0];
	u[1] = sqrt(p[0]-lam*abs3p)*chi[1];
	u[2] = sqrt(p[0]+lam*abs3p)*chi[0];
	u[3] = sqrt(p[0]+lam*abs3p)*chi[1];

	return;
}


// v spinor - tested
void getV(char lam, std::complex<double> v[4],double p[4])
{
	double abs3p = abs3(p);
	std::complex<double> mchi[2];
	getChi(-lam,mchi,p);

	v[0] = -lam*sqrt(p[0]+lam*abs3p)*mchi[0];
	v[1] = -lam*sqrt(p[0]+lam*abs3p)*mchi[1];
	v[2] = lam*sqrt(p[0]-lam*abs3p)*mchi[0];
	v[3] = lam*sqrt(p[0]-lam*abs3p)*mchi[1];

	return;
}


// spinor bar operation - tested
void sbar(std::complex<double> sbar[4], std::complex<double> spinor[4])
{
	int i;

	for (i=0; i<4; i++) sbar[i] = conj(spinor[(i+2)%4]);

	return;
}


/* get the specified spinor type for momentum mom of helicity lam. - needs testing
// now also BARRED!!!
void getSpinor(char lam, char type, char bar, double complex s[4],double mom[4])
type = 0 for U spinor, 1 for V spinor
bar = 0 for non barred, 1 for barred spinor. 

{
	double complex chi[2];
	getChi(lam*type,chi,mom);
	double complex abs3p = abs3(mom);

//	the following is a bit ugly...
	s[0] = -lam*sqrt(mom[0]-lam*(1-2*type)*abs3p)*chi[0];
	s[1] = -lam*sqrt(mom[0]-lam*(1-2*type)*abs3p)*chi[1];
	s[2] = lam*sqrt(mom[0]+lam*(1-2*type)*abs3p)*chi[0];
	s[3] = lam*sqrt(mom[0]+lam*(1-2*type)*abs3p)*chi[1];

	return;
}
lets not do that mkay?*/


