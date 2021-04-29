/*  polarization vector representation borrowed from: MURAYAMA et al.: HELAS, 1992; Appendix A.1              */

/* massive bosons not included. 10.57GeV too low to go on-shell. */
/*getEps(char lam, double k[4], double complex eps[4])
{
    double kT = sqrt(k[1]*k[1]+k[2]*k[2]), abs3k = sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);
    
    if (lam == 0)
    {
        eps[0] = abs3k/mk;
        eps[1] = k[1]*k[0]/mk/abs3k;
        eps[2] = k[2]*k[0]/mk/abs3k;
        eps[3] = k[3]*k[0]/mk/abs3k;
    }
    else 
    {
        eps[0] = 0;
        eps[1] = (-lam*k[1]*k[3]/abs3k + I*k[2])/kT;
        eps[2] = (-lam*k[2]*k[3]/abs3k - I*k[1])/kT;
        eps[3] = (lam*kT*kT/abs3k)/kT;
    }
    return;
}*/


/* polarization vector for massless bosons. */
void getEps(char lam, double k[4], double complex eps[4])
{
    double kT = sqrt(k[1]*k[1]+k[2]*k[2]), abs3k = sqrt(k[1]*k[1] + k[2]*k[2] + k[3]*k[3]);
    
    eps[0] = 0;
    eps[1] = (-lam*k[1]*k[3]/abs3k + I*k[2])/kT;
    eps[2] = (-lam*k[2]*k[3]/abs3k - I*k[1])/kT;
    eps[3] = (lam*kT*kT/abs3k)/kT;
    return;
}
