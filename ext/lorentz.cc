// Vector Algebra in Minkowski space
// Minkowski metric - tested
double metric(double vec1[4], double vec2[4])
{
	return vec1[0]*vec2[0] - vec1[1]*vec2[1] - vec1[2]*vec2[2] - vec1[3]*vec2[3];
}


// Add Lorentz vectors - tested
void addLV(double sum[4], double vec1[4], double vec2[4])
{
	int i;

	for(i=0;i<4;i++) sum[i] = vec1[i] + vec2[i];
	
	return ;
}


// multiply LV with complex number - tested
void prefV(double result[4], double pref, double vec[4])
{
	int i;

	for(i=0;i<4;i++) result[i] = pref * vec[i];
	
	return;
}

// "slash" a lorentz vector.
void slash(std::complex<double> pSlash[4][4], double p[4])
{
    std::complex<double> 
        gamma0[4][4] = {{0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}, {0, 1, 0, 0}},
        gamma1[4][4] = {{0, 0, 0, 1}, {0, 0, 1, 0}, {0, -1, 0, 0}, {-1, 0, 0, 0}},
        gamma2[4][4] = {{0, 0, 0, -I}, {0, 0, I, 0}, {0, I, 0, 0}, {-I, 0, 0, 0}},
        gamma3[4][4] = {{0, 0, 1, 0}, {0, 0, 0, -1}, {-1, 0, 0, 0}, {0, 1, 0, 0}};
        
    int i,j;
        
    for(i=0;i<4;i++) {
        for(j=0;j<4;j++) 
            pSlash[i][j] = gamma0[i][j]*p[0] - gamma1[i][j]*p[1] 
                         - gamma2[i][j]*p[2] - gamma3[i][j]*p[3];
    }
    
    return;
}

// Absolute value of 3-momentum - tested
double abs3(double vec[4])
{
	return sqrt(vec[1]*vec[1] + vec[2]*vec[2] + vec[3]*vec[3]);
}

// Scalar product in euclidian space
double dot3(double vec1[4], double vec2[4])
{
	return (vec1[1]*vec2[1] + vec1[2]*vec2[2] + vec1[3]*vec2[3]);
}

// Calculate cos(theta) between two vectors in the same frame.
double cosTheta(double pIn[4], double pOut[4]) 
{
    return dot3(pIn,pOut)/abs3(pIn)/abs3(pOut);
}

// Lorentz-Boost along the designated axis.
void boost(double vecIn[4], double vecOut[4], int dir, double beta, double gamma)
{
	int i;
    
    for(i=0;i<4;i++)
        vecOut[i] = vecIn[i];
    
    vecOut[0] = gamma*(vecIn[0] - beta*vecIn[dir]);
    vecOut[dir] = gamma*(vecIn[dir] - beta*vecIn[0]);
    
    return;
}

// Calculates the Kallen function lam(x,y,z) = x^2 + y^2 + z^2 - 2xy - 2xz - 2yz
double kallen(double x, double y, double z)
{
	return x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z;
}

// Depending on what environment is used calculates the dilogarithm
double Li2(double x)
{
//     return 0.0;
	return gsl_sf_dilog(x);
}
