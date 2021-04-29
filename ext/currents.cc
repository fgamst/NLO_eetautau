// calculate the "standard" currents - tested
// I have some ideas how to bring this down to a few less lines but this is low prio
// double complex jnew(unsigned char spinor, unsigned char slash, char lam, char lam1, char lam2,
//                  double p[4][4])
// /*using the following convention for the lambdas:
// first string				second string
// vbar(lam2,p2).?.w(lam).u(lam1,p1)	ubar(lam2,p4).?.w(lam).v(lam1,p3)*/
// {	
//     int i;
// //     double sum = 0.0;
// 	if (spinor == 1) { // spinor string one i.e. vbar(p2).?.u(p1)
// 		double complex u1m[4],u1p[4],u1[4],v2m[4],v2p[4],vb2m[4],vb2p[4],vb2[4];
// 
// // 		get the spinors with correct helicity
// 		getU(-1,u1m,p[0]);
// 		getU(1,u1p,p[0]);
// 		getV(-1,v2m,p[1]);
// 		getV(1,v2p,p[1]);
// 		sbar(vb2m,v2m);
// 		sbar(vb2p,v2p);
// 
//         for(i=0;i<4;i++) {
//             u1[i] = u1m[i]+u1p[i];
//             vb2[i] = vb2m[i]+vb2p[i];
//         }
//         
//         
// 		if (slash == 0) { // j1(0,lambda)=vbar(p2).omega(lam).u(p1) 
// 			if (lam == -1) {
// 				return vb2[0]*u1[0]+vb2[1]*u1[1];
// 			}
// 			else {
// 				return vb2[2]*u1[2]+vb2[3]*u1[3];
// 			}
// 		}
// 		else{ // j1(1,lambda)=vbar(p2).slash(p3).omega(lambda).u(p1)
// 			if (lam == -1) {
// 				return
// 				vb2[2]*((p[2][0]+p[2][3])*u1[0]+(p[2][1]-I*p[2][2])*u1[1])+
// 				vb2[3]*((p[2][1]+I*p[2][2])*u1[0]+(p[2][0]-p[2][3])*u1[1]);
// 				
// 			}
// 			else {
// 				return
// 				vb2[0]*((p[2][0]-p[2][3])*u1[2]+(-p[2][1]+I*p[2][2])*u1[3])+
// 				vb2[1]*((-p[2][1]-I*p[2][2])*u1[2]+(p[2][0]+p[2][3])*u1[3]);
// 			}
// 		}
// 	}
// 	else { // spinor string two i.e. ubar(p4).?.v(p3)
// 		double complex u4m[4],u4p[4],ub4m[4],ub4p[4],ub4[4],v3[4];
// 
//  //		get the spinors with correct helicity
// 		getV(lam1,v3,p[2]);
// 		getU(lam2,u4,p[3]);
// 		sbar(ub4,u4);
// 
// 		if (slash == 0) { // j2(0,lambda)=(ubar(p4).omega(lambda).v(p3)
// 			if (lam == -1) {
// 				return ub4[0]*v3[0]+ub4[1]*v3[1];
// 			}
// 			else {
// 				return ub4[2]*v3[2]+ub4[3]*v3[3];
// 			}
// 		}
// 		else { // j2(1,lambda)=ubar(p4).slash(p1).omega(lambda).v(p3)
// 			if (lam == -1) {
// 				return
// 				ub4[2]*((p[0][0]+p[0][3])*v3[0]+(p[0][1]-I*p[0][2])*v3[1])+
// 				ub4[3]*((p[0][1]+I*p[0][2])*v3[0]+(p[0][0]-p[0][3])*v3[1]);
// 			}
// 			else {
// 				return
// 				ub4[0]*((p[0][0]-p[0][3])*v3[2]+(-p[0][1]+I*p[0][2])*v3[3])+
// 				ub4[1]*((-p[0][1]-I*p[0][2])*v3[2]+(p[0][0]+p[0][3])*v3[3]);
// 
// 			}
// 		}
// 	}
// 	return 0;
// }

std::complex<double> curr(unsigned char spinor, unsigned char slash,
                          char lam, char lam1, char lam2, double p[4][4])
/*using the following convention for the lambdas:
first string				second string
vbar(lam2,p2).?.w(lam).u(lam1,p1)	ubar(lam2,p4).?.w(lam).v(lam1,p3)*/
{	
	if (spinor == 1) { // spinor string one i.e. vbar(p2).?.u(p1)
		std::complex<double> u1[4],v2[4],vb2[4];

// 		get the spinors with correct helicity
		getU(lam1,u1,p[0]);
		getV(lam2,v2,p[1]);
		sbar(vb2,v2);     

		if (slash == 0) { // j1(0,lambda)=vbar(p2).omega(lam).u(p1) 
			if (lam == -1) {
				return vb2[0]*u1[0]+vb2[1]*u1[1];
			}
			else {
				return vb2[2]*u1[2]+vb2[3]*u1[3];
			}
		}
		else{ // j1(1,lambda)=vbar(p2).slash(p3).omega(lambda).u(p1)
			if (lam == -1) {
				return
				vb2[2]*((p[2][0]+p[2][3])*u1[0]+(p[2][1]-I*p[2][2])*u1[1])+
				vb2[3]*((p[2][1]+I*p[2][2])*u1[0]+(p[2][0]-p[2][3])*u1[1]);
				
			}
			else {
				return
				vb2[0]*((p[2][0]-p[2][3])*u1[2]+(-p[2][1]+I*p[2][2])*u1[3])+
				vb2[1]*((-p[2][1]-I*p[2][2])*u1[2]+(p[2][0]+p[2][3])*u1[3]);
			}
		}
	}
	else { // spinor string two i.e. ubar(p4).?.v(p3)
		std::complex<double> u4[4],ub4[4],v3[4];

 //		get the spinors with correct helicity
		getV(lam1,v3,p[2]);
		getU(lam2,u4,p[3]);
		sbar(ub4,u4);

		if (slash == 0) { // j2(0,lambda)=(ubar(p4).omega(lambda).v(p3)
			if (lam == -1) {
				return ub4[0]*v3[0]+ub4[1]*v3[1];
			}
			else {
				return ub4[2]*v3[2]+ub4[3]*v3[3];
			}
		}
		else { // j2(1,lambda)=ubar(p4).slash(p1).omega(lambda).v(p3)
			if (lam == -1) {
				return
				ub4[2]*((p[0][0]+p[0][3])*v3[0]+(p[0][1]-I*p[0][2])*v3[1])+
				ub4[3]*((p[0][1]+I*p[0][2])*v3[0]+(p[0][0]-p[0][3])*v3[1]);
			}
			else {
				return
				ub4[0]*((p[0][0]-p[0][3])*v3[2]+(-p[0][1]+I*p[0][2])*v3[3])+
				ub4[1]*((-p[0][1]-I*p[0][2])*v3[2]+(p[0][0]+p[0][3])*v3[3]);

			}
		}
	}
	return 0;
}

// calculates the currents.
void getCurr(double pIn[2][4], double pOut[2][4], complex j1[4][4], complex j2[4][4])
{
    int i,j;
    double p[4][4];
           
    for(i=0;i<2;i++) {
        for(j=0;j<4;j++) {
            p[i][j] = pIn[i][j];
            p[i+2][j] = pOut[i][j];
        }
    }
//         curr(unsigned char spinor, unsigned char slash,
//                           char lam, char lam1, char lam2, double p[4][4])
    
    //j1(0,-1) has only 1 nonvanishing element.
    j1[0][0] = curr(1,0,-1,-1,-1,p);
    j1[0][1] = /*curr(1,0,-1,1,-1,p) = */0.0;
    j1[0][2] = /*curr(1,0,-1,-1,1,p) = */0.0;
    j1[0][3] = /*curr(1,0,-1,1,1,p) = */0.0;
                        
    //j1(0,1) also only has 1 nonvanishing element.
    j1[1][0] = /*curr(1,0,1,-1,-1,p) = */0.0;
    j1[1][1] = /*curr(1,0,1,1,-1,p) = */0.0;
    j1[1][2] = /*curr(1,0,1,-1,1,p) = */0.0;
    j1[1][3] = /*curr(1,0,1,1,1,p) = */-j1[0][0];
    
    //j1(1,-1) has only 1 nonvanishing element.
    j1[2][0] = /*curr(1,1,-1,-1,-1,p) = */0.0;
    j1[2][1] = /*curr(1,1,-1,1,-1,p) = */0.0;
    j1[2][2] = curr(1,1,-1,-1,1,p);
    j1[2][3] = /*curr(1,1,-1,1,1,p) = */0.0;
                        
    //j1(1,1) also only has 1 nonvanishing element.
    j1[3][0] = /*curr(1,1,1,-1,-1,p) = */0.0;
    j1[3][1] = /*curr(1,1,1,1,-1,p) = */conj(j1[2][2]);
    j1[3][2] = /*curr(1,1,1,-1,1,p) = */0.0;
    j1[3][3] = /*curr(1,1,1,1,1,p) = */0.0;
    
    //j2(0,-1) the alternating spin parts vanish.
    j2[0][0] = curr(2,0,-1,-1,-1,p);
    j2[0][1] = /*curr(2,0,-1,1,-1,p) = */0.0;
    j2[0][2] = /*curr(2,0,-1,-1,1,p) = */0.0;
    j2[0][3] = curr(2,0,-1,1,1,p);
    
    //j2(0,1) the alternating spin parts vanish aswell.
    j2[1][0] = /*curr(2,0,1,-1,-1,p) = */ -conj(j2[0][3]);
    j2[1][1] = /*curr(2,0,1,1,-1,p) = */0.0;
    j2[1][2] = /*curr(2,0,1,-1,1,p) = */0.0;
    j2[1][3] = /*curr(2,0,1,1,1,p) = */ -conj(j2[0][0]);
    
    //j2(1,-1)
    j2[2][0] = curr(2,1,-1,-1,-1,p);
    j2[2][1] = curr(2,1,-1,1,-1,p);
    j2[2][2] = curr(2,1,-1,-1,1,p);
    j2[2][3] = curr(2,1,-1,1,1,p);
    
    //j2(1,1) all values coincide with some version of j2(1,-1)
    j2[3][0] = /*curr(2,1,1,-1,-1,p) = */-conj(j2[2][3]);
    j2[3][1] = /*curr(2,1,1,1,-1,p) = */j2[2][2];
    j2[3][2] = /*curr(2,1,1,-1,1,p) = */j2[2][1];
    j2[3][3] = /*curr(2,1,1,1,1,p) = */-conj(j2[2][0]);
    
    return;
}




