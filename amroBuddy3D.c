
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#include <stdbool.h>
typedef enum { false, true } bool;

#define DIM 100
#define ACC 10e-7
#define DIMXY 7
#define DIMZ 50

double t   =  1.;
double tp  = -0.12;
double tpp =  0.06;
double tz  =  0.07;
double mu  = -0.85;
double tau =  25;

// double B = 0.1;
// double theta = 0;
// double phi = 0;

double Bx = 0.01;//B*sin(theta)*cos(phi);
double By = 0.02;//B*sin(theta)*sin(phi);
double Bz = 0.05;//B*cos(theta);
double h  = 0.1;// 
// int    Nt  =  500;


bool wrapBrillouin=true;

double calculateDispersion(double kx, double ky, double kz){
	double  dispersion=0; 
	if (tz!=0){
		dispersion = -2.*tz*cos(kz/2.);
		dispersion *= (cos(kx)-cos(ky))*(cos(kx)-cos(ky));
		dispersion *=  cos(kx/2.)*cos(ky/2.);
	}
	dispersion += -mu;
	dispersion += -2.*t   * (cos(kx) + cos(ky));
	dispersion += -4.*tp  *  cos(kx)*cos(ky);
	dispersion += -2.*tpp * (cos(2.*kx) + cos(2.*ky));
	return  dispersion;
}

double *calculateVelocity(double kx, double ky, double kz){
	static double velocity[3];
	velocity[0]=velocity[1]=velocity[2] = 0;
	if (tz!=0){
		velocity[0]  = -2.*(cos(kx)-cos(ky))*sin(kx)*cos(kx/2.)*cos(ky/2.);
		velocity[0] += -(cos(kx)-cos(ky))*(cos(kx)-cos(ky))*sin(kx/2.)*cos(ky/2.) /2.;
		velocity[0] *=  -2.*tz*cos(kz/2.);
		
		velocity[1]  =  2*(cos(kx)-cos(ky))*sin(ky)*cos(kx/2.)*cos(ky/2.);
		velocity[1] += -(cos(kx)-cos(ky))*(cos(kx)-cos(ky))*cos(kx/2.)*sin(ky/2.) /2.;
		velocity[1] *= -2.*tz*cos(kz/2.);
		
		velocity[2]  =  tz*sin(kz/2.);
		velocity[2] *=  (cos(kx)-cos(ky))*(cos(kx)-cos(ky));
		velocity[2] *=  cos(kx/2.)*cos(ky/2.);
	}
	velocity[0] +=  2.*t    *  sin(kx);
	velocity[0] +=  4.*tp   *  sin(kx)*cos(ky);
	velocity[0] +=  4.*tpp  *  sin(2.*kx);
	velocity[1] +=  2.*t    *  sin(ky);
	velocity[1] +=  4.*tp   *  cos(kx)*sin(ky);
	velocity[1] +=  4.*tpp  *  sin(2.*ky);
	return  velocity;
}

double* calculateForce(double kx, double ky, double kz){
	static double force[3];
	double *velocity; 
	velocity = calculateVelocity(kx,ky,kz);
	force[0]  =  velocity[1]*Bz;
	force[0] += -velocity[2]*By;
	force[1]  = -velocity[0]*Bz;
	force[1] +=  velocity[2]*Bx;
	force[2]  =  velocity[0]*By;
	force[2] += -velocity[1]*Bx;
	return  force;
}


void calculateFS(double fermiSurface[DIMXY][DIMZ][3]) {
	for (int zz = 0; zz < DIMZ; zz++) {
		double kFz = -2*M_PI + (zz/((double)DIMZ-1))*4*M_PI;
		
		// theta is clockwise, from (pi,pi), with theta =0 being the 45deg axis (0,0)-(pi,pi)
		// the max value of angle depend on the position of the FS for kx = 0.
		double maxTheta = M_PI/4.; 
		double maxKy = M_PI;
		do {
			maxTheta -= 0.01;
			maxKy = M_PI*(1+tan(maxTheta-M_PI/4.));
		} while (calculateDispersion(0, maxKy, kFz) > 0.);
		
		for (int tt = 0; tt < DIMXY; tt++) {			
			double theta = (tt/((double)DIMXY-1)) * maxTheta;
			double kx1 = 0   , ky1 = M_PI*(1+tan(theta-M_PI/4.));
			double kx2 = M_PI, ky2 = M_PI;
			double xi1 = calculateDispersion(kx1, ky1, kFz);
			double xi2 = calculateDispersion(kx2, ky2, kFz);

			double xi, kFx = NAN, kFy = NAN;		

			printf("%f\n", xi1*xi2);fflush(stdout);	
			if (xi1*xi2 < 0) {
				do {
					kFx = (kx2+kx1)/2.;
					kFy = (ky2+ky1)/2.;
					xi = calculateDispersion(kFx, kFy, kFz);
					if (xi > 0) {
						kx2 = kFx;
						ky2 = kFy;
					} else if (xi < 0) {
						kx1 = kFx;
						ky1 = kFy;
					}
				} while (fabs(xi) > ACC);
			}
			fermiSurface[tt][zz][0] = kFx;
			fermiSurface[tt][zz][1] = kFy;
			fermiSurface[tt][zz][2] = kFz;
		}
	}
}

void printFS(double fermiSurface[DIMXY][DIMZ][3]){
	FILE *fileOut = fopen("FS.dat","w");
	fprintf(fileOut, "  kx         ky         kz\n");
	for (int zz = 0; zz < DIMZ; zz++) {
		double* kvec;
		for (int tt = DIMXY-1; tt >= 0; tt--) { 
			kvec = &fermiSurface[tt][zz][0]; 
			fprintf(fileOut,"  %2.6f   %2.6f   %2.6f  \n",  kvec[1],  kvec[0], kvec[2]);
		} for (int tt = 0; tt < DIMXY; tt++) { 
			kvec = &fermiSurface[tt][zz][0];
			fprintf(fileOut,"  %2.6f   %2.6f   %2.6f  \n",  kvec[0],  kvec[1], kvec[2]);
			if ( (fabs(fabs(kvec[1])-M_PI) < 0.1) || (fabs(fabs(kvec[2])-M_PI) < 0.01) ) fprintf(fileOut,"\n");
		} for (int tt = DIMXY-1; tt >= 0; tt--){
			kvec = &fermiSurface[tt][zz][0];
			fprintf(fileOut,"  %2.6f   %2.6f   %2.6f  \n", -kvec[0],  kvec[1], kvec[2]);
		} for (int tt = 0; tt < DIMXY; tt++) { 
			kvec = &fermiSurface[tt][zz][0]; 
			fprintf(fileOut,"  %2.6f   %2.6f   %2.6f  \n", -kvec[1],  kvec[0], kvec[2]);
			if ( (fabs(fabs(kvec[1])-M_PI) < 0.1) || (fabs(fabs(kvec[2])-M_PI) < 0.01) ) fprintf(fileOut,"\n");
		} for (int tt = DIMXY-1; tt >= 0; tt--) { 
			kvec = &fermiSurface[tt][zz][0];
			fprintf(fileOut,"  %2.6f   %2.6f   %2.6f  \n", -kvec[1], -kvec[0], kvec[2]);
		} for (int tt = 0; tt < DIMXY; tt++) { 
			kvec = &fermiSurface[tt][zz][0];
			fprintf(fileOut,"  %2.6f   %2.6f   %2.6f  \n", -kvec[0], -kvec[1], kvec[2]);
			if ( (fabs(fabs(kvec[1])-M_PI) < 0.1) || (fabs(fabs(kvec[2])-M_PI) < 0.01) ) fprintf(fileOut,"\n");
		} for (int tt = DIMXY-1; tt >= 0; tt--) { 
			kvec = &fermiSurface[tt][zz][0];
			fprintf(fileOut,"  %2.6f   %2.6f   %2.6f  \n",  kvec[0], -kvec[1], kvec[2]);
		} for (int tt = 0; tt < DIMXY; tt++) { 
			kvec = &fermiSurface[tt][zz][0]; 
			fprintf(fileOut,"  %2.6f   %2.6f   %2.6f  \n",  kvec[1], -kvec[0], kvec[2]);
			if ( (fabs(fabs(kvec[1])-M_PI) < 0.1) || (fabs(fabs(kvec[2])-M_PI) < 0.01) ) fprintf(fileOut,"\n");
		}
		kvec = &fermiSurface[DIMXY-1][zz][0]; fprintf(fileOut,"  %2.6f   %2.6f   %2.6f  \n",  kvec[1],  kvec[0], kvec[2]);
		fprintf(fileOut,"\n\n");
	}
}



double* vbar(double Kx, double Ky, double Kz) {
	static double vbar[3];
	vbar[0] = vbar[1] = vbar[2] = 0;

	//FILE *fileOut = fopen("amro.dat","w");
	//fprintf(fileOut, "            t            kx            ky            kz         vbarx         vbary         vbarz\n");

	int nn=0; 
	while (exp(-nn*h/tau) > ACC) {
		nn++;
		// compute next point with Runge-Kutta
		double* force;
		force = calculateForce( Kx, Ky, Kz);
		double k1x = +h*force[0];
		double k1y = +h*force[1];
		double k1z = +h*force[2];
		force = calculateForce(Kx+k1x/2., Ky+k1y/2., Kz+k1z/2.);
		double k2x = +h*force[0];
		double k2y = +h*force[1];
		double k2z = +h*force[2];
		force = calculateForce(Kx+k2x/2., Ky+k2y/2., Kz+k2z/2.);
		double k3x = +h*force[0];
		double k3y = +h*force[1];
		double k3z = +h*force[2];
		force = calculateForce(Kx+k3x, Ky+k3y, Kz+k3z);
		double k4x = +h*force[0];
		double k4y = +h*force[1];
		double k4z = +h*force[2];
		
		double newKx = Kx + k1x/6. + k2x/3. + k3x/3. + k4x/6.;
		double newKy = Ky + k1y/6. + k2y/3. + k3y/3. + k4y/6.;
		double newKz = Kz + k1z/6. + k2z/3. + k3z/3. + k4z/6.;
		if (wrapBrillouin){
			if (newKx >    M_PI) {newKx = newKx - 2*M_PI; fprintf(fileOut,"\n");}
			if (newKx <   -M_PI) {newKx = newKx + 2*M_PI; fprintf(fileOut,"\n");}
			if (newKy >    M_PI) {newKy = newKy - 2*M_PI; fprintf(fileOut,"\n");}
			if (newKy <   -M_PI) {newKy = newKy + 2*M_PI; fprintf(fileOut,"\n");}
			if (newKz >  2*M_PI) {newKz = newKz - 4*M_PI; fprintf(fileOut,"\n");}
			if (newKz < -2*M_PI) {newKz = newKz + 4*M_PI; fprintf(fileOut,"\n");}
		}
		double* newVelocity; 
		newVelocity = calculateVelocity(newKx, newKy, newKz);

		vbar[0] += newVelocity[0]*exp(-nn*h/tau);
		vbar[1] += newVelocity[1]*exp(-nn*h/tau);
		vbar[2] += newVelocity[2]*exp(-nn*h/tau);

		//fprintf(fileOut,"% 13f % 13f % 13f % 13f % 13f % 13f % 13f \n", nn*h, newKx, newKy, newKz, vbar[0]/(double)nn, vbar[1]/(double)nn, vbar[2]/(double)nn);
		//printf("iteration %5i  --  t = %4f , vbar = (%4f,%4f,%4f)\n", nn, nn*h, vbar[0]/(double)nn, vbar[1]/(double)nn, vbar[2]/(double)nn);//fflush(stdout);
		
		Kx = newKx;
		Ky = newKy;
		Kz = newKz;
	}
	vbar[0] *= 1/(double)nn;
	vbar[1] *= 1/(double)nn;
	vbar[2] *= 1/(double)nn;

	return vbar;
}









int res = DIM;
double step = 1/(double)DIM;
bool visited[DIM][DIM][DIM];
double treshold = 1.4;
int iter=0;

double spectral(int ii, int jj, int kk){
	double Kx = 2* M_PI* ii* step;
	double Ky = 2* M_PI* jj* step;
	double Kz = 2* M_PI* kk* step;
	double xi = calculateDispersion(Kx,Ky,Kz);
	return 0.05/(xi*xi+0.025);
}

void visitPlane(int kk){

}

void evalPoint(int ii, int jj, int kk){
	if (ii >= res) ii = ii-res;
	if (jj >= res) jj = jj-res;
	if (kk >= res) kk = kk-res;
	if (ii < 0)    ii = ii+res;
	if (jj < 0)    jj = jj+res;
	if (kk < 0)    kk = kk+res;
	
	iter++;
	if ( !(visited[ii][jj][kk]) && (spectral(ii,jj,kk) > treshold)) {
		visited[ii][jj][kk] = true;
		evalPoint(ii-1,jj  ,kk  );
		evalPoint(ii+1,jj  ,kk  );
		evalPoint(ii  ,jj-1,kk  );
		evalPoint(ii  ,jj+1,kk  );
		evalPoint(ii  ,jj  ,kk-1);
		evalPoint(ii  ,jj  ,kk+1);
	}

	return;
}

int main(int argc, const char * argv[]) {
	printf("\namroBuddy starting\n\n");


	double FS[DIMXY][DIMZ][3];
	calculateFS(FS);
	printFS(FS);

	//vbar(FS[10][15][0], FS[10][15][1], FS[10][15][2]);

	printf("\n\namroBuddy over.\n");
	return 0;
}




// 	//// PRE-CALCULATE
//     double sink[DIMXY]; double sin2k[DIMXY];  double sink_2[DIMXY];
//     double cosk[DIMXY]; double cos2k[DIMXY];  double cosk_2[DIMXY];
    
// 	int ii=0; for(ii=0; ii<DIMXY; ii++){
//         double k = M_PI*(ii*1.0/nK);
//         sink[ii] = sin(k); sin2k[ii] = sin(2.*k); sink_2[ii] = sin(k/2.);
//         cosk[ii] = cos(k); cos2k[ii] = cos(2.*k); cosk_2[ii] = cos(k/2.);
//     }

//     double coskz[DIMZ]; double coskz_2[DIMZ]; double sinkz_2[DIMZ];
    
// 	int kk=0; for(kk=0; kk<DIMZ; kk++){
//         double kz = M_PI*(kk*2.0/(nKz)); // 0 to 2Pi assumes parity (period of Markie is 4Pi)
//         coskz[kk] = cos(kz); coskz_2[kk] = cos(kz/2.); sinkz_2[kk] = sin(kz/2.);
//     }

// 	double           kOnGrid [DIMXY][DIMXY][DIMZ][3];
// 	double  dispersionOnGrid [DIMXY][DIMXY][DIMZ];
// 	double    velocityOnGrid [DIMXY][DIMXY][DIMZ][3];

// 	int ii=0; for(ii=0; ii<DIMXY; ii++){
//         double k = M_PI*(ii*1.0/nK);
//         sink[ii] = sin(k); sin2k[ii] = sin(2.*k); sink_2[ii] = sin(k/2.);
//         cosk[ii] = cos(k); cos2k[ii] = cos(2.*k); cosk_2[ii] = cos(k/2.);
//     }
// 	int kk=0; for(kk=0; kk<DIMZ; kk++){
//         double kz = M_PI*(kk*2.0/(nKz)); // 0 to 2Pi assumes parity (period of Markie is 4Pi)
//         coskz[kk] = cos(kz); coskz_2[kk] = cos(kz/2.); sinkz_2[kk] = sin(kz/2.);
//     }


// 	if (tz!=0){
// 		dispersion = -2.*tz*cos(kz/2.);
// 		dispersion *= (cos(kx)-cos(ky))*(cos(kx)-cos(ky));
// 		dispersion *= cos(kx/2.)*cos(ky/2.);
// 	}
// 	dispersion += -mu;
// 	dispersion += -2.*t   * (cos(kx) + cos(ky));
// 	dispersion += -4.*tp  *  cos(kx)*cos(ky);
// 	dispersion += -2.*tpp * (cos(2.*kx) + cos(2.*ky));
// 	return  dispersion;

// 	static double velocity[3];
// 	if (tz!=0){
// 		velocity[0] =  -2.*tz*cos(kz/2.);
// 		velocity[0] *= -2.*(cos(kx)-cos(ky))*sin(kx)*cos(kx/2.)*cos(ky/2.);
// 		velocity[0] *= (cos(kx)-cos(ky))*(cos(kx)-cos(ky))*sin(kx/2.)*cos(ky/2.) /2.;
// 		velocity[1]  = -2.*tz*cos(kz/2.);
// 		velocity[1] *=  2*(cos(kx)-cos(ky))*sin(ky)*cos(kx/2.)*cos(ky/2.);
// 		velocity[1] *=    (cos(kx)-cos(ky))*(cos(kx)-cos(ky))*cos(kx/2.)*sin(ky/2.) /2.;
// 		velocity[2]  = -2.*tz*sin(kz/2.);
// 		velocity[2] *=  (cos(kx)-cos(ky))*(cos(kx)-cos(ky));
// 		velocity[2] *=  cos(kx/2.)*cos(ky/2.);
// 	}
// 	velocity[0] +=  2.*t    *  sin(kx);
// 	velocity[0] +=  4.*tp   *  sin(kx)*cos(ky);
// 	velocity[0] +=  4.*tpp  *  sin(2.*kx);
// 	velocity[1] +=  2.*t    *  sin(ky);
// 	velocity[1] +=  4.*tp   *  sin(ky)*cos(kx);
// 	velocity[1] +=  4.*tpp  *  sin(2.*ky);
// 	return  velocity;
// }
















	// FILE *fileOut = fopen("mdc.out","w");

	// for (int ii = 0; ii<res; ii++){
	// 	for (int jj = 0; jj<res; jj++){
	// 		fprintf(fileOut,"  %5f  ", spectral(ii,jj,0));
	// 	}
	// 	fprintf(fileOut,"\n");
	// }

	// for (int ii = 0; ii<res; ii++){
	// 	for (int jj = 0; jj<res; jj++){
	// 		for (int kk = 0; kk<res; kk++){
	// 			visited[ii][jj][kk]=0;
	// 		}
	// 	}
	// }

	// evalPoint(21,80,0);

	// for (int ii = 0; ii<res; ii++){
	// 	for (int jj = 0; jj<res; jj++){
	// 		kk=0;//for (int kk = 0; kk<res; kk++){
	// 			if (visited[ii][jj][kk]) fprintf(fileOut,"%i\t%i\t%i\n", ii,jj,kk);
	// 		//}
	// 	}
	// }
	// long int total=0;
	// for (int ii = 0; ii<10000; ii++){
	// 	for (int jj = 0; jj<10000; jj++){
	// 		for (int kk = 0; kk<10000; kk++){
	// 			total++;
	// 		}
	// 	}
	// }
	// printf("%li",total);
// 	return 0;
// }
