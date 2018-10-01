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
double mu  = -0.81;
double tau =  25;

// double B = 0.006;
// double theta = 0;
// double phi = 0;

// double Bx = 0.01;//B*sin(theta)*cos(phi);
// double By = 0.01;//B*sin(theta)*sin(phi);
// double Bz = 0.05;//B*cos(theta);
//double h  = 0.1;// 
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

double* calculateForce(double kx, double ky, double kz, double Bx, double By, double Bz){
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
		
		//// Only an eight of the FS is computed. theta is clockwise, from (pi,pi), which means that
		//// theta =0 is the axis from (0,0) to (pi,pi) and =pi/4 is the axis from (0,0) to (0,pi)
		//// the max value of angle is determined for the FS for kx = 0.
		double maxTheta = M_PI/4.; 
		double maxKy = M_PI;
		do {
			maxTheta -= ACC;
			maxKy = M_PI*(1+tan(maxTheta-M_PI/4.));
		} while (calculateDispersion(0, maxKy, kFz) > 0.);
		
		for (int tt = 0; tt < DIMXY; tt++) {			
			double theta = (tt/((double)DIMXY-1)) * maxTheta;
			double kx1 = 0   , ky1 = M_PI*(1+tan(theta-M_PI/4.));
			double kx2 = M_PI, ky2 = M_PI;
			double xi1 = calculateDispersion(kx1, ky1, kFz);
			double xi2 = calculateDispersion(kx2, ky2, kFz);
			double xi, kFx = NAN, kFy = NAN;		

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



double* vbar(double Kx, double Ky, double Kz, double h, double B[3], bool printTrigger) {
	static double vbar[3];
	vbar[0] = vbar[1] = vbar[2] = 0;

	FILE *fileOut;
	if (printTrigger) fileOut = fopen("trajectory.dat","w");
	if (printTrigger) fprintf(fileOut, "            t            kx            ky            kz         vbarx         vbary         vbarz\n");

	int nn=0; 
	while (exp(-nn*h/tau) > ACC) {
		nn++;
		//// compute next point with Runge-Kutta
		double* force;
		force = calculateForce( Kx, Ky, Kz , B[0],B[1],B[2]);
		double k1x = +h*force[0];
		double k1y = +h*force[1];
		double k1z = +h*force[2];
		force = calculateForce(Kx+k1x/2., Ky+k1y/2., Kz+k1z/2., B[0],B[1],B[2]);
		double k2x = +h*force[0];
		double k2y = +h*force[1];
		double k2z = +h*force[2];
		force = calculateForce(Kx+k2x/2., Ky+k2y/2., Kz+k2z/2., B[0],B[1],B[2]);
		double k3x = +h*force[0];
		double k3y = +h*force[1];
		double k3z = +h*force[2];
		force = calculateForce(Kx+k3x, Ky+k3y, Kz+k3z, B[0],B[1],B[2]);
		double k4x = +h*force[0];
		double k4y = +h*force[1];
		double k4z = +h*force[2];
		
		double newKx = Kx + k1x/6. + k2x/3. + k3x/3. + k4x/6.;
		double newKy = Ky + k1y/6. + k2y/3. + k3y/3. + k4y/6.;
		double newKz = Kz + k1z/6. + k2z/3. + k3z/3. + k4z/6.;
		if (wrapBrillouin && printTrigger){
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

		if (printTrigger) fprintf(fileOut,"% 13f % 13f % 13f % 13f % 13f % 13f % 13f \n", nn*h, newKx, newKy, newKz, vbar[0]/(double)nn, vbar[1]/(double)nn, vbar[2]/(double)nn);
		Kx = newKx;
		Ky = newKy;
		Kz = newKz;
	}
	vbar[0] *= 1/(double)nn;
	vbar[1] *= 1/(double)nn;
	vbar[2] *= 1/(double)nn;

	return vbar;
}

double distance(double v1[3], double v2[3]){
	return sqrt((v2[0]-v1[0])*(v2[0]-v1[0]) + (v2[1]-v1[1])*(v2[1]-v1[1]) + (v2[2]-v1[2])*(v2[2]-v1[2]));
}

int main(int argc, const char * argv[]) {
	printf("\namroBuddy starting\n\n");

	double FS[DIMXY][DIMZ][3];
	calculateFS(FS);
	printFS(FS);

	double hh = 5.;
	double integral=0;
	double integralRef = 0;
	double Bamp = 0.03;
	double Bphi = 0;//M_PI/6.;
	int fieldSamples = 31;

	FILE *fileOut = fopen("amro.dat","w");
	fprintf(fileOut,"B           theta       phi         sigma_zz    sigma_ref   \n");fflush(fileOut);

	for (int bb=0; bb<fieldSamples; bb++) {
		double Btheta = bb*M_PI/((double)fieldSamples-1);
		double Bfield[3];
		Bfield[0] = Bamp*sin(Btheta)*cos(Bphi);
		Bfield[1] = Bamp*sin(Btheta)*sin(Bphi);
		Bfield[2] = Bamp*cos(Btheta);
		
		integral =0;
		for (int zz = 0; zz < DIMZ; zz++) {
			
			double* kvec;
			for (int tt = 0; tt < DIMXY; tt++) {
				if (!isnan(FS[tt][zz][0])) {
					kvec = &FS[tt][zz][0];
					//// WARNING: using the following definition can lead to segmentation faults if tt=0 or tt=DIMXY 
					double* next = &FS[(tt+1)][zz][0];
					double* prev = &FS[(tt-1)][zz][0];
					
					double len = sqrt((M_PI-kvec[0])*(M_PI-kvec[0])+(M_PI-kvec[1])*(M_PI-kvec[1]))*2*sin(M_PI/8./(double)DIMXY);
					if (tt == 0 || isnan(*prev)) {
						if (!isnan(*next)) len = distance(prev,kvec)/2.;
					} else if (tt == DIMXY-1 || isnan(*next)) {
						if (!isnan(*prev)) len = distance(prev,kvec)/2.;
					} else {                   
						len = distance(prev,kvec)/2.+distance(kvec,next)/2.;
					}
					
					double vz = calculateVelocity(kvec[1],  kvec[0], kvec[2])[2];
					bool printTrigger = (bb==fieldSamples/4 && zz==DIMZ/2 && tt==DIMXY/2); 
					integral += vz * vbar( kvec[1],  kvec[0], kvec[2] , hh , Bfield, printTrigger)[2];
					integral += vz * vbar( kvec[0],  kvec[1], kvec[2] , hh , Bfield, false)[2];
					integral += vz * vbar(-kvec[0],  kvec[1], kvec[2] , hh , Bfield, false)[2];
					integral += vz * vbar(-kvec[1],  kvec[0], kvec[2] , hh , Bfield, false)[2];
					integral += vz * vbar(-kvec[1], -kvec[0], kvec[2] , hh , Bfield, false)[2];
					integral += vz * vbar(-kvec[0], -kvec[1], kvec[2] , hh , Bfield, false)[2];
					integral += vz * vbar( kvec[0], -kvec[1], kvec[2] , hh , Bfield, false)[2];
					integral += vz * vbar( kvec[1], -kvec[0], kvec[2] , hh , Bfield, false)[2];
					//integral *= len;
					//printf("%f\n",len);
				}
			}
			printf("bb=%i/%i -- zz=%i/%i -- integral = %f\n",bb,51, zz, DIMZ, integral); fflush(stdout);
		}
		integral /= DIMZ*8.;
		if (bb==0) {integralRef = integral;}
		fprintf(fileOut,"%2.9f %2.9f %2.9f %2.9f %2.9f\n",Bamp, Btheta, Bphi, integral, integralRef);fflush(fileOut);
	}
	fprintf(fileOut,"\n");

	printf("\n\namroBuddy over.\n");
	return 0;
}
