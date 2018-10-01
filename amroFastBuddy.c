#define _USE_MATH_DEFINES
#include <math.h> 
#include <complex.h> 
#include <sys/stat.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 


typedef enum { false, true } bool;

#define DIM 100
#define ACC 10e-7
#define DIMXY 7
#define DIMZ 50
#define NNMAX 100000

double t   =  1.;
double tp  = -0.12;
double tpp =  0.06;
double tz  =  0.07;
double mu  = -0.41;

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

double* lorentzForce(double kx, double ky, double kz, double Bx, double By, double Bz){
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


double distance(double v1[3], double v2[3]){
	return sqrt((v2[0]-v1[0])*(v2[0]-v1[0]) + (v2[1]-v1[1])*(v2[1]-v1[1]) + (v2[2]-v1[2])*(v2[2]-v1[2]));
}


int main(int argc, const char * argv[]) {
	printf("\nfastBuddy starting\n\n");
	
	int NZ=401;
	double dt = 2.;
	double tau =  25;
	
	double integral=0;
	double integralRef = 0;
	FILE *integralOut;
	integralOut = fopen("fastAMRO.dat","w");
	fprintf(integralOut,"B           theta       phi         sigma_zz    sigma_ref   \n");fflush(integralOut);
	
	FILE *fileOut;
	fileOut = fopen("fastTrajectory.dat","w");
	fprintf(fileOut, "            t            kx            ky            kz\n");
			

	//// precalculate lifetime
	double lifeProbability[NNMAX];
	int lifeSize = 0;
	int nn=1; lifeProbability[0] = 1;
	while (lifeProbability[nn-1] > ACC) {
		lifeProbability[nn] = exp(-nn*dt/tau);
		lifeSize++;
		nn++;
		if (nn == NNMAX) {printf("exceeded maximum enabled length for lifetime\n\n"); exit(0);}
	}

	double Bamp = 0.03;
	double Bphi = 0;//M_PI/6.;
	int fieldSamples = 31;
	for (int bb=0; bb<fieldSamples; bb++) {
		double Btheta = bb * M_PI/ ((double)fieldSamples-1);
		printf("Field at theta = %5.2f\n", 180*Btheta/2./M_PI );fflush(stdout);

		double B[3];
		B[0] = Bamp*sin(Btheta)*cos(Bphi);
		B[1] = Bamp*sin(Btheta)*sin(Bphi);
		B[2] = Bamp*cos(Btheta);
		
		for (int zz=0; zz < NZ; zz++) {
			double dkz = 4*M_PI/(double)NZ; 
			double kz0 = -2*M_PI + zz * dkz; 
			double kx0=NAN; //
			double ky0=NAN; // needs to be the Fermi surface
			
			//// finds the Fermi surface
			double kxBelow = 0, kxAbove = M_PI;
			double kyBelow = 0, kyAbove = M_PI;
			double xi;
			do {
				kx0 = (kxAbove + kxBelow) /2.;
				ky0 = (kyAbove + kyBelow) /2.;
				xi = calculateDispersion(kx0, ky0, kz0);
				if (xi > 0) {
					kxAbove = kx0;
					kyAbove = ky0;
				} else if (xi < 0) {
					kxBelow = kx0;
					kyBelow = ky0;
				}
			} while (fabs(xi) > ACC);

			//// build trajectory from k0
			int nn = 0;
			int trajectorySize = 0;
			double kTrajectory[NNMAX][3];
			double Kx = kx0;
			double Ky = ky0;
			double Kz = kz0;
			//double angularProgression = 0;
			bool   EndOfOrbit = false;
			
			kTrajectory[nn][0] = Kx;
			kTrajectory[nn][1] = Ky;
			kTrajectory[nn][2] = Kz;

			while (!EndOfOrbit) {
				nn++;
				trajectorySize++;
				if (nn == NNMAX) {printf("exceeded maximum enabled length for trajectory\n\n"); exit(0);}
				
				
				//// next point in trajectory from Runge-Kutta
				double* force; //=dk/dt
				force = lorentzForce(Kx       , Ky       , Kz       , B[0],B[1],B[2]);
				double k1x = +dt * force[0];
				double k1y = +dt * force[1];
				double k1z = +dt * force[2];
				force = lorentzForce(Kx+k1x/2., Ky+k1y/2., Kz+k1z/2., B[0],B[1],B[2]);
				double k2x = +dt * force[0];
				double k2y = +dt * force[1];
				double k2z = +dt * force[2];
				force = lorentzForce(Kx+k2x/2., Ky+k2y/2., Kz+k2z/2., B[0],B[1],B[2]);
				double k3x = +dt * force[0];
				double k3y = +dt * force[1];
				double k3z = +dt * force[2];
				force = lorentzForce(Kx+k3x   , Ky+k3y   , Kz+k3z   , B[0],B[1],B[2]);
				double k4x = +dt * force[0];
				double k4y = +dt * force[1];
				double k4z = +dt * force[2];
				
				Kx = Kx + k1x/6. + k2x/3. + k3x/3. + k4x/6.;
				Ky = Ky + k1y/6. + k2y/3. + k3y/3. + k4y/6.;
				Kz = Kz + k1z/6. + k2z/3. + k3z/3. + k4z/6.;
				
				kTrajectory[nn][0] = Kx;
				kTrajectory[nn][1] = Ky;
				kTrajectory[nn][2] = Kz;

				double wrappedK[3]= {Kx,Ky,Kz};;
				if (Kx >    M_PI) wrappedK[0] = Kx - 2*M_PI;
				if (Kx <   -M_PI) wrappedK[0] = Kx + 2*M_PI;
				if (Ky >    M_PI) wrappedK[1] = Ky - 2*M_PI;
				if (Ky <   -M_PI) wrappedK[1] = Ky + 2*M_PI;
				if (Kz >  2*M_PI) wrappedK[2] = Kz - 4*M_PI;
				if (Kz < -2*M_PI) wrappedK[2] = Kz + 4*M_PI;

				// end of orbit?
				EndOfOrbit = (distance(wrappedK, kTrajectory[0]) < distance(kTrajectory[1], kTrajectory[0]));
			}
		
			// print one of the trajectories
			if (bb == fieldSamples/3 && zz == NZ/2){
			for (nn=0; nn<trajectorySize; nn++){
				fprintf(fileOut,"% 13f % 13f % 13f % 13f\n", nn*dt, kTrajectory[nn][0], kTrajectory[nn][1], kTrajectory[nn][2]);
				//if (distance(kTrajectory[nn+1],kTrajectory[nn]) > 0.5) fprintf(fileOut, "\n");
			}
			fprintf(fileOut, "\n\n");
			}

			//// precalculate velocity and lenght
			double velocityTrajectory[trajectorySize][3];
			double dlengthTrajectory[trajectorySize];
			double totalLenght = 0;

			double* velocity;
			for (nn=0; nn<trajectorySize; nn++){
				
				velocity = calculateVelocity(kTrajectory[nn][0],kTrajectory[nn][1],kTrajectory[nn][2]);
				velocityTrajectory[nn][0] = velocity[0];
				velocityTrajectory[nn][1] = velocity[1];
				velocityTrajectory[nn][2] = velocity[2];

				int nPrev = (nn-1)%trajectorySize;
				int nNext = (nn+1)%trajectorySize;
				dlengthTrajectory[nn] = distance(kTrajectory[nNext], kTrajectory[nPrev]);
				totalLenght += dlengthTrajectory[nn];
			}
	
			//// MAIN INTEGRAL the zz resistivity is chosen by velo[2] * velo [2]
			integral = 0; 
			for (nn=0; nn < trajectorySize; nn++){
				for (int mm=0; mm < lifeSize; mm++) {
					int ll = mm % trajectorySize;
					integral += velocityTrajectory[nn][2] * lifeProbability[mm] * velocityTrajectory[ll][2] * dlengthTrajectory[ll] / totalLenght / NZ;
				}
			}
		}
	if (bb==0) {integralRef = integral;}
	fprintf(integralOut,"%2.9f %2.9f %2.9f %2.9f %2.9f\n",Bamp, Btheta, Bphi, integral, integralRef);fflush(integralOut);	
	}
	// fprintf(fileOut,"\n");

	printf("\n\nfastBuddy over.\n");
	return 0;
}




//////////// TRASH

// double dkx1 = kTrajectory[nn-1][0] - kTrajectory[nn-2][0];
// double dky1 = kTrajectory[nn-1][1] - kTrajectory[nn-2][1];
// double dkz1 = kTrajectory[nn-1][2] - kTrajectory[nn-2][2];
// double dkx2 = Kx - kTrajectory[nn-1][0];
// double dky2 = Ky - kTrajectory[nn-1][1];
// double dkz2 = Kz - kTrajectory[nn-1][2];

// double dot    =  dkx1*dkx2 + dky1*dky2 + dkz1*dkz2;
// int dotsign   = (dot >= 0)? +1:-1; // small angle =+1, large=-1
// double crossX =  dky1*dkz2 - dkz1*dky2;
// double crossY = -dkx1*dkz2 + dkz1*dkx2;
// double crossZ =  dkx1*dky2 - dky1*dkx2;
// double crosslen = (crossX*B[0]+crossY*B[1]+crossZ*B[2])/Bamp; // since it is perp to B
// int crossign  = (crosslen >= 0)? +1:-1; // righthand movevement =+1, lefthand=-1

// double angularChange = crossign * atan(crosslen/dot);

