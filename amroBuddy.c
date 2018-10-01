#define _USE_MATH_DEFINES
#include <math.h>
#include <complex.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef enum { false, true } bool;

#define ACCURACY 10e-7



///////////////////////////////////// FUNCTIONS FOR READING FILES
void readDouble(FILE * file, char * name,  double * value) {
    rewind(file);
    char tempbuff[200];
    while(!feof(file))
    {
        if (fgets(tempbuff,200,file))
        {
            char tmpstr1[50];
            char tmpstr2[50];
            sscanf(tempbuff, "%49s %49s\n", tmpstr1, tmpstr2);
            if (strcmp(tmpstr1,name)==0) { *value = atof(tmpstr2); return;}
        }
    }
    printf("\ncannot find the %s parameter in 'model.dat'", name);
    exit(1);
}

void readInt(FILE * file, char * name,  int * value) {
    rewind(file);
    char tempbuff[200];
    while(!feof(file))
    {
        if (fgets(tempbuff,200,file))
        {
            char tmpstr1[50];
            char tmpstr2[50];
            sscanf(tempbuff, "%49s %49s\n", tmpstr1, tmpstr2);
            if (strcmp(tmpstr1,name)==0) { *value = atoi(tmpstr2); return;}
        }
    }
    printf("\ncannot find the %s parameter in 'model.dat'", name);
    exit(1);
}
///////////////////////////////////////////////////////////////



double calculateDispersion(double bp[5], double kx, double const ky, double kz){
	double mu  = bp[0];
	double t   = bp[1];
	double tp  = bp[2];
	double tpp = bp[3];
	double tz  = bp[4];
	
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

double *calculateVelocity(double bp[5], double kx, double ky, double kz){
	static double velocity[3];
	double t   = bp[1];
	double tp  = bp[2];
	double tpp = bp[3];
	double tz  = bp[4];

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

double* calculateForce(double bp[5], double kx, double ky, double kz, double Bx, double By, double Bz){
	static double force[3];
	double *velocity; 
	velocity = calculateVelocity(bp, kx,ky,kz);
	force[0]  =  velocity[1]*Bz;
	force[0] += -velocity[2]*By;
	force[1]  = -velocity[0]*Bz;
	force[1] +=  velocity[2]*Bx;
	force[2]  =  velocity[0]*By;
	force[2] += -velocity[1]*Bx;
	return  force;
}


double distance(const double v1[3], const double v2[3]){
	return sqrt((v2[0]-v1[0])*(v2[0]-v1[0]) + (v2[1]-v1[1])*(v2[1]-v1[1]) + (v2[2]-v1[2])*(v2[2]-v1[2]));
}

double* vbar(double bp[5], double Kx, double Ky, double Kz, double dt, double B[3], double tau, double acc, bool printTrigger) {
	static double vbar[3];
	vbar[0] = vbar[1] = vbar[2] = 0;

	FILE *fileOut;
	if (printTrigger) fileOut = fopen("trajectory.dat","w");
	if (printTrigger) fprintf(fileOut, "            t            kx            ky            kz         vbarx         vbary         vbarz\n");

	int nn=0; 
	while (exp(-nn*dt/tau) > acc) {
		nn++;
	// int nn=0;
	// for (nn=0; nn<500; nn++) {
	
		//// compute next point witdt Runge-Kutta
		double* force;
		force = calculateForce(bp,  Kx, Ky, Kz , B[0],B[1],B[2]);
		double k1x = +dt*force[0];
		double k1y = +dt*force[1];
		double k1z = +dt*force[2];
		force = calculateForce(bp, Kx+k1x/2., Ky+k1y/2., Kz+k1z/2., B[0],B[1],B[2]);
		double k2x = +dt*force[0];
		double k2y = +dt*force[1];
		double k2z = +dt*force[2];
		force = calculateForce(bp, Kx+k2x/2., Ky+k2y/2., Kz+k2z/2., B[0],B[1],B[2]);
		double k3x = +dt*force[0];
		double k3y = +dt*force[1];
		double k3z = +dt*force[2];
		force = calculateForce(bp, Kx+k3x, Ky+k3y, Kz+k3z, B[0],B[1],B[2]);
		double k4x = +dt*force[0];
		double k4y = +dt*force[1];
		double k4z = +dt*force[2];
		
		double newKx = Kx + k1x/6. + k2x/3. + k3x/3. + k4x/6.;
		double newKy = Ky + k1y/6. + k2y/3. + k3y/3. + k4y/6.;
		double newKz = Kz + k1z/6. + k2z/3. + k3z/3. + k4z/6.;
		
		double* newVelocity; 
		newVelocity = calculateVelocity(bp, newKx, newKy, newKz);
		vbar[0] += newVelocity[0]*exp(-nn*dt/tau);
		vbar[1] += newVelocity[1]*exp(-nn*dt/tau);
		vbar[2] += newVelocity[2]*exp(-nn*dt/tau);
		Kx = newKx;
		Ky = newKy;
		Kz = newKz;
		
		if (printTrigger){
			if (newKx >    M_PI) {newKx = newKx - 2*M_PI; fprintf(fileOut,"\n");}
			if (newKx <   -M_PI) {newKx = newKx + 2*M_PI; fprintf(fileOut,"\n");}
			if (newKy >    M_PI) {newKy = newKy - 2*M_PI; fprintf(fileOut,"\n");}
			if (newKy <   -M_PI) {newKy = newKy + 2*M_PI; fprintf(fileOut,"\n");}
			if (newKz >  2*M_PI) {newKz = newKz - 4*M_PI; fprintf(fileOut,"\n");}
			if (newKz < -2*M_PI) {newKz = newKz + 4*M_PI; fprintf(fileOut,"\n");}
			fprintf(fileOut,"% 13f % 13f % 13f % 13f % 13f % 13f % 13f \n", nn*dt, newKx, newKy, newKz, vbar[0]/(double)nn, vbar[1]/(double)nn, vbar[2]/(double)nn);
		}
	}
	vbar[0] *= 1*dt;
	vbar[1] *= 1*dt;
	vbar[2] *= 1*dt;

	return vbar;
}


int main(int argc, const char * argv[]) {
	printf("\namroBuddy starting\n\n");

	//// PARAMETERS
	double t   =  1.;
	double tp  = -0.12;
	double tpp =  0.06;
	double tz  =  0.07;
	double mu  = -0.71;
	double Bamp = 0.03;
	double Bphi = 0;//M_PI/6.;
	double Btheta_min = 0;
	double Btheta_max = 160;
	double Btheta_step = 5;
	double tau =  25;
	double accuracy = 10e-9;
	double dt = 1.;
	int Nz = 30;
	int Nxy = 20;
	
	//// READING PARAMETERS
    FILE * file = fopen("model.dat", "rt");
    if(file == NULL) {printf("file %s not found", "model.dat"); exit(1);}
    printf("reading parameters from model.dat\n\n") ;
    readDouble(file, "mu",              &mu);
    readDouble(file, "t",               &t);
    readDouble(file, "tp",              &tp);
    readDouble(file, "tpp",             &tpp);
    readDouble(file, "tz",              &tz);
    readDouble(file, "Bamp",            &Bamp);
    readDouble(file, "Bphi",            &Bphi);
    readDouble(file, "Btheta_min",      &Btheta_min);
    readDouble(file, "Btheta_max",      &Btheta_max);
    readDouble(file, "Btheta_step",     &Btheta_step);
    readDouble(file, "tau",             &tau);
    readDouble(file, "accuracy",        &accuracy);
    readDouble(file, "dt",              &dt);
    readInt(file, "Nxy",    &Nxy);
    readInt(file, "Nz",     &Nz);
    fclose(file);

	double bandParameters[5] = {mu,t,tp,tpp,tz};
	int Ntheta = (Btheta_max - Btheta_min)/Btheta_step;
	

	//////////////////////////////// CALCULATE FS
	double FS[Nxy][Nz][3];
	for (int zz = 0; zz < Nz; zz++) {
		double kFz = -2*M_PI + (zz/((double)Nz-1))*4*M_PI;
		
		//// Only an eight of the FS is computed. theta is clockwise, from (pi,pi), which means that
		//// theta =0 is the axis from (0,0) to (pi,pi) and =pi/4 is the axis from (0,0) to (0,pi)
		//// the max value of angle is determined for the FS for kx = 0.
		double maxTheta = M_PI/4.; 
		double maxKy = M_PI;
		do {
			maxTheta -= ACCURACY;
			maxKy = M_PI*(1+tan(maxTheta-M_PI/4.));
		} while (calculateDispersion(bandParameters, 0, maxKy, kFz) > 0.);
		
		for (int tt = 0; tt < Nxy; tt++) {			
			double theta = (tt/((double)Nxy-1)) * maxTheta;
			double kx1 = 0   , ky1 = M_PI*(1+tan(theta-M_PI/4.));
			double kx2 = M_PI, ky2 = M_PI;
			double xi1 = calculateDispersion(bandParameters, kx1, ky1, kFz);
			double xi2 = calculateDispersion(bandParameters, kx2, ky2, kFz);
			double xi, kFx = NAN, kFy = NAN;		

			if (xi1*xi2 < 0) {
				do {
					kFx = (kx2+kx1)/2.;
					kFy = (ky2+ky1)/2.;
					xi = calculateDispersion(bandParameters, kFx, kFy, kFz);
					if (xi > 0) {
						kx2 = kFx;
						ky2 = kFy;
					} else if (xi < 0) {
						kx1 = kFx;
						ky1 = kFy;
					}
				} while (fabs(xi) > ACCURACY);
			}
			FS[tt][zz][0] = kFx;
			FS[tt][zz][1] = kFy;
			FS[tt][zz][2] = kFz;
		}
	}
	//////////////////////////////////////////////////

	///////////////////////////// PRINT FS
	FILE *fsFileout = fopen("FS.dat","w");
	fprintf(fsFileout, "  kx         ky         kz\n");
	for (int zz = 0; zz < Nz; zz++) {
		double* kvec;
		for (int tt = Nxy-1; tt >= 0; tt--) { 
			kvec = &FS[tt][zz][0]; 
			fprintf(fsFileout,"  %2.6f   %2.6f   %2.6f  \n",  kvec[1],  kvec[0], kvec[2]);
		} for (int tt = 0; tt < Nxy; tt++) { 
			kvec = &FS[tt][zz][0];
			fprintf(fsFileout,"  %2.6f   %2.6f   %2.6f  \n",  kvec[0],  kvec[1], kvec[2]);
			if ( (fabs(fabs(kvec[1])-M_PI) < 0.1) || (fabs(fabs(kvec[2])-M_PI) < 0.01) ) fprintf(fsFileout,"\n");
		} for (int tt = Nxy-1; tt >= 0; tt--){
			kvec = &FS[tt][zz][0];
			fprintf(fsFileout,"  %2.6f   %2.6f   %2.6f  \n", -kvec[0],  kvec[1], kvec[2]);
		} for (int tt = 0; tt < Nxy; tt++) { 
			kvec = &FS[tt][zz][0]; 
			fprintf(fsFileout,"  %2.6f   %2.6f   %2.6f  \n", -kvec[1],  kvec[0], kvec[2]);
			if ( (fabs(fabs(kvec[1])-M_PI) < 0.1) || (fabs(fabs(kvec[2])-M_PI) < 0.01) ) fprintf(fsFileout,"\n");
		} for (int tt = Nxy-1; tt >= 0; tt--) { 
			kvec = &FS[tt][zz][0];
			fprintf(fsFileout,"  %2.6f   %2.6f   %2.6f  \n", -kvec[1], -kvec[0], kvec[2]);
		} for (int tt = 0; tt < Nxy; tt++) { 
			kvec = &FS[tt][zz][0];
			fprintf(fsFileout,"  %2.6f   %2.6f   %2.6f  \n", -kvec[0], -kvec[1], kvec[2]);
			if ( (fabs(fabs(kvec[1])-M_PI) < 0.1) || (fabs(fabs(kvec[2])-M_PI) < 0.01) ) fprintf(fsFileout,"\n");
		} for (int tt = Nxy-1; tt >= 0; tt--) { 
			kvec = &FS[tt][zz][0];
			fprintf(fsFileout,"  %2.6f   %2.6f   %2.6f  \n",  kvec[0], -kvec[1], kvec[2]);
		} for (int tt = 0; tt < Nxy; tt++) { 
			kvec = &FS[tt][zz][0]; 
			fprintf(fsFileout,"  %2.6f   %2.6f   %2.6f  \n",  kvec[1], -kvec[0], kvec[2]);
			if ( (fabs(fabs(kvec[1])-M_PI) < 0.1) || (fabs(fabs(kvec[2])-M_PI) < 0.01) ) fprintf(fsFileout,"\n");
		}
		kvec = &FS[Nxy-1][zz][0]; fprintf(fsFileout,"  %2.6f   %2.6f   %2.6f  \n",  kvec[1],  kvec[0], kvec[2]);
		fprintf(fsFileout,"\n\n");
	}


	///////////////////////////////////////////////// INTEGRAL


	FILE *fileOut = fopen("amro.dat","w");
	fprintf(fileOut,"B           theta       phi         sigma_zz    sigma_ref   \n");fflush(fileOut);


	double integral = 0;
	double integralRef = 0;
	for (int bb=0; bb<Ntheta; bb++) {
		double Btheta = bb*M_PI/((double)Ntheta-1);
		double Bfield[3];
		Bfield[0] = Bamp*sin(Btheta)*cos(Bphi);
		Bfield[1] = Bamp*sin(Btheta)*sin(Bphi);
		Bfield[2] = Bamp*cos(Btheta);
		
		integral =0;
		for (int zz = 0; zz < Nz; zz++) {
			
			const double* kvec;
			for (int tt = 0; tt < Nxy; tt++) {
				if (!isnan(FS[tt][zz][0])) {
					kvec = &FS[tt][zz][0];
					
					// WARNING: using the following definition can lead to segmentation faults if tt=0 or tt=Nxy 
					double* next = &FS[(tt+1)][zz][0];
					double* prev = &FS[(tt-1)][zz][0];
					
					double len = sqrt((M_PI-kvec[0])*(M_PI-kvec[0])+(M_PI-kvec[1])*(M_PI-kvec[1]))*2*sin(M_PI/8./(double)Nxy);
					if (tt == 0 || isnan(*prev)) {
						if (!isnan(*next) && (tt != Nxy-1) ) len = distance(kvec,next)/2.;
					} else if (tt == Nxy-1 || isnan(*next)) {
						if (!isnan(*prev)) len = distance(prev,kvec)/2.;
					} else {                   
						len = distance(prev,kvec)/2.+distance(kvec,next)/2.;
					}
					if (len > 10) {
						printf("zz = %6i, kk= %6i\n", zz,tt);
						printf("prev = %6f,%6f,%6f \n", *prev,prev[1],prev[2]);
						printf("kvec = %6f,%6f,%6f \n", kvec[0],kvec[1],kvec[2]);
						printf("kvec = %6f,%6f,%6f \n", FS[tt][zz][0],FS[tt][zz][1],FS[tt][zz][2]);
						printf("next = %6f,%6f,%6f \n", next[0],next[1],next[2]);
						printf("\n\n");
					}
					//double fac = 1/(double)(Nxy*Nz);
					
					double fac = len / (double)Nz;
					
					double vz = calculateVelocity(bandParameters, kvec[1],  kvec[0], kvec[2])[2];
					bool printTrigger = (bb==Ntheta/4 && zz==Nz/2 && tt==Nxy/2); 
					integral += fac * vz * vbar(bandParameters,  kvec[1],  kvec[0], kvec[2], dt , Bfield, tau, accuracy, printTrigger)[2];
					integral += fac * vz * vbar(bandParameters,  kvec[0],  kvec[1], kvec[2], dt , Bfield, tau, accuracy, false       )[2];
					integral += fac * vz * vbar(bandParameters, -kvec[0],  kvec[1], kvec[2], dt , Bfield, tau, accuracy, false       )[2];
					integral += fac * vz * vbar(bandParameters, -kvec[1],  kvec[0], kvec[2], dt , Bfield, tau, accuracy, false       )[2];
					integral += fac * vz * vbar(bandParameters, -kvec[1], -kvec[0], kvec[2], dt , Bfield, tau, accuracy, false       )[2];
					integral += fac * vz * vbar(bandParameters, -kvec[0], -kvec[1], kvec[2], dt , Bfield, tau, accuracy, false       )[2];
					integral += fac * vz * vbar(bandParameters,  kvec[0], -kvec[1], kvec[2], dt , Bfield, tau, accuracy, false       )[2];
					integral += fac * vz * vbar(bandParameters,  kvec[1], -kvec[0], kvec[2], dt , Bfield, tau, accuracy, false       )[2];
					;
				}
			}
			//printf("bb=%i/%i -- zz=%i/%i -- integral = %f\n",bb,Ntheta, zz, Nz, integral); //fflush(stdout);
		}
		printf("bb=%i/%i -- sigma_zz = %f\n",bb,Ntheta, integral); //fflush(stdout);
		
		integral /= Nz*8.;
		if (bb==0) {integralRef = integral;}
		fprintf(fileOut,"%2.9f %2.9f %2.9f %2.9f %2.9f\n",Bamp, Btheta, Bphi, integral, integralRef); // fflush(fileOut);
	}
	fprintf(fileOut,"\n");

	printf("\n\namroBuddy over.\n");
	return 0;
}
