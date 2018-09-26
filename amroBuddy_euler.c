
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double t   =  1.;
double tp  = -0.3;
double tpp =  0.01;
double mu  = -0.9;
double Bz  =  0.1;
double h   =  0.05;
int    Nt  =  10000;

double kx0 =  M_PI/3.;
double ky0 =  M_PI/2.;

int main(int argc, const char * argv[]) {

	printf("amroBuddy starting\n\n");
	FILE *fileOut = fopen("amro.dat","w");
	fprintf(fileOut, "            t            kx            ky\n");
	fprintf(fileOut,"% 13f % 13f % 13f \n", 0., kx0, ky0);
		
	for (int nn=0; nn<Nt; nn++){
		printf("= %f / %f\n",nn*h,Nt*h);//fflush(stdout);
		
		// double xik  = -mu;
		//        xik += -2.*t   * (cos(kx0) + cos(ky0));
		//        xik += -4.*tp  *  cos(kx0)*cos(ky0);
		//        xik += -2.*tpp * (cos(2.*kx0) + cos(2.*ky0));
		double  vx  =  0;
		        vx +=  2.*t    *  sin(kx0);
		        vx +=  4.*tp   *  sin(kx0)*cos(ky0);
		        vx +=  4.*tpp  *  sin(2.*kx0);
		double  vy  =  0;
		        vy +=  2.*t    *  sin(ky0);
		        vy +=  4.*tp   *  sin(ky0)*cos(kx0);
		        vy +=  4.*tpp  *  sin(2.*ky0);
		double Fx = vy*Bz;
		double Fy = -vx*Bz;

		double kx1 = kx0 + h*Fx;
		double ky1 = ky0 + h*Fy;

		if (kx1 > M_PI)  kx1=kx1-2.*M_PI;
		if (kx1 < -M_PI) kx1=kx1+2.*M_PI;
		if (ky1 > M_PI)  ky1=ky1-2.*M_PI;
		if (ky1 < -M_PI) ky1=ky1+2.*M_PI;

		fprintf(fileOut,"% 13f % 13f % 13f \n", nn*h, kx1, ky1);
		kx0 = kx1;
		ky0 = ky1;
	}
}