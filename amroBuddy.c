
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
double h   =  2.;
int    Nt  =  100;
double tau =  25;

double Kx =  M_PI;
double Ky =  M_PI/10.;

double xi(double kx, double ky){
	double xi  = -mu;
		   xi += -2.*t   * (cos(Kx) + cos(Ky));
		   xi += -4.*tp  *  cos(Kx)*cos(Ky);
		   xi += -2.*tpp * (cos(2.*Kx) + cos(2.*Ky));
	return xi;
}		

double vx(double kx,double ky){
	double  vx  =  0;
			vx +=  2.*t    *  sin(kx);
			vx +=  4.*tp   *  sin(kx)*cos(ky);
			vx +=  4.*tpp  *  sin(2.*kx);
	return vx;
}

double vy(double kx,double ky){
	double  vy  =  0;
			vy +=  2.*t    *  sin(ky);
			vy +=  4.*tp   *  sin(ky)*cos(kx);
			vy +=  4.*tpp  *  sin(2.*ky);
	return vy;
}


int main(int argc, const char * argv[]) {

	printf("amroBuddy starting\n\n");
	FILE *fileOut = fopen("amro.dat","w");
	fprintf(fileOut, "            t            kx            ky         vbarx         vbary\n");
	fprintf(fileOut,"% 13f % 13f % 13f % 13f % 13f\n", 0., Kx, Ky, vx(Kx, Ky), vy(Kx, Ky) );
		
	double vbarx = 0;
	double vbary = 0; 
	for (int nn=0; nn<Nt; nn++){
		

		double k1x = +h*vy( Kx          , Ky          )*Bz;
		double k1y = -h*vx( Kx          , Ky          )*Bz;
		// compute next point with Runge-Kutta
		double k2x = +h*vy( Kx + k1x/2. , Ky + k1y/2. )*Bz;
		double k2y = -h*vx( Kx + k1x/2. , Ky + k1y/2. )*Bz;
		double k3x = +h*vy( Kx + k2x/2. , Ky + k2y/2. )*Bz;
		double k3y = -h*vx( Kx + k2x/2. , Ky + k2y/2. )*Bz;
		double k4x = +h*vy( Kx + k3x    , Ky + k1y    )*Bz;
		double k4y = -h*vx( Kx + k3x    , Ky + k1y    )*Bz;
		double newKx = Kx + k1x/6. + k2x/3. + k3x/3. + k4x/6.;
		double newKy = Ky + k1y/6. + k2y/3. + k3y/3. + k4y/6.;
		
		vbarx += vx(newKx, newKy)*exp(-nn*h/tau);
		vbary += vy(newKx, newKy)*exp(-nn*h/tau);

		if (newKx >  M_PI) newKx = newKx - 2*M_PI;
		if (newKx < -M_PI) newKx = newKx + 2*M_PI;
		if (newKy >  M_PI) newKy = newKy - 2*M_PI;
		if (newKy < -M_PI) newKy = newKy + 2*M_PI;

		fprintf(fileOut,"% 13f % 13f % 13f % 13f % 13f \n", nn*h, newKx, newKy, vbarx, vbary);
		Kx = newKx;
		Ky = newKy;


		printf("t = % 4f / % 4f , vbar = (% 4f,% 4f)\n", nn*h, Nt*h, vbarx, vbary);//fflush(stdout);
	}
}