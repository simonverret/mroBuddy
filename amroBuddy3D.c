
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double t   =  1.;
double tp  = -0.12;
double tpp =  0.06;
double tz  =  0.07;
double mu  = -0.9;
double h   =  0.5;
int    Nt  =  500;
double tau =  25;

double Kx =  M_PI/1.;
double Ky =  M_PI/45.;
double Kz =  0.1;

double B = 0.1;
double theta = 0;
double phi = 0;

double Bx = 0.06;//B*sin(theta)*cos(phi);
double By = 0.02;//B*sin(theta)*sin(phi);
double Bz = 0.05;//B*cos(theta);

int wrapBrillouin=1;

double xi(double kx, double ky, double kz){
	double  xi  = -2.*tz*cos(kz/2.);
			xi *= (cos(kx)-cos(ky))*(cos(kx)-cos(ky));
			xi *= cos(kx/2.)*cos(ky/2.);
			// above this line: 3d, below this line: 2d
			xi += -mu;
			xi += -2.*t   * (cos(kx) + cos(ky));
			xi += -4.*tp  *  cos(kx)*cos(ky);
			xi += -2.*tpp * (cos(2.*kx) + cos(2.*ky));
	return  xi;
}
double vx(double kx, double ky, double kz){
	double  vx =  -2.*tz*cos(kz/2.);
			vx *= -2.*(cos(kx)-cos(ky))*sin(kx)*cos(kx/2.)*cos(ky/2.);
			vx *= (cos(kx)-cos(ky))*(cos(kx)-cos(ky))*sin(kx/2.)*cos(ky/2.) /2.;
			// above this line: 3d, below this line: 2d
			vx +=  2.*t    *  sin(kx);
			vx +=  4.*tp   *  sin(kx)*cos(ky);
			vx +=  4.*tpp  *  sin(2.*kx);
	return  vx;
}
double vy(double kx,double ky, double kz){
	double  vy  = -2.*tz*cos(kz/2.);
			vy *=  2*(cos(kx)-cos(ky))*sin(ky)*cos(kx/2.)*cos(ky/2.);
		    vy *=    (cos(kx)-cos(ky))*(cos(kx)-cos(ky))*cos(kx/2.)*sin(ky/2.) /2.;
			// above this line: 3d, below this line: 2d
			vy +=  2.*t    *  sin(ky);
			vy +=  4.*tp   *  sin(ky)*cos(kx);
			vy +=  4.*tpp  *  sin(2.*ky);
	return  vy;
}
double vz(double kx, double ky, double kz){
	double  vz  = -2.*tz*sin(kz/2.);
			vz *=  (cos(kx)-cos(ky))*(cos(kx)-cos(ky));
			vz *=  cos(kx/2.)*cos(ky/2.);
	return vz;
}

double lorentzFx(double kx, double ky, double kz){
	double  forcex  =  vy(kx,ky,kz)*Bz;
			forcex += -vz(kx,ky,kz)*By;
	return  forcex;
}
double lorentzFy(double kx, double ky, double kz){
	double  forcey  = -vx(kx,ky,kz)*Bz;
			forcey +=  vz(kx,ky,kz)*Bx;
	return  forcey;
}
double lorentzFz(double kx, double ky, double kz){
	double  forcez  =  vx(kx,ky,kz)*By;
			forcez += -vy(kx,ky,kz)*Bx;
	return  forcez;
}

void vbar() {
	printf("amroBuddy starting\n\n");

	FILE *fileOut = fopen("amro.dat","w");
	fprintf(fileOut, "            t            kx            ky            kz         vbarx         vbary         vbarz\n");
	fprintf(fileOut,"% 13f % 13f % 13f % 13f % 13f % 13f % 13f\n", 0., Kx, Ky, Kz, vx(Kx, Ky, Kz), vy(Kx, Ky, Kz), vz(Kx, Ky, Kz) );

	double vbarx = 0;
	double vbary = 0;
	double vbarz = 0;

	double KxInit=Kx;
	double KyInit=Ky;
	for (double KzInit=0; KzInit<0.1; KzInit+=0.1){ 
		Kx = KxInit;	
		Ky = KyInit;
		Kz = KzInit;
		
		int nn=0; 
		while (exp(-nn*h/tau)>0.00001) {
			nn++;

			double k1x = +h*lorentzFx( Kx          , Ky         , Kz);
			double k1y = +h*lorentzFy( Kx          , Ky         , Kz);
			double k1z = +h*lorentzFz( Kx          , Ky         , Kz);
			
			// compute next point with Runge-Kutta
			double k2x = +h*lorentzFx( Kx + k1x/2. , Ky + k1y/2., Kz + k1z/2. );
			double k2y = +h*lorentzFy( Kx + k1x/2. , Ky + k1y/2., Kz + k1z/2. );
			double k2z = +h*lorentzFz( Kx + k1x/2. , Ky + k1y/2., Kz + k1z/2. );
			
			double k3x = +h*lorentzFx( Kx + k2x/2. , Ky + k2y/2., Kz + k2z/2. );
			double k3y = +h*lorentzFy( Kx + k2x/2. , Ky + k2y/2., Kz + k2z/2. );
			double k3z = +h*lorentzFz( Kx + k2x/2. , Ky + k2y/2., Kz + k2z/2. );
			
			double k4x = +h*lorentzFx( Kx + k3x    , Ky + k3y   , Kz + k3z    );
			double k4y = +h*lorentzFy( Kx + k3x    , Ky + k3y   , Kz + k3z    );
			double k4z = +h*lorentzFz( Kx + k3x    , Ky + k3y   , Kz + k3z    );
			
			double newKx = Kx + k1x/6. + k2x/3. + k3x/3. + k4x/6.;
			double newKy = Ky + k1y/6. + k2y/3. + k3y/3. + k4y/6.;
			double newKz = Kz + k1z/6. + k2z/3. + k3z/3. + k4z/6.;
			
			vbarx += vx(newKx, newKy, newKz)*exp(-nn*h/tau);
			vbary += vy(newKx, newKy, newKz)*exp(-nn*h/tau);
			vbarz += vz(newKx, newKy, newKz)*exp(-nn*h/tau);

			if (wrapBrillouin){
				if (newKx >  M_PI) {newKx = newKx - 2*M_PI; fprintf(fileOut,"\n");}
				if (newKx < -M_PI) {newKx = newKx + 2*M_PI; fprintf(fileOut,"\n");}
				if (newKy >  M_PI) {newKy = newKy - 2*M_PI; fprintf(fileOut,"\n");}
				if (newKy < -M_PI) {newKy = newKy + 2*M_PI; fprintf(fileOut,"\n");}
				if (newKz >  2*M_PI) {newKz = newKz - 4*M_PI; fprintf(fileOut,"\n");}
				if (newKz < -2*M_PI) {newKz = newKz + 4*M_PI; fprintf(fileOut,"\n");}
			}

			fprintf(fileOut,"% 13f % 13f % 13f % 13f % 13f % 13f % 13f \n", nn*h, newKx, newKy, newKz, vbarx, vbary, vbarz);
			Kx = newKx;
			Ky = newKy;
			Kz = newKz;

			printf("t = % 4f / % 4f , vbar = (% 4f,% 4f,% 4f)\n", nn*h, Nt*h, vbarx, vbary, vbarz);//fflush(stdout);
		}
	fprintf(fileOut,"\n\n");
	}
}









int res = 100;
double step = 1/100.;
int visited[100][100][100];
double treshold = 1.4;
int iter=0;

double spectral(int ii, int jj, int kk){
	double Kx = 2* M_PI* ii* step;
	double Ky = 2* M_PI* jj* step;
	double Kz = 2* M_PI* kk* step;
	
	return 0.05/(xi(Kx,Ky,Kz)*xi(Kx,Ky,Kz)+0.025);
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
	if (visited[ii][jj][kk]) return;
	
	iter++;
	visited[ii][jj][kk] = 1;

	if (spectral(ii,jj,kk) > treshold)
	{
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

	FILE *fileOut = fopen("mdc.out","w");

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
	long int total=0;
	for (int ii = 0; ii<10000; ii++){
		for (int jj = 0; jj<10000; jj++){
			for (int kk = 0; kk<10000; kk++){
				total++;
			}
		}
	}
	printf("%li",total);
	return 0;
}
