//
//  main.cpp
//  transportBuddy for 3D cuprates
//  by Simon Verret, Universite de Sherbrooke
//
//  code partially adapted from afmCond.c
//  by Maxime Charlebois, Universite de Sherbrooke
//  under MIT license (september 2017)
//

#define _USE_MATH_DEFINES
#include <math.h>
#include <complex.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



// FUNCTIONS FOR READING FILES

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




//// MAIN

int main(int argc, const char * argv[]) {
    
    
    //// SETTING PARAMETERS
    
    
    int markiewicz = 1; //bool
    double t;
    double tp;
    double tpp;
    double tppp;
    double tz;
    
    double ETA;
    double ETAell;
    double ETAdos;
    double ETAaFL;
    double ETAbFL;
    double ETAaPL;
    int nOmega;
    int nK;
    int nKz;
    int nMu; double muMin; double muMax;
    int nT; double Tmin; double Tmax;
    int logT=0; //bool
    double amplitudeCutoff;
    
    
    //// SETTING PARAMETERS
    FILE * file = fopen("model.dat", "rt");
    if(file == NULL) {printf("file %s not found", "model.dat"); exit(1);}
    printf("reading parameters from model.dat\n\n") ;
    
    readDouble(file, "t",               &t);
    readDouble(file, "tp",              &tp);
    readDouble(file, "tpp",             &tpp);
    readDouble(file, "tppp",            &tppp);
    readDouble(file, "tz",              &tz);
    
    readDouble(file, "muMin",           &muMin);
    readDouble(file, "muMax",           &muMax);
    readDouble(file, "Tmin",            &Tmin);
    readDouble(file, "Tmax",            &Tmax);
    readDouble(file, "ETA",             &ETA);
    readDouble(file, "ETAell",          &ETAell);
    readDouble(file, "ETAdos",          &ETAdos);
    readDouble(file, "ETAaFL",          &ETAaFL);
    readDouble(file, "ETAbFL",          &ETAbFL);
    readDouble(file, "ETAaPL",          &ETAaPL);
    readDouble(file, "amplitudeCutoff", &amplitudeCutoff);
    
    readInt(file, "nMu",    &nMu);
    readInt(file, "logT",   &logT); //to act as a bool
    readInt(file, "nT",     &nT);
    readInt(file, "nK",     &nK);
    readInt(file, "nKz",    &nKz);
    readInt(file, "nOmega", &nOmega);
    
    fclose(file);
    
    

    
    //// COMPUTING
    
    printf("transportBuddy starting\n\n");
    FILE *fileOut = fopen("transport_vs_T.dat","w");
    
    
    //// PRE-CALCULATE TEMPERATURE VARIABLES
    double T[nT]; double beta[nT]; double omega[nT][nOmega];
    double energyCutoff[nT]; double dfermiDirac_dw[nT][nOmega]; double dfermiDirac_dT[nT][nOmega];
    
    int nn=0; for(nn=0; nn<nT; nn++) {
        T[nn] = Tmin;
        if ((Tmin != Tmax) && (nT != 1)) {
            if(logT==1) T[nn] *= exp(nn*log(Tmax/Tmin)/(nT-1));
            else T[nn] += nn*(Tmax-Tmin)/(nT-1);
        }
        beta[nn] = 1./T[nn];

        energyCutoff[nn] = 2.*2.*acosh(0.25*sqrt(beta[nn]/amplitudeCutoff)) /beta[nn];
        int n=0; for(n=0; n<nOmega; n++)
        {
            omega[nn][n]= -energyCutoff[nn] + 2.*energyCutoff[nn]*n/(nOmega-1);
            double expBw = exp(beta[nn]*omega[nn][n]);
            dfermiDirac_dw[nn][n]= -beta[nn]*expBw/((expBw+1.)*(expBw+1.));
            dfermiDirac_dT[nn][n]= -beta[nn]*omega[nn][n] * dfermiDirac_dw[nn][n];
        }
    }
    
    //// PRE-CALCULATE SINES AND COSINES
    double sink[nK]; double sin2k[nK];  double sink_2[nK];
    double cosk[nK]; double cos2k[nK];  double cosk_2[nK];
    int ii=0; for(ii=0; ii<nK; ii++){
        double k = M_PI*(ii*1.0/nK);
        sink[ii] = sin(k); sin2k[ii] = sin(2.*k); sink_2[ii] = sin(k/2.);
        cosk[ii] = cos(k); cos2k[ii] = cos(2.*k); cosk_2[ii] = cos(k/2.);
    }
    double coskz[nKz]; double coskz_2[nKz]; double sinkz_2[nKz];
    int kk=0; for(kk=0; kk<nKz; kk++){
        double kz = M_PI*(kk*2.0/(nKz)); // 0 to 2Pi assumes parity (period of Markie is 4Pi)
        coskz[kk] = cos(kz); coskz_2[kk] = cos(kz/2.); sinkz_2[kk] = sin(kz/2.);
    }
    
    
    //// RESULTS RECORDER
    double dos[nMu];
    double density0[nMu];
    double density[nMu][nT];
    
    //// LOOP ON mu
    for(int iMu=0; iMu<nMu; iMu++)
    {
        double mu = muMin;
        if (muMin!=muMax) mu += iMu*(muMax-muMin)/(nMu-1);
        
        //// RESULTS
        dos[iMu] = 0.;
        density0[iMu] = 0.;
        int nn=0; for(nn=0; nn<nT; nn++)
        {
            density[iMu][nn]=0.;
        }
        
        //// INTEGRAL ON k
        for(int ii=0; ii<nK; ii++) {
            printf("muu= %i / %i, k = %i / %i \n",iMu,nMu,ii,nK);//fflush(stdout);
            
            for(int jj=0; jj<nK; jj++) {
                
                double epsilon_k   = -mu;
                epsilon_k         += -2.*t    * (cosk[ii] + cosk[jj]);
                epsilon_k         += -4.*tp   *  cosk[ii]*cosk[jj];
                epsilon_k         += -2.*tpp  * (cos2k[ii] + cos2k[jj]);
                epsilon_k         += -4.*tppp * (cos2k[ii]*cosk[jj] + cosk[ii]*cos2k[jj]);
                double depsilon_dkx = 0;
                depsilon_dkx      +=  2.*t    *  sink[ii];
                depsilon_dkx      +=  4.*tp   *  sink[ii]*cosk[jj];
                depsilon_dkx      +=  4.*tpp  *  sin2k[ii];
                depsilon_dkx      +=     tppp * (8.*sin2k[ii]*cosk[jj] + 4.*sink[ii]*cos2k[jj]);
                double depsilon_dky = 0;
                depsilon_dky      +=  2.*t    *  sink[jj];
                depsilon_dky      +=  4.*tp   *  sink[jj]*cosk[ii];
                depsilon_dky      +=  4.*tpp  *  sin2k[jj];
                depsilon_dky      +=     tppp * (8.*sin2k[jj]*cosk[ii] + 4.*sink[jj]*cos2k[ii]);
                double depsilon_dkx_dky = 0;
                depsilon_dkx_dky  += -4.*tp   *  sink[ii]*sink[jj];
                depsilon_dkx_dky  += -8.*tppp * (sin2k[ii]*sink[jj] + sink[ii]*sin2k[jj]);
                double depsilon_dky_dky = 0;
                depsilon_dky_dky  +=  2.*t    *  cosk[jj];
                depsilon_dky_dky  +=  4.*tp   *  cosk[jj]*cosk[ii];
                depsilon_dky_dky  +=  8.*tpp  *  cos2k[jj];
                depsilon_dky_dky  +=     tppp * (16.*cos2k[jj]*cosk[ii] + 8.*cosk[jj]*cos2k[ii]);
                
                
                for(int kk=0; kk<nKz; kk++) {
                    
                    double epsilonz_k         = -2.*tz*coskz[kk];
                    double depsilonz_dkx      =  0.;
                    double depsilonz_dky      =  0.;
                    double depsilonz_dkx_dky  =  0.;
                    double depsilonz_dky_dky  =  0.;
                    double depsilonz_dkz      =  2.*tz*sink[ii];
                    if (markiewicz) {
                        epsilonz_k             = -2.*tz*coskz_2[kk];
                            epsilonz_k        *=  (cosk[ii]-cosk[jj])*(cosk[ii]-cosk[jj]);
                            epsilonz_k        *=   cosk_2[ii]*cosk_2[jj];
                        depsilonz_dkx          = -2.*tz*coskz_2[kk];
                            depsilonz_dkx     *= -2*(cosk[ii]-cosk[jj])*sink[ii]*cosk_2[ii]*cosk_2[jj];
                            depsilonz_dkx     *=    (cosk[ii]-cosk[jj])*(cosk[ii]-cosk[jj])*sink_2[ii]*cosk_2[jj] /2.;
                        depsilonz_dky          = -2*(cosk[ii]-cosk[jj])*sink[jj]*cosk_2[ii]*cosk_2[jj];
                            depsilonz_dky     *=    (cosk[ii]-cosk[jj])*(cosk[ii]-cosk[jj])*cosk_2[ii]*sink_2[jj] /2.;
                        depsilonz_dkx_dky      = -2.*tz*coskz_2[kk];
                            depsilonz_dkx_dky *=-2*(cosk[ii]+sink[jj])*sink[ii]*cosk_2[ii]*cosk_2[jj];
                            depsilonz_dkx_dky *=  2*(cosk[ii]-cosk[jj])*sink[ii]*cosk_2[ii]*sink_2[jj] /2.;
                            depsilonz_dkx_dky *= -2*(cosk[ii]-cosk[jj])*sink[jj]*sink_2[ii]*cosk_2[jj] /2.;
                            depsilonz_dkx_dky *=   -(cosk[ii]-cosk[jj])*(cosk[ii]-cosk[jj])*sink_2[ii]*sink_2[jj] /4.;
                        depsilonz_dky_dky      = -2.*tz*coskz_2[kk];
                            depsilonz_dky_dky *= -2*(cosk[ii]+sink[jj])*sink[jj]*cosk_2[ii]*cosk_2[jj];
                            depsilonz_dky_dky *= -2*(cosk[ii]-cosk[jj])*cosk[jj]*cosk_2[ii]*cosk_2[jj];
                            depsilonz_dky_dky *=  2*(cosk[ii]-cosk[jj])*sink[jj]*cosk_2[ii]*sink_2[jj] /2.;
                            depsilonz_dky_dky *= -2*(cosk[ii]-cosk[jj])*sink[jj]*cosk_2[ii]*sink_2[jj] /2.;
                            depsilonz_dky_dky *=    (cosk[ii]-cosk[jj])*(cosk[ii]-cosk[jj])*cosk_2[ii]*cosk_2[jj] /4.;
                        depsilonz_dkz          = -2.*tz*sinkz_2[kk];
                            depsilonz_dkz     *=  (cosk[ii]-cosk[jj])*(cosk[ii]-cosk[jj]);
                            depsilonz_dkz     *=   cosk_2[ii]*cosk_2[jj];
                    }

                    double ep_k        =  epsilon_k   +  epsilonz_k;
                    double dep_dkx     = depsilon_dkx + depsilonz_dkx;
                    double dep_dky     = depsilon_dky + depsilonz_dky;
                    double dep_dkz     = depsilonz_dkz;
                    double dep_dkx_dky = depsilon_dkx_dky + depsilonz_dkx_dky;
                    double dep_dky_dky = depsilon_dky_dky + depsilonz_dky_dky;
                    

                    double Gamma = ETA;
                    double normV_k = sqrt(dep_dkx*dep_dkx + dep_dky*dep_dky + dep_dkz*dep_dkz);
                    
                    Gamma += ETAell*normV_k;
                    // Gamma += ETAdos/normV_k;
                    // double Gamma = ETA + //// Fermi Liquid 

                    //// A
                    double Ak0   = -(1./M_PI)*cimag(1.0/ (I*Gamma - ep_k));
                    double kernel_xx = dep_dkx*dep_dkx;
                    double kernel_zz = dep_dkz*dep_dkz;
                    double kernel_xy = -(2./3.)*(dep_dkx * (dep_dkx * dep_dky_dky - dep_dky * dep_dkx_dky));
                    
                    
                    dos[iMu]      += Ak0;
                    density0[iMu] += 1.0/(1.0+exp(1000*ep_k));
                    
                    
                    //// LOOP ON TEMPERATURES
                    
                    int nn=0; for(nn=0; nn<nT; nn++)
                    {
                        density[iMu][nn]    += 1.0/(1.0+exp(beta[nn]*ep_k));
                        Gamma += ETAbFL*T[nn]*T[nn];
                        //// INTEGRAL IN ENERGY
                        int n=0; for(n=0; n<nOmega; n++)
                        {
                            Gamma += ETAaFL*omega[nn][n]*omega[nn][n];
                            double complex z = omega[nn][n] + Gamma * I;
                            double A_k = -(1./M_PI)*cimag(1.0/ (z-ep_k) );
                        }
                    }
                }
            }
        }
        
        double f0 = 2.0/nK/nK/nKz;
        printf("% 4.8f % 4.8f % 4.8f ", mu, 1.0-f0*density0[iMu], f0*dos[iMu]);
        printf("% 4.8f % 4.8f % 4.8f ", f0*sigmaxx0[iMu], f0*sigmaxy0[iMu], f0*sigmazz0[iMu]);
        printf("\n");
        fflush(stdout);
        

        //// FILE GROUPED BY T
        fprintf(fileOut, "           mu            p0           dos ");
        fprintf(fileOut, "     averageV          minV          maxV ");
        fprintf(fileOut, "     sigmaxx0      sigmaxy0      sigmazz0 ");
        fprintf(fileOut, "            T             p            Cv ");
        fprintf(fileOut, "      sigmaxx       sigmaxy       sigmazz ");
        fprintf(fileOut, "      alphaxx       alphaxy       alphazz ");
        fprintf(fileOut, "       betaxx        betaxy        betazz\n");
        
        //// LOOP ON TEMPERATURES
        nn=0; for(nn=0; nn<nT; nn++)
        {
            double f = (2.*energyCutoff[nn]) /(nOmega)*f0;
            
            fprintf(fileOut,"% 13f % 13f % 13f ", mu, 1.0-f0*density0[iMu], f0*dos[iMu]);
            fprintf(fileOut,"% 13f % 13f % 13f ", aver_Vk0[iMu]/dos[iMu], min_Vk0[iMu], max_Vk0[iMu]);
            fprintf(fileOut,"% 13f % 13f % 13f ", f0*sigmaxx0[iMu], f0*sigmaxy0[iMu], f0*sigmazz0[iMu]);
            
            fprintf(fileOut,"% 13f % 13f % 13f ", T[nn], 1.0-f0*density[iMu][nn], f*Cv[iMu][nn]);
            fprintf(fileOut,"% 13f % 13f % 13f ", f*sigma_xx[iMu][nn], f*sigma_xy[iMu][nn], f*sigma_zz[iMu][nn] );
            fprintf(fileOut,"% 13f % 13f % 13f ", f*alpha_xx[iMu][nn], f*alpha_xy[iMu][nn], f*alpha_zz[iMu][nn] );
            fprintf(fileOut,"% 13f % 13f % 13f ", f*beta_xx[iMu][nn],  f*beta_xy[iMu][nn],  f*beta_zz[iMu][nn] );
            fprintf(fileOut, "\n");
        }
        fprintf(fileOut, "\n\n");
        fflush(fileOut);
        
        
    }
    fclose(fileOut);
    
    


    //// FILE GROUPED BY mu
    FILE *fileOutMu = fopen("transport_vs_mu.dat","w");
    
    nn=0; for(nn=0; nn<nT; nn++)
    {
        
        double f0 = 2.0/nK/nK/nKz;
        double f = (2.*energyCutoff[nn])/(nOmega)*f0;
        
        
        fprintf(fileOutMu, "           mu            p0           dos ");
        fprintf(fileOutMu, "     averageV          maxV          minV ");
        fprintf(fileOutMu, "     sigmaxx0      sigmaxy0      sigmazz0 ");
        fprintf(fileOutMu, "            T             p            Cv ");
        fprintf(fileOutMu, "      sigmaxx       sigmaxy       sigmazz ");
        fprintf(fileOutMu, "      alphaxx       alphaxy       alphazz ");
        fprintf(fileOutMu, "       betaxx        betaxy        betazz\n");
        
        for(int iMu=0; iMu<nMu; iMu++)
        {
            double mu = muMin + iMu*(muMax-muMin)/(nMu-1);
        
            fprintf(fileOutMu,"% 13f % 13f % 13f ", mu, 1.0-f0*density0[iMu], f0*dos[iMu]);
            fprintf(fileOutMu,"% 13f % 13f % 13f ", aver_Vk0[iMu]/dos[iMu], min_Vk0[iMu], max_Vk0[iMu]);
            fprintf(fileOutMu,"% 13f % 13f % 13f ", f0*sigmaxx0[iMu], f0*sigmaxy0[iMu], f0*sigmazz0[iMu]);
            fprintf(fileOutMu,"% 13f % 13f % 13f ", T[nn], 1.0-f0*density[iMu][nn], f*Cv[iMu][nn]);
            fprintf(fileOutMu,"% 13f % 13f % 13f ", f*sigma_xx[iMu][nn], f*sigma_xy[iMu][nn], f*sigma_zz[iMu][nn] );
            fprintf(fileOutMu,"% 13f % 13f % 13f ", f*alpha_xx[iMu][nn], f*alpha_xy[iMu][nn], f*alpha_zz[iMu][nn] );
            fprintf(fileOutMu,"% 13f % 13f % 13f ", f*beta_xx[iMu][nn],  f*beta_xy[iMu][nn],  f*beta_zz[iMu][nn] );
            fprintf(fileOutMu, "\n");
        
        }
        fprintf(fileOutMu, "\n\n");
    }
    
    printf("\n transportBuddy over.\n");
    return 0;
}





