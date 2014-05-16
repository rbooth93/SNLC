# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# define tau56ni 751680.0 //sec (8.7 d)
# define tau56co 9616320.0 //sec 111.3d
# define R0_norm 10.0
# define M56ni_norm 0.0
# define M_norm 3
# define E51 1.0
# define c 3e10
# define kappa 0.1
# define kappagamma 0.030
# define tend 100.

/* DESCRIPTION ----------------------------------
-------------------------------------------------
COMPUTES BOLOMETRIC LIGHT CURVE in 56CO POWERED
EJECTA FOLLOWING ARNETT1982 APPROXIMATIONS.
-------------------------------------------------
------------------------------------------------*/

int main(void) {
  
int it, jt, tgridpoints;
 double taum, tau0, rho, Rout, vout, taugamma, G, Dgamma, M56ni, M, E, tprim, dtprim, P, x, R0, Eth0;

M56ni = M56ni_norm*2e33; // ok
M = M_norm*2e33;    // ok
E = E51*1e51; // ok
R0 = R0_norm*6.96e10;  // to cm

tgridpoints = 300;  // Number of days to go to NOTE does NOT determine resolution of time integral

double *L = malloc(tgridpoints*sizeof(double));
double *Llarge = malloc(tgridpoints*sizeof(double));
double *time = malloc(tgridpoints*sizeof(double));

FILE *finout;
finout = fopen("out/arnett_2.txt", "w" );

vout = pow(10./3.*E/1.0/M, 0.5); // cm/s

printf("Vout : %g km/s\n", vout/1e5);

taum = 1.05/(pow(13.7*c,0.5))*pow(kappa,0.5)*pow(M,0.75)*pow(E,-0.25);  // Diffusion time Inserra 2013 Eq D2

 printf("taum %g\n",taum/(24*3600));

tau0 = kappa*M/(13.7*c*R0);   // ARnett 1982 Eq (22)

 printf("tau0 %g\n",tau0/(24*3600));

Eth0 = E/2.0; // Equipartition

printf("Eth0 %g\n",Eth0);

printf("Diffusion time Taum %g days\n", taum/(24*3600));

for (it = 0; it <= tgridpoints-1; it++) {  // Compute a value for each day...

  L[it] = 0.;  // Initialize
  Llarge[it] = 0;

  time[it] = ((double)it)/((double)tgridpoints)*tend*24*3600; //seconds

  tprim = 0;
  dtprim = 1000; // seconds
  while (tprim < time[it]) {

    P = 7.8e43*(M56ni/2e33)*exp(-tprim/tau56ni) + 1.4e43*(M56ni/2e33)*(exp(-tprim/tau56co) - exp(-tprim/tau56ni))/(1. - tau56ni/tau56co); // Inserra 2013 eq D4 and D5  erg/sec
    //P = 7.8e43*(M56ni/2e33)*exp(-tprim/tau56ni); // TEMP
   

    //printf("%g %g\n", time[it], P);

    // Compute deposition function Dgamma
    Rout = R0 + vout*tprim; // cm
    rho = M/(4*M_PI/3.*pow(Rout,3)); // g/cm3
    taugamma = kappagamma*rho*Rout; 
    G = taugamma/(taugamma + 1.6);    // Arnett 1982 Eq 51  
    Dgamma = G*(1 + 2*G*(1-G)*(1-0.75*G));    // Arnett 1982 eq 50 
 
    L[it] = L[it] + P*(0.97*Dgamma+0.03)*2*tprim/taum*exp( pow(tprim/taum, 2) - pow(time[it]/taum,2))*dtprim/taum; // Eq D1 in Inserra 2013

    Llarge[it] = Llarge[it] + P*(0.97*Dgamma+0.03)*(taum/tau0 + 2*tprim/taum)*exp( pow(tprim/taum, 2) - pow(time[it]/taum,2))*dtprim/taum;
                             // Arnett 1982 Eq  47- 48 NOT WORKING YET

    // printf("%g\t%g\n",

    tprim = tprim + dtprim;

    }

    Llarge[it] = Llarge[it] + exp(-(time[it]/tau0 + pow(time[it],2)/pow(taum,2)))*Eth0/tau0;  // Arnett 1982 Eq  47- 48 NOT WORKING YET

 }


for ( it = 0; it <= tgridpoints-1; it = it+1 )
  {
			
    fprintf(finout, "%g\t%g\t%g\n", time[it]/(24*3600), L[it], Llarge[it]);
   
  }
	
fclose( finout );

}



