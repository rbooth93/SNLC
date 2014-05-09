# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <sys/types.h>
# include <unistd.h>

/* DESCRIPTION : 
Computes bolometric light curve of homologously expanding sphere.
Parameters : Mass, energy, 56NI mass.
Grid of parameters is defined in here, output is written as files
with parameters encoded into filename.
Authors : J. Haughey, A. Jerkstrand
*/

// PARAMETERS :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// Grid:
# define TMAX_days 50.    // Final light curve time in days.
# define mgridpoints 101  // The number of mass grid points to use

// Fixed model parameters:
# define r0_rsun 0.01     // Initial radius in Rsun
# define kappa 0.1        // UVOIR opacity cm2/g
# define kappagamma 0.030 // gamma ray opacity cm2/g 
# define frac_kin 1.0   // Fraction of explosion energy used to compute velocities. Use 1 for compact progenitors, 0.5 for extended.
# define deplim 0.5    // Fraction of ejecta to do energy deposition in (takes from inner edge)

// Model parameter grid
# define MassNi_min 0.07
# define MassNi_max 0.08
# define MassNi_step 0.01
# define Emin 0.5
# define Emax 0.5
# define Estep 0.1
# define Mej_min 1.5
# define Mej_max 1.5
# define Mej_step 0.1

// END PARAMETERS :::::::::::::::::::::::::::::::::::::::::::::::::::

// Physical constants:
# define a 7.5657E-15	 // radiation constant
# define c 3E10		 // speed of light
# define tauNi 751680       // decay time scale for 56Ni in seconds
# define tauCo 9590400      // decay time scale for 56Co in seconds

int main(void) {
    
  // INITIALISING EVERYTHING/////////////////////////////////////////////
  int tgridpoints=4000000; // some quantites are saved at each time step, this is size of array to hold these valyues
  int nm, nt, i, ntfinal; // loop indeces
  double Tinitial, deltam, MMAX, TMAX, Eexplosion, r0, v, uinitial, deltat, tadiabatic, R, fluxlimiter, radius, A1, A4, Amax, dudt, MMAX_norm, Eexplosion_norm, MassNi, radius_outer, Ltot, taugamma, G, Dgamma;   
  char name[256];    
  FILE *namefile;
  namefile = fopen("out/modellist.txt", "w");      // puts all filenames into separate file so leastsquares.c can analyse data
  double *rho = malloc(tgridpoints*sizeof(double));
  double *tcourant = malloc(tgridpoints*sizeof(double));
  double *L = malloc(tgridpoints*sizeof(double));
  double *timesec = malloc(tgridpoints*sizeof(double));
  double *vvec = malloc(mgridpoints*sizeof(double));
  double *rinit = malloc(mgridpoints*sizeof(double));
  double *LNi = malloc(tgridpoints*sizeof(double));
  double *LCo = malloc(tgridpoints*sizeof(double));
  double *Ldecay = malloc(tgridpoints*sizeof(double));
  double *utimesrhobis = malloc(mgridpoints*sizeof(double));
  double *utimesrhoprim = malloc(mgridpoints*sizeof(double));
  double *u = malloc(mgridpoints*sizeof(double));
  double *unext = malloc(mgridpoints*sizeof(double));
  // double x;

  time_t t0, t1; // For timing purposes
  clock_t c0, c1;
  t0 = time(NULL);
  c0 = clock();
    
  for( MassNi = MassNi_min; MassNi <= MassNi_max; MassNi = MassNi + MassNi_step)     // loop over Nickel mass
    {
    
      for(Eexplosion_norm = Emin; Eexplosion_norm <= Emax; Eexplosion_norm = Eexplosion_norm + Estep)      // loop over explosion energy
	{
    
	  for( MMAX_norm = Mej_min; MMAX_norm <= Mej_max; MMAX_norm = MMAX_norm + Mej_step)     // loop over ejecta mass
	    {
    
	      FILE *finout;
        
	      snprintf( name, sizeof(name), "Ni_%g_E_%g_M_%g.txt", MassNi, Eexplosion_norm, MMAX_norm); // Create filename with parameters baked in..gives for example Ni_0.07_E0_0.5_mass_1.5.txt
    
	      finout = fopen( name, "w+" );
    
	      Eexplosion = Eexplosion_norm * 1E51; // Convert to erg
	      MMAX = MMAX_norm * 2E33;  // Convert to gram
	      r0 = r0_rsun * 7E10;  // Convert to cm
	      TMAX = TMAX_days * 86400; // Convert to seconds

	      Ltot = 0.0; // Initialize
    
	      deltam = MMAX/((float)(mgridpoints-1));  // Mass step
       
	      printf ("mass = %g\t\tEexplosion = %g\t\tM(56Ni) = %g\nmgridpoints =  %d\n",MMAX_norm, Eexplosion_norm, MassNi, mgridpoints);
       
	      // INITIAL CONDITIONS///////// full description is missing in H13
    	      
	      Tinitial =  pow(( ( 3 * (Eexplosion/2.0) )/ ( 4 * M_PI * a * pow( r0, 3 ) ) ),0.25); // ok  assume equipartition, half explosion energy goes to internal energy
    
	      rho[0] = ( 3 * MMAX ) / ( 4 * M_PI * pow( r0, 3 ) ); // ok  uniform sphere expression
    
	      uinitial = a*pow(Tinitial,4)/rho[0]; //ok  assume radiation energy dominated
    
	      v = pow( (10 * (frac_kin*Eexplosion) )/ ( 3 * MMAX ), 0.5 );	//ok	
	      
	      for (nm = 0; nm <= mgridpoints-1; nm++ ) // Loop over mass grid
		{
		  rinit[nm] = r0*pow( (float)nm/((float)(mgridpoints-1)), 1.0/3.0); //  homology : make v proportional to r (not m)

		  vvec[nm] = rinit[nm]/r0*v; //
    
		  unext[nm] = uinitial;	   //SPECIFY INITIAL CONDITIONS HERE
		}
           
	      timesec[0] = 0;

	      nt = 0; // Time step index
    
	      while(timesec[nt] <= TMAX) // Loop until final time is reached
		{
        
		  //u = unext;  This doesnt work well
		  //u = &unext[0];
		  for (nm = 0; nm <= mgridpoints-1; nm++)
		    {
		      u[nm] = unext[nm];  
		    }

		  //printf(" %g\n", u[mgridpoints-2]-u[mgridpoints-3]);
		          
		  radius_outer = r0 + v*timesec[nt]; // ok
        
		  rho[nt] = ( 3 * MMAX ) / ( 4 * M_PI * pow(radius_outer , 3 ) ); // ok
        
		  Amax = 16*M_PI*M_PI*c/(3*kappa)*pow(radius_outer, 4);  // ok  Eq 2.5 in H13..typo an extra rho there
        
		  tcourant[nt] = 0.5*pow(deltam,2)/(Amax*rho[nt]); //  Time-step over which numerical instability arises  Eq 2.4 in H13
        
		  tadiabatic = (r0 + v*timesec[nt])/v; // IMPORTANT Time scale over which density changes significantly - cannot exceed this
        
		  deltat = 0.25*tcourant[nt]; // Take a time-step smaller than the courant step
        
		  if (deltat > 0.2*tadiabatic) deltat = 0.2*tadiabatic; // ..but if adiabatiuc time scale is shorter use that.
        
		  LNi[nt] = 7.8E43 * MassNi * exp(-timesec[nt]/tauNi);  //56Ni decay luminosity..H13 eq 1.22

		  LCo[nt] =  1.4E43 * MassNi * ((exp(-timesec[nt]/tauCo) - exp(-timesec[nt]/tauNi))/ ( 1 - (tauNi/tauCo)));  //56Co luminosity..H13 eq 1.23

		  Ldecay[nt] = LNi[nt] + LCo[nt];  // Total decay luminosity
        
		  for ( nm = 1; nm <= mgridpoints-2; nm ++ )		// for each point on the mass grid (except two borders for which 2nd derivates cannot be done)..
		    {
            
		      radius = rinit[nm] + vvec[nm]*timesec[nt]; // ok
            
		      utimesrhoprim[nm] = rho[nt]*(u[nm+1] - u[nm])/deltam;  // forward derivative
            
		      utimesrhobis[nm] = rho[nt]*( u[nm+1] - 2 * u[nm] + u[nm-1] ) / (pow(deltam, 2.0)); // centered second derivate
            
		      A1 = 64*M_PI*M_PI*c/3.0*pow(radius, 1)/kappa/(4*M_PI*rho[nt]);  // Leading terms in Eq 1.43 in H13..(a*dT4/dm is outside)
            
		      A4 = 16*M_PI*M_PI*c/3.0*pow(radius, 4)/kappa; // ok  Leading terms in Eq 1.46 in H13
            
		      R = 4*M_PI*pow( radius, 2)/(kappa*u[nm])*fabs((u[nm] - u[nm-1])/deltam);  // Eq 1.19 in H13
            
		      fluxlimiter = (6 + 3*R)/(6 + 3*R + R*R); // fluxlimiter is "beta" in H13 writeup. eq 1.19
            
		      // Gamma ray Deposition
		      taugamma = kappagamma*rho[nt]*radius_outer;  // Eq 1.26 in H13
            
		      G = taugamma/(taugamma + 1.6);    // Arnett 1982 Eq 51. Eq 1.25 in H13
            
		      Dgamma = G*(1 + 2*G*(1-G)*(1-0.75*G));    // Arnett 1982 eq 50. Eq 1.24 in H13
		               
		      if( (float)nm/(float)mgridpoints <= deplim)  // Is the mass coordinate inside limit for deposition?
                
			{ 
			  dudt = (A1*utimesrhoprim[nm] + A4*utimesrhobis[nm])*fluxlimiter - u[nm]*(vvec[nm]/(rinit[nm]+vvec[nm]*timesec[nt])) + (1./deplim)*(Ldecay[nt]/MMAX)*(0.97*Dgamma + 0.03); // add fluxlimiter and use rinit[nm]                
			}
		      
		      else
			{                
			  dudt = (A1*utimesrhoprim[nm] + A4*utimesrhobis[nm])*fluxlimiter- u[nm]*(vvec[nm]/(rinit[nm]+vvec[nm]*timesec[nt]));// NEW add fluxlimiter and use rinit[nm]
			}
            
		      unext[nm] = u[nm] + dudt*deltat;  // Do the explicit time-derivate step
            
		    }
        
		  // Boundary conditions
		  unext[0] = unext[1];		//Neuman inner boundary condition (symmetry)

		  unext[mgridpoints - 1] = unext[mgridpoints-2]*(1.0/(1.0 + deltam/(4*M_PI*radius_outer*radius_outer*4/(3.0*kappa*2.0)))); // H13 Eq 2.13 This outer BC makes a smooth L(m) function close to surface
        		  
		  R = 4*M_PI*pow( radius_outer, 2)/(kappa*u[mgridpoints-2])*fabs((u[mgridpoints-2] - u[mgridpoints-3])/deltam);  // H13 eq 1.19
        
		  fluxlimiter = (6 + 3*R)/(6 + 3*R + R*R); // Bersten Eq 3.6. H13 Eq 1.18
        
		  L[nt] = -16*M_PI*M_PI*c/3.0*pow(radius_outer, 4)/kappa*fluxlimiter*rho[nt]*(u[mgridpoints-2] - u[mgridpoints-3])/deltam; // Bersten Eq 3.5. H13 Eq 1.48
        
		  //L[nt]= 4 * M_PI * pow(radius_outer, 2) * c /4.0 * u[nt][mgridpoints-2] * rho[nt] *2; // sigma T^4 = ac/4 T^4 = c/4 u*rho. But then Teff^4 = 2*T^4. But using this formula does not seem to work well...the fl
        
		  Ltot = Ltot + L[nt]*deltat;
        
		  timesec[nt+1] = timesec[nt] + deltat;

		  nt++;

		  if ( nt >= tgridpoints )
		    {
		      printf("WARNING out of time grid points before reaching final time. Light curve cut.\n");
		      break;

		    }

		}
   

	      t1 = time(NULL);
	      c1 = clock();
	      printf("Wall clock Time %d\n", t1-t0);
	      printf("CPU time %f\n", (float) (c1-c0)/CLOCKS_PER_SEC); 

	      ntfinal = nt;
	      printf("Ltot = %g\nntfinal = %d\n:\n", Ltot, nt);
	      fprintf( namefile, "%s\n", name);     // prints all new text file names to file - important
	      
	      for ( nt = 0; nt <= ntfinal-1; nt = nt+20 )
		{
        
		  fprintf(finout, "%g\t%g\n", timesec[nt]/(86400.0), L[nt]);
        
		}
	
        
	      fclose( finout );
	
	    }
	}   
    }
    
    return EXIT_SUCCESS;
}








