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
# define tmax_days 50.0    // Final light curve time in days.
# define mgridpoints 101  // The number of mass grid points to use

// Fixed model parameters:
# define r0_rsun 1.0     // Initial radius in Rsun
# define kappa 0.1        // UVOIR opacity cm2/g
# define kappagamma 0.030 // gamma ray opacity cm2/g 
# define frac_kin 1.0   // Fraction of explosion energy used to compute velocities. Use 1 for compact progenitors, 0.5 for extended.
# define deplim 0.5    // Fraction of ejecta to do energy deposition in (takes from inner edge)

// Model parameter grid
# define MassNi_int_msun_min 0.2
# define MassNi_int_msun_max 0.2
# define MassNi_int_msun_step 0.01
# define Emin_E51 1.0
# define Emax_E51 1.0
# define Estep 0.1
# define Mej_msun_min 5.0
# define Mej_msun_max 5.0
# define Mej_msun_step 0.1

// END PARAMETERS :::::::::::::::::::::::::::::::::::::::::::::::::::

// Physical constants:
# define a 7.5657E-15	 // radiation constant
# define c 3E10		 // speed of light
# define tauNi 751680       // decay time scale for 56Ni in seconds
# define tauCo 9590400      // decay time scale for 56Co in seconds

int main(void) {
    
  // INITIALISING EVERYTHING/////////////////////////////////////////////
  int tgridpoints=400000000; // some quantites are saved for each time step, this is maximum size of array to hold these values
  int nm, nt, i, ntfinal; 
  double Tinitial, deltam, Mej_gram, tmax_sec, Eexplosion_erg, r0_cm, v, uinitial, deltat, tadiabatic, R, fluxlimiter, radius, A1, A4, Amax;
  double dudt, Mej_msun, Eexplosion_E51, MassNi_int_msun, radius_outer, Ltot, taugamma, G, Dgamma;   
  char name[256];    
  FILE *namefile;
  namefile = fopen("out/modellist.txt", "w");      // puts all output filenames into a list.
  double *rho = malloc(tgridpoints*sizeof(double)); // Scalar quantities for each time step
  double *tcourant = malloc(tgridpoints*sizeof(double));
  double *L = malloc(tgridpoints*sizeof(double));
  double *timesec = malloc(tgridpoints*sizeof(double));
  double *LNi = malloc(mgridpoints*sizeof(double));
  double *LCo = malloc(mgridpoints*sizeof(double));
  double *Ldecay = malloc(mgridpoints*sizeof(double));
  double *vvec = malloc(mgridpoints*sizeof(double)); // Mass-coordinate dependent quantities, overwrite each time step
  double *rinit = malloc(mgridpoints*sizeof(double));
  double *utimesrhobis = malloc(mgridpoints*sizeof(double));
  double *utimesrhoprim = malloc(mgridpoints*sizeof(double));
  double *u = malloc(mgridpoints*sizeof(double));
  double *unext = malloc(mgridpoints*sizeof(double));
  double *MassNi_msun = malloc(mgridpoints*sizeof(double));
  
  time_t t0, t1; // For timing purposes
  clock_t c0, c1;
  t0 = time(NULL);
  c0 = clock();
  
  
   for( MassNi_int_msun = MassNi_int_msun_min; MassNi_int_msun <= MassNi_int_msun_max; MassNi_int_msun = MassNi_int_msun + MassNi_int_msun_step)  // loop over initial 56Ni mass
    {
     for(Eexplosion_E51 = Emin_E51; Eexplosion_E51 <= Emax_E51; Eexplosion_E51 = Eexplosion_E51 + Estep)      // loop over explosion energy
	  {
       for( Mej_msun = Mej_msun_min; Mej_msun <= Mej_msun_max; Mej_msun = Mej_msun + Mej_msun_step)     // loop over ejecta mass
	    {
    
	      FILE *finout;
        
	      snprintf( name, sizeof(name), "out/Ni_%g_E_%g_M_%g_R0_%g.txt", MassNi_int_msun, Eexplosion_E51, Mej_msun, r0_rsun); // Create filename with parameters baked in..gives for example Ni_0.07_E0_0.5_mass_1.5.txt
    
	      finout = fopen( name, "w+" );
    
	      Eexplosion_erg = Eexplosion_E51 * 1E51; // Convert to erg

	      Mej_gram = Mej_msun * 2E33;  // Convert to gram

	      r0_cm = r0_rsun * 7E10;  // Convert to cm

	      tmax_sec = tmax_days * 86400; // Convert to seconds

	      Ltot = 0.; // Initialize
    
	      deltam = Mej_gram/((float)(mgridpoints-1));  // Mass step
       
	      printf ("mass = %g\t\tEexplosion_erg = %g\t\tM(56Ni) = %g\nmgridpoints =  %d\n",Mej_msun, Eexplosion_E51, MassNi_int_msun, mgridpoints);
       
	      // INITIAL CONDITIONS///////// note full description is missing in H13, assume uniform (radiation dominated) energy density distribution
    	      
	      Tinitial =  pow(( ( 3 * (Eexplosion_erg/2.0) )/ ( 4 * M_PI * a * pow( r0_cm, 3 ) ) ),0.25); // ok  assume equipartition, half explosion energy goes to internal energy
    
	      rho[0] = ( 3 * Mej_gram ) / ( 4 * M_PI * pow( r0_cm, 3 ) ); // ok  uniform sphere expression
    
	      uinitial = a*pow(Tinitial,4)/rho[0]; //ok  assume radiation energy dominated
    
	      v = pow( (10 * (frac_kin*Eexplosion_erg) )/ ( 3 * Mej_gram ), 0.5 );	//ok outer velocity of homologously expanding uniform sphere, use fraction frac_kin of explosion energy	

	      for (nm = 0; nm <= mgridpoints-1; nm++ ) // Loop over mass grid
		{
		  rinit[nm] = r0_cm*pow( (float)nm/((float)(mgridpoints-1)), 1./3.); //  compute r as function of m, uniform sphere

		  vvec[nm] = rinit[nm]/r0_cm*v; //  homology : make v proportional to r 
    
		  unext[nm] = uinitial;	   //
		}
           
	      timesec[0] = 0;

	      nt = 0; // Time step index
    
	      while(timesec[nt] <= tmax_sec) // Loop until final time is reached
		{
        
		  //u = unext;  This doesnt work well
		  //u = &unext[0];
		  for (nm = 0; nm <= mgridpoints-1; nm++)
		    {
		      u[nm] = unext[nm];  
		    }

		  radius_outer = r0_cm + v*timesec[nt]; // ok
        
		  rho[nt] = ( 3 * Mej_gram ) / ( 4 * M_PI * pow(radius_outer , 3 ) ); // ok
        
		  Amax = 16*M_PI*M_PI*c/(3*kappa)*pow(radius_outer, 4);  // ok  Eq 2.5 in H13..typo an extra rho there
        
		  tcourant[nt] = 0.1*pow(deltam,2)/(Amax*rho[nt]); //  Time-step over which numerical instability arises  Eq 2.4 in H13
        
		  tadiabatic = (r0_cm + v*timesec[nt])/v; // IMPORTANT Time scale over which density changes significantly - should not exceed this
        
		  deltat = 0.1*tcourant[nt]; // Take a time-step a few times smaller than the courant step
        
		  if (deltat > 0.25*tadiabatic) deltat = 0.25*tadiabatic; // ..but if adiabatic time scale is shorter use that.

		  //if (timesec[nt] < 86400) deltat = 100.;  // First 24 h, do very short steps to resolve TEMP
              
		  for ( nm = 1; nm <= mgridpoints-2; nm ++ )		// for each point on the mass grid (except two borders for which 2nd derivates cannot be done and will be replaced by boundary conditions)..
		    {
                      
              radius = rinit[nm] + vvec[nm]*timesec[nt]; // ok
            
		      utimesrhoprim[nm] = rho[nt]*(u[nm+1] - u[nm])/deltam;  // forward derivative du/dx
            
		      utimesrhobis[nm] = rho[nt]*( u[nm+1] - 2*u[nm] + u[nm-1] ) / (pow(deltam, 2.0)); // centered second derivate d^2/dx^2
            
		      A1 = 64*M_PI*M_PI*c/3.0*pow(radius, 1)/kappa/(4*M_PI*rho[nt]);  // Leading terms of term 1 (Eq 1.43) in H13..(a*dT4/dm is outside)
            
		      A4 = 16*M_PI*M_PI*c/3.0*pow(radius, 4)/kappa; // ok  Leading terms of term 4 (Eq 1.46) in H13
            
		      R = 4*M_PI*pow( radius, 2)/(kappa*u[nm])*fabs((u[nm] - u[nm-1])/deltam);  // Eq 1.19 in H13
            
		      fluxlimiter = (6 + 3*R)/(6 + 3*R + R*R); // fluxlimiter is "beta" in H13 writeup. eq 1.19
            
		      // Gamma ray Deposition
		      taugamma = kappagamma*rho[nt]*radius_outer;  // Eq 1.26 in H13
            
		      G = taugamma/(taugamma + 1.6);    // Arnett 1982 Eq 51. Eq 1.25 in H13
            
		      Dgamma = G*(1 + 2*G*(1-G)*(1-0.75*G));    // Arnett 1982 eq 50. Eq 1.24 in H13
		      
		      MassNi_msun[nm] = MassNi_int_msun * exp(-vvec[nm]/v);
		      
		      LNi[nm] = 7.8E43 * MassNi_msun[nm] * exp(-timesec[nt]/tauNi);  //56Ni decay luminosity..H13 eq 1.22

		      LCo[nm] =  1.4E43 * MassNi_msun[nm] * ((exp(-timesec[nt]/tauCo) - exp(-timesec[nt]/tauNi))/ ( 1 - (tauNi/tauCo)));  //56Co luminosity..H13 eq 1.23

		      Ldecay[nm] = LNi[nm] + LCo[nm];  // Total decay luminosity (56Ni + 56Co)

		      if( (float)nm/(float)mgridpoints <= deplim)  // Is the mass coordinate inside limit for deposition?
                
			{ 
			  dudt = (A1*utimesrhoprim[nm] + A4*utimesrhobis[nm])*fluxlimiter - u[nm]*(vvec[nm]/(rinit[nm]+vvec[nm]*timesec[nt])) + (1./deplim)*(Ldecay[nm]/Mej_gram)*(0.97*Dgamma + 0.03); // add fluxlimiter and use rinit[nm]                
			}

		      else  // no energy source term
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
		      printf("WARNING - Out of time grid points before reaching final time. Light curve cut. Increase tgridpoints to solve. \n");
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
