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
# define max_theta 50 //Number of angles used for gamma ray deposition
# define max_ns 50 //Number of steps used for gamma ray deposition

// Fixed model parameters:
# define r0_rsun 1.0     // Initial radius in Rsun
# define kappa 0.1        // UVOIR opacity cm2/g
# define kappagamma 0.030 // gamma ray opacity cm2/g 
# define frac_kin 1.0   // Fraction of explosion energy used to compute velocities. Use 1 for compact progenitors, 0.5 for extended.
# define deplim 0.5    // Fraction of ejecta to do energy deposition in (takes from inner edge)

// Model parameter grid
# define MassNi_int_msun_min 0.07
# define MassNi_int_msun_max 0.07
# define MassNi_int_msun_step 0.01
# define Emin_E51 1.0
# define Emax_E51 1.0
# define Estep 0.1
# define Mej_msun_min 3.0
# define Mej_msun_max 3.0
# define Mej_msun_step 0.1

// END PARAMETERS :::::::::::::::::::::::::::::::::::::::::::::::::::

// Physical constants:
# define a 7.5657E-15	 // radiation constant
# define c 3E10		 // speed of light
# define tauNi 751680       // decay time scale for 56Ni in seconds
# define tauCo 9590400      // decay time scale for 56Co in seconds
# define pi 3.14159265359	//Pi

//FUNCTIONS
double calculate_Ldecay(double mass_Ni, double time);
double calculate_normalisation_factor();
double calculate_ds(double theta, double rinit, double r0, double radius_outer);
double calculate_velocity(int ns, double ds, double theta, double radius, double radius_outer, double v);
double calculate_dep(double rho, double ds);
int find_value_nm(double vvec[], double velocity);
double gamma_deposition(double Ldeposition[], double Ldecay[], double radius[], double rho_nm[], double vvec[], double radius_outer, double v, double r0_cm, double rinit[], double normalisation_factor);

int main(void){
    
  // INITIALISING EVERYTHING/////////////////////////////////////////////
  int tgridpoints=400000000; // some quantites are saved for each time step, this is maximum size of array to hold these values
  int nm, nm2, nt, ntfinal, row=99, col=2, i, j, ns, value_nm, ntheta; 
  double Tinitial, deltam, Mej_gram, tmax_sec, Eexplosion_erg, r0_cm, v, uinitial, deltat, tadiabatic, R, fluxlimiter, A1, A4, Amax;
  double dudt, Mej_msun, Eexplosion_E51, MassNi_int_msun, radius_outer, Ltot, taugamma, G, Dgamma, mass_cell;
  double sum_Ldep, sum_Ldecay, tau, normalisation_factor, tauwant;
  char name[256];    
  FILE *namefile;
  FILE *infile;
  FILE *outfile;
  outfile = fopen("tau.txt", "w");
  namefile = fopen("out/modellist.txt", "w");      // puts all output filenames into a list.
  infile = fopen("read_file.txt", "r"); //Opens file to read values
  double *rho = malloc(tgridpoints*sizeof(double)); // Scalar quantities for each time step
  double *tcourant = malloc(tgridpoints*sizeof(double)); 
  double *L = malloc(tgridpoints*sizeof(double));
  double *timesec = malloc(tgridpoints*sizeof(double));
  double *LNi = malloc(mgridpoints*sizeof(double));
  double *LCo = malloc(mgridpoints*sizeof(double));
  double *Ldecay = malloc(mgridpoints*sizeof(double));
  double *vvec = malloc(mgridpoints*sizeof(double)); // Mass-coordinate dependent quantities, overwrite each time step
  double *rinit = malloc(mgridpoints*sizeof(double));
  double *radius = malloc(mgridpoints*sizeof(double));
  double *utimesrhobis = malloc(mgridpoints*sizeof(double));
  double *utimesrhoprim = malloc(mgridpoints*sizeof(double));
  double *u = malloc(mgridpoints*sizeof(double));
  double *unext = malloc(mgridpoints*sizeof(double));
  double *MassNi_msun = malloc(mgridpoints*sizeof(double));
  double *mass_coord = malloc(mgridpoints*sizeof(double));
  double *rho_nm = malloc(mgridpoints*sizeof(double)); 
  double *Ldeposition = malloc(mgridpoints*sizeof(double));  
  double masscoord_MassNi[row][col];
  
  time_t t0, t1; // For timing purposes
  clock_t c0, c1;
  t0 = time(NULL);
  c0 = clock();
  
  for(i=0; i<row; i++){
    for(j=0; j<col; j++){
      fscanf(infile, " %lf", &masscoord_MassNi[i][j]); //Reading in values for vvec and corresponding Ni mass into a 2D array
    }
  } 
  
  normalisation_factor = calculate_normalisation_factor();
  //printf("%g\n", normalisation_factor); //DC prints normalisation factor
    
   for( MassNi_int_msun = MassNi_int_msun_min; MassNi_int_msun <= MassNi_int_msun_max; MassNi_int_msun = MassNi_int_msun + MassNi_int_msun_step){  // loop over initial 56Ni mass
     for(Eexplosion_E51 = Emin_E51; Eexplosion_E51 <= Emax_E51; Eexplosion_E51 = Eexplosion_E51 + Estep){      // loop over explosion energy
       for( Mej_msun = Mej_msun_min; Mej_msun <= Mej_msun_max; Mej_msun = Mej_msun + Mej_msun_step){     // loop over ejecta mass
    
	      FILE *finout;
        
	      snprintf( name, sizeof(name), "out/Ni_%g_E_%g_M_%g_R0_%g.txt", MassNi_int_msun, Eexplosion_E51, Mej_msun, r0_rsun); // Create filename with parameters baked in..gives for example Ni_0.07_E0_0.5_mass_1.5.txt
    
	      finout = fopen( name, "w+" );
    
	      Eexplosion_erg = Eexplosion_E51 * 1E51; // Convert to erg

	      Mej_gram = Mej_msun * 2E33;  // Convert to gram
	      
	      for(i=0; i<row; i++){
            masscoord_MassNi[i][0] = masscoord_MassNi[i][0] * Mej_gram; // Converting normalised mass coordinate into mass coordinate
          }
  
          for(nm=0; nm<mgridpoints-1; nm++){ // Loop over number of mass gridpoints
            mass_coord[nm] = ((float)nm/((float)mgridpoints-1)) * Mej_gram; //Calculating mass coordinate for a value of nm
	        for(i=0; i<row; i++){ // Search through array rows
              if(masscoord_MassNi[i][0] > mass_coord[nm]){
                MassNi_msun[nm] = masscoord_MassNi[i-1][1]; //If the two values match then allocate the corresponding 56 Ni mass
                break;
                }else{
                  continue;
              }
            }
          }       
           
	      mass_cell = Mej_gram / (mgridpoints-1); // Calculating the mass of one cell DC 
	      
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
	      
	      for (nm = 0; nm <= mgridpoints-1; nm++ ){ // Loop over mass grid
		
		  	rinit[nm] = r0_cm*pow( (float)nm/((float)(mgridpoints-1)), 1./3.); //  compute r as function of m, uniform sphere
		  	//printf("%d\t%g\n", nm, rinit[nm]); //DC rinit does vary with nm

		  	vvec[nm] = rinit[nm]/r0_cm*v; //  homology : make v proportional to r
		  	//printf("nm vvec %d\t%g\n", nm, vvec[nm]);
    
		  	unext[nm] = uinitial;
		  }
           
	      timesec[0] = 0;

	      nt = 0; // Time step index
    
	      while(timesec[nt] <= tmax_sec){ // Loop until final time is reached
		
		  	//u = unext;  This doesnt work well
		  	//u = &unext[0];
		  	for (nm = 0; nm <= mgridpoints-1; nm++){
		      u[nm] = unext[nm];  
		    }
          
          	radius_outer = r0_cm + v*timesec[nt]; // ok
          	//printf("aaa %g\t%g\t%g\n", v, timesec[nt], radius_outer);
             
          	rho[nt] = ( 3 * Mej_gram ) / ( 4 * M_PI * pow(radius_outer , 3 ) );
        
		  	Amax = 16*M_PI*M_PI*c/(3*kappa)*pow(radius_outer, 4);  // ok  Eq 2.5 in H13..typo an extra rho there
        
		  	tcourant[nt] = 0.1*pow(deltam,2)/(Amax*rho[nt]); //  FIX Time-step over which numerical instability arises  Eq 2.4 in H13
        
		  	tadiabatic = (r0_cm + v*timesec[nt])/v; // IMPORTANT Time scale over which density changes significantly - should not exceed this
        
		  	deltat = 0.25*tcourant[nt]; // Take a time-step a few times smaller than the courant step
        
		  	if (deltat > 0.25*tadiabatic) deltat = 0.25*tadiabatic; // ..but if adiabatic time scale is shorter use that.

		  	//if (timesec[nt] < 86400) deltat = 100.;  // First 24 h, do very short steps to resolve TEMP
		  
		  	for(nm=0; nm<mgridpoints-1; nm++){
		    	rho_nm[nm] = rho[nt]; //Setting the density for each mass co-ordinate as the density at that time.
		  	}
		  	//printf("ddd %g\t%g\t%g\t%g\n", rho_nm[0], r0_cm, radius_outer, rho_nm[1]);
		  
		    for (nm = 0; nm < mgridpoints-1; nm++){	
		      radius[nm] = rinit[nm] + vvec[nm]*timesec[nt];	//Calculating radius of each mass shell as time progresses.
		      Ldecay[nm] = calculate_Ldecay(MassNi_msun[nm], timesec[nt]); 	//Calculating Ldecay for each mass shell.
		    }
		    		    
		    tau = kappagamma * rho[nt] * radius_outer;
		    
		    /*tauwant = 10.0;
		    
			for(nm=0; nm<mgridpoints-1; nm++){ // TEMPORARY
		    	rho_nm[nm] = rho[nt]*tauwant/tau; //Setting the density for each mass co-ordinate as the density at that time.
		  	}*/ //Used for testing purposes
		    
		    if(tau>100){
		      for(nm = 0; nm < mgridpoints-1; nm++){
		      Ldeposition[nm] = Ldecay[nm]; //If optical depth is high then gamma rays are absorbed in the shell that they are produced
		      }
		      //printf("tau >100\t%g\n", timesec[nt]);
		    }else{  //If optical depth is low then gamma rays are able to travel
		      
		      gamma_deposition(Ldeposition, Ldecay, radius, rho_nm, vvec, radius_outer, v, r0_cm, rinit, normalisation_factor); //Calculates deposition in each mass shell
		      
		      //fprintf(outfile, "nm\tLdeposition\n");
		      /*for (nm = 0; nm < mgridpoints-1; nm++){
		         fprintf(outfile, "%d\t%g\n", nm, Ldeposition[nm]);
		      }
		      exit(0);*/
		    
		      sum_Ldep = 0;
		      sum_Ldecay = 0;
		    
		      for(nm=0; nm<mgridpoints-1; nm++){
		        sum_Ldecay = sum_Ldecay + Ldecay[nm];
		        sum_Ldep = sum_Ldep + Ldeposition[nm];
		      }
		      //printf("tau <100\t%g\n", timesec[nt]);
		      printf("Energy %g\t%g\n", sum_Ldep,  sum_Ldecay); //DC Energy conserved with both positrons and gamma rays
		    }
		    
		   for ( nm = 1; nm <= mgridpoints-2; nm ++ ){		// for each point on the mass grid (except two borders for which 2nd derivates cannot be done and will be replaced by boundary conditions)..
            
		      utimesrhoprim[nm] = rho[nt]*(u[nm+1] - u[nm])/deltam;  // forward derivative du/dx
            
		      utimesrhobis[nm] = rho[nt]*( u[nm+1] - 2*u[nm] + u[nm-1] ) / (pow(deltam, 2.0)); // centered second derivate d^2/dx^2
            
		      A1 = 64*M_PI*M_PI*c/3.0*pow(radius[nm], 1)/kappa/(4*M_PI*rho[nt]);  // Leading terms of term 1 (Eq 1.43) in H13..(a*dT4/dm is outside)
            
		      A4 = 16*M_PI*M_PI*c/3.0*pow(radius[nm], 4)/kappa; // ok  Leading terms of term 4 (Eq 1.46) in H13
            
		      R = 4*M_PI*pow(radius[nm], 2)/(kappa*u[nm])*fabs((u[nm] - u[nm-1])/deltam);  // Eq 1.19 in H13
            
		      fluxlimiter = (6 + 3*R)/(6 + 3*R + R*R); // fluxlimiter is "beta" in H13 writeup. eq 1.19
		          
		      for(nm=0; nm<mgridpoints-1; nm++){
		      		      
			  	dudt = (A1*utimesrhoprim[nm] + A4*utimesrhobis[nm])*fluxlimiter - u[nm]*(vvec[nm]/(rinit[nm]+vvec[nm]*timesec[nt])) + (Ldeposition[nm]/mass_cell); // add fluxlimiter and use rinit[nm]
            
		      	unext[nm] = u[nm] + dudt*deltat;  // Do the explicit time-derivate step
		                  
		      }
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

		  if ( nt >= tgridpoints ){
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

	    for ( nt = 0; nt <= ntfinal-1; nt = nt+20 ){
		  fprintf(finout, "%g\t%g\n", timesec[nt]/(86400.0), L[nt]);
		}

	    fclose( finout );
	    
	   }
     }   
   }

return EXIT_SUCCESS;
}

//FUNCTIONS
double calculate_Ldecay(double mass_Ni, double time){
  double LNi, LCo, Ldecay;
  
  LNi = 7.8E43 * mass_Ni * exp(-time/tauNi);  //56Ni decay luminosity..H13 eq 1.22

  LCo =  1.4E43 * mass_Ni * ((exp(-time/tauCo) - exp(-time/tauNi))/ ( 1 - (tauNi/tauCo)));  //56Co luminosity..H13 eq 1.23

  Ldecay = LNi + LCo;  // Total decay luminosity (56Ni + 56Co)
  
  return Ldecay;   
}

double calculate_normalisation_factor(){
  int ntheta, nm;
  double Li, sum_Li, theta, theta_rad, dtheta, normalisation_factor;
   
    theta = 0;
    
    sum_Li=0;
    
    for(ntheta=0; ntheta<max_theta; ntheta++){
      
      theta_rad = theta * (pi/180);
      
      dtheta = pi/(max_theta);
      
      Li = (sin(theta_rad)*(dtheta)*(1.0/2.0));
      
      sum_Li = sum_Li + Li;
      
      theta = theta + (180 / (double)(max_theta-1));
    }
    
    normalisation_factor = sum_Li;
  
  return normalisation_factor;
  
}

double calculate_ds(double theta, double rinit, double r0, double radius_outer){
  double theta_1, path_length, ds;
  
  theta_1 = 180 - theta;
		        
  theta_1 = theta_1 * (pi/180.0); //Convert to Radians
		        
  path_length = ((2*rinit*cos(theta_1))+pow((pow((2*rinit*cos(theta_1)),2)-(4*-(pow(r0,2)-pow(rinit,2)))),0.5))/2*(radius_outer/r0); // Length from current mass shell to surface in a particular direction
		        
  ds = path_length/max_ns; // Divide path length into segments
  
  return ds;
}

double calculate_velocity(int ns, double ds, double theta, double radius, double radius_outer, double v){
  double length, r2, theta_1, velocity;
  
  theta_1 = 180 - theta;
		        
  theta_1 = theta_1 * (pi/180.0); //Convert to Radians
  
  length = ns * ds; //Calculates how far along the path length
		          
  r2 = (pow((pow(radius,2)+pow(length,2)-(2*radius*length*cos(theta_1))),0.5)); //Calculates the radius of point that gamma ray is at from centre
		          
  velocity = v*(r2/radius_outer); //Calculate velocity at this radius
		          
  return velocity;
  
}

double calculate_dep(double rho, double ds){
  double tau, deposition;
  
  tau = kappagamma * rho * ds; // Calculating the optical depth
		         
  deposition = 1 - exp(-tau); //Calculating amount of gamma rays absorbed in segment
  
  return deposition;
}

int find_value_nm(double vvec[], double velocity){
  int value_nm, nm2;
  
  //printf("vvecmax velocity %g\t%g\n",vvec[
  
  for(nm2=0; nm2<=mgridpoints-1; nm2++){
    if((float)vvec[nm2] > (float)velocity ){ //Comparing vvec to velocity
      value_nm = nm2 - 1; //Allocate value of mass coordinate
      break;
	  }else{
	    if((float)velocity > (float)vvec[100]){
	      printf("Velocity is greater than vvec[100]\n");
	      value_nm = mgridpoints - 1;
	    }else{
	      continue;
	  }
	}
	
  }
  return value_nm;
}

double gamma_deposition(double Ldeposition[], double Ldecay[], double radius[], double rho_nm[], double vvec[], double radius_outer, double v, double r0_cm, double rinit[], double normalisation_factor){
  int nm, ntheta, ns, nm2, value_nm;
  double theta, Li, ds, velocity, deposition, theta_rad, dtheta, Li_init, L_escape=0;
  
  for (nm = 0; nm < mgridpoints-1; nm++){	
    Ldeposition[nm] = 0;	// Initialize to zero before each new calculation
  }
  
  for(nm=0; nm<mgridpoints-1; nm++){    	
	//Ldeposition[nm] = Ldeposition[nm] + (0.03 * Ldecay[nm]); //3% of energy goes to positrons which are absorbed in that shell
		    	
	theta=0; //Initialise theta
	
	if(Ldecay[nm] == 0){
	    continue;
	  }
		    		  
	for(ntheta=0; ntheta<max_theta; ntheta++){
	          
      theta_rad = theta * (pi/180.0); //Convert to Radians
      
      dtheta = pi / (max_theta);
	  
	  Li = ((1 * Ldecay[nm]) * (sin(theta_rad)*(dtheta)*(1.0/2.0))) / normalisation_factor; //97% of energy goes to gamma rays which are free to move
	
	  Li_init = Li;
			  
	  ds = calculate_ds(theta, rinit[nm], r0_cm, radius_outer); //Calculate the path length and divide into 100 segments - ds
		      	 	
	  for(ns=0; ns<max_ns; ns++){    	  
	    velocity = calculate_velocity(ns, ds, theta, radius[nm], radius_outer, v); //Calculate the velocity of material at that point
		                     
		value_nm = find_value_nm(vvec, velocity); //Compare velocity to vvec to find what mass shell we are currently in
		
		//printf("value_nm %d\n",value_nm);
		         
		deposition = calculate_dep(rho_nm[value_nm], ds); //Calculate deposition
		         
		Ldeposition[value_nm] = Ldeposition[value_nm] + (deposition*Li); //Adding amount of gamma ray absorbed to the total amount absorbed in that mass coordinate      		         
		         
		Li = Li - (deposition * Li); //Energy of gamma beam is decreased by the amount absorbed in segment
		
		//printf("%d\t%d\t%d\n", nm, ntheta, ns);
		            
		if(Li < (0.0001*Li_init)){
		  //printf("Li < 0.0001 initial\n");
		  break;
		}else{
		  continue;
		}
		               
	  }
	  L_escape = L_escape + Li;
	        
	  theta = theta + (180.0 / (double)(max_theta-1)); //Increase theta
		                   
	}  
  }
  //printf("Lescape %g\n", L_escape);
  return *Ldeposition;
}
