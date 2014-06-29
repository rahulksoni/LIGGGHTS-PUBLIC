/* ----------------------------------------------------------------------
LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
Transfer Simulations

www.liggghts.com | www.cfdem.com
Christoph Kloss, christoph.kloss@cfdem.com

LIGGGHTS is based on LAMMPS
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulatorrad_b
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov
insert
Copyright (2003) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.

See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/*************************For reading data from files*************************************/
#include "assert.h"
#include <iostream>
#include "stdio.h"
#include <fstream>
#include <limits>
#include <cmath>
#include <cstdlib>		//for accessing rand number generator function
#include <ctime>		//for accessing current time to create a seed value for the random number generator
/*****************************************************************************************/

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"										//added//
#include "atom_vec.h"									//added//
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "domain.h"
#include "random_park.h"								//added//
#include "memory.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_particledistribution_discrete.h"			//added//
#include "fix_template_multiplespheres.h"				//added//
#include "particleToInsert.h"							//added//

/***************************************************************/
#include "fix_breakparticle_force.h"
#include "breakage_coordinates.h"
/***************************************************************/

using namespace LAMMPS_NS;
using namespace FixConst; //needed for END_OF_STEP//

/* ---------------------Function to jump to a specific line of file------------------------- */
using namespace std;
std::fstream& GotoLine(std::fstream& myfile, unsigned int num)
{
    myfile.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i)
    {
        myfile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return myfile;
}	

/* --------------------------Function to write an array to a file---------------------------- */
std::fstream& writefile(std::fstream& file, int nparticles, int ntimestep, double *array, int array_length, int array_type)
{
	if(file.is_open()) 
		{
			file<<"\n";
			if(array_type == 1) 
				{
					file<<ntimestep;	
					file<<"\t \t";
					file<<nparticles;
					file<<"\t \t";
				}	
			for(int i = 0; i < array_length; i++)
				{
					file<<array[i];
					file<<"\t";
				}	
		} 
}	
/*---------------------------------------------------------------------------------------------*/

FixBreakparticleForce::FixBreakparticleForce(LAMMPS *lmp, int narg, char **arg) :
  FixInsert(lmp, narg, arg)   
	{
		if(screen) fprintf(screen ,"\n \n \n \n Entering function FixBreakparticleForce(). \n");
		if(logfile) fprintf(logfile ,"\n \n \n \n Entering function FixBreakparticleForce(). \n");
				// set defaults first, then parse args
		init_defaults();      	
				//function defined in this file itself. Does following
						
		bool hasargs = true;
		ECS_flag = 0; 			
				//narg is the no of total arguments and iarg is an element in that array//
		while(iarg < narg && hasargs)
			{
				hasargs = false;		
					//if both ifs below fails then it will have its role in terminating the while loop//
				if (strcmp(arg[iarg],"volumefraction") == 0)    
					//if comparison passes//Please strcmp returns 0 when strings are matched//
					{
						//if we are here then **arg's iargth argument is volumefraction. That means iargth+1 argument will be its numeric value//
						if (iarg+2 > narg) error->all(FLERR,"Illegal fix breakparticle/force command");   
							//that means volumefraction must be always the last argument in the **arg list//
						volumefraction = atof(arg[iarg+1]);     
							//atof converts numeric value to its float form//if enter item is character string then it fails and resturns 0//here it assigns volumefraction's numeric value to variable volumefraction//
						if(volumefraction < 0. || volumefraction >= 1.)		
							//checking if volumefraction has a legitimate value//
							error->all(FLERR,"Illegal fix breakparticle/force command, 0 < 'volumefraction' < 1 required");
						iarg += 2;		
							//increasing iarg value by +2//if we were in this loop then we were already working with the last argument. Now increasing this by +2 will turn next while condition into false and will stop further attempts//
						hasargs = true;
					}else if (strcmp(arg[iarg],"force_threshold") == 0)  
						//if comparison passes//Means ith argument is force_thresold//
					{
						if(screen) fprintf(screen ,"\n Instructed to break based on force_thresold. \n \n Extracting the value of force_thresold. \n");
						if(logfile) fprintf(logfile ,"\n Instructed to break based on force_thresold. \n \n Extracting the value of force_thresold. \n");
						if (iarg+2 > narg) error->all(FLERR,"Illegal fix breakparticle/force command");		
							//force_thresold can never be an argument before last//
						f_break = atof(arg[iarg+1]); 	
							//assigning value of force_thresold to variable f_break//
						if(f_break <= 0. ) error->all(FLERR,"Illegal fix breakparticle/force command, 'force_treshold' must be > 0"); 	
							//checking the legitimate value of force_thresold//
						iarg += 2;	
							//condition to declare that our work is done and while lopp must be terminated now//
						hasargs = true;
						if(screen) fprintf(screen ,"\n f_thresold extracted as %f Newton \n",f_break);
						if(logfile) fprintf(logfile ,"\n f_thresold extracted as %f Newton \n",f_break);
						//In all it means that the last argument could be either volumefraction or force_thresold. Whatever it is we have to assign its value to an appropriate variable//
					}	
					/**********Above like similar condition is created here for ECS_thresold***********/
					else if (strcmp(arg[iarg],"ECS_threshold") == 0)  
						//if comparison passes//Means ith argument is ECS_thresold//
					{
						if(screen) fprintf(screen ,"\n Instructed to break based on ECS. \n Extracting the value of ECS_thresold. \n");
						if(logfile) fprintf(logfile ,"\n Instructed to break based on ECS. \n Extracting the value of ECS_thresold. \n");
						
						if (iarg+2 > narg) error->all(FLERR,"Illegal fix breakparticle/force command");		
							//ECS_thresold can never be an argument before last//
						ECS_break = atof(arg[iarg+1]); 	
							//assigning value of ECS_thresold to variable ECS_break//
						if(ECS_break <= 0. ) error->all(FLERR,"Illegal fix breakparticle/force command, 'ECS_treshold' must be > 0"); 	
							//checking the legitimate value of force_thresold//
						iarg += 2;	
							//condition to declare that our work is done and while lopp must be terminated now//
						hasargs = true;
						
						ECS_flag = 1;
												
						if(screen) fprintf(screen ,"\n ECS_thresold extracted as %f kWh/ton and turned ECS flag to 1. \n \n",ECS_break);
						if(logfile) fprintf(logfile ,"\n ECS_thresold extracted as %f kWh/ton and turned ECS flag to 1. \n \n",ECS_break);
						//This flag will help in determining whether breakage is govern by force_tresold or ECS_tresold//
						//In all it means that the last argument could be either volumefraction or force_thresold or ECS_tresold. Whatever it is we have to assign its value to an appropriate variable//	
						/***************************************************************************************/					
					} else error->all(FLERR,"Illegal fix breakparticle/force command, unknown keyword");
			}
		
		if(ECS_flag == 1)
			{					
					double sieves_series[] = {0.000106, 0.000125, 0.00015, 0.00018, 0.000212, 0.00025, 0.0003, 0.000355, 0.000425, 0.0005, 0.0006, 0.00071, 0.00085, 0.001, 0.00118, 0.0017, 0.00236, 0.00335, 0.00475, 0.0056, 0.0067, 0.008, 0.0095, 0.0112, 0.0125, 0.0132, 0.016, 0.019, 0.0224, 0.025, 0.0265, 0.0315, 0.0375, 0.045, 0.05, 0.053, 0.063, 0.075, 0.09, 0.1};			
						//Sieve sizes obtained from http://mltest.com/PDF/astm_chart_wstyler.pdf
						
					sieves_series_length = 40;
										
					fstream parameters;
					
					parameters.open("parametersfile", std::ios::in);

						/*--------------------------------------------------------------------------------------------*/
						GotoLine(parameters,5);

						parameters>>JK_parameter_A;
						parameters>>JK_parameter_b;
						parameters>>WI;
						parameters>>C_GM;		//coefficient preciding the size
						parameters>>C_WI;		//coefficient preciding the work index
						parameters>>C_intercept;		//contact in the linear model
						/*--------------------------------------------------------------------------------------------*/					
						GotoLine(parameters,9);
						
						parameters>>min_parent_size_to_break;		//below this won't broke to save computational cost
						parameters>>min_daughter_size;			
						parameters>>max_feed_particle_size;				
						/*--------------------------------------------------------------------------------------------*/								
						GotoLine(parameters,13);
						
						parameters>>ECS_max;		//thresold energy above which particle breaks	
						parameters>>con_fac_joule;		//factors to be adapted for force to ECS conversion
						parameters>>con_fac_force;
						/*--------------------------------------------------------------------------------------------*/
						GotoLine(parameters,17);
						
						parameters>>x_min;		//region dimensions to regenerate broken particles bigger than the min_parent_size_to_break
						parameters>>x_max;
						parameters>>y_min;
						parameters>>y_max;
						parameters>>z_min;
						parameters>>z_max;
					//	parameters>>mill_axis;		//to be used to rotate the insertion box for adjustment
						parameters>>regenerate_in_bin_switch;		
						/*-------------------------------------------------------------------------------------------*/
						GotoLine(parameters,20);
						
						parameters>>massdis_nevery;			//frequency in steps for writing the cumulative and differential mass in the files
						/*-------------------------------------------------------------------------------------------*/
						GotoLine(parameters,24);
						
						parameters>>bin_x_min;		//region dimensions to regenerate broken particles lesser than grate size or min_parent_size_to_break
						parameters>>bin_x_max;
						parameters>>bin_y_min;
						parameters>>bin_y_max;
						parameters>>bin_z_min;
						parameters>>bin_z_max;
						parameters>>passing_truncation_switch;		//a boolean operator to decide whether smaller sizes than min_parent_size_to_break will be regenerated or not
						/*-------------------------------------------------------------------------------------------*/

					parameters.close();
					
					/*--------------------------------------------------------------------------------------------------*/
						
					/*********************************************************************************************/					
					conversion_factor = (con_fac_joule / con_fac_force) * 2.77777778 / 10000000.0;	// = 7.3491e-10	
							//Estimating EI on particle based on force comparison with figure 6 of Nikhil's paper//In KWh//
							//1757.08 N is force back-calculated for EI = 4.65 Joule from Nikhil's paper//
							//"E.T. Tuzcu, N. Dhawan and R.K. Rajamani. Coarse particle fracture with the ultrafast load cell. Minerals & Metallurgical Processing, 2011, Vol. 28, No. 4, pp. 176-186"//
							
					alpha_intercept = C_WI * WI + C_intercept;		//WI is a constant for the system therefore will give the constant as shown
										
					if(screen) fprintf(screen ,"JK_parameter_A = %f \t JK_parameter_b = %f \t WI = %f \t C_GM = %f \t C_WI = %f \t C_intercept = %f \t alpha_intercept = %f \n \n", JK_parameter_A,JK_parameter_b,WI,C_GM,C_WI,C_intercept,alpha_intercept);
					if(logfile) fprintf(logfile ,"JK_parameter_A = %f \t JK_parameter_b = %f \t WI = %f \t C_GM = %f \t C_WI = %f \t C_intercept = %f \t alpha_intercept = %f \n \n", JK_parameter_A,JK_parameter_b,WI,C_GM,C_WI,C_intercept,alpha_intercept);

					if(screen) fprintf(screen , "min_parent_size_to_break = %f \t min_daughter_size = %f \t max_feed_particle_size = %f \n \n ",min_parent_size_to_break,min_daughter_size,max_feed_particle_size);
					if(logfile) fprintf(logfile , "min_parent_size_to_break = %f \t min_daughter_size = %f \t max_feed_particle_size = %f \n \n ",min_parent_size_to_break,min_daughter_size,max_feed_particle_size);
						
					if(screen) fprintf(screen , "ECS_max = %f \t conversion_factor_energy = %f J \t conversion_factor_force = %f N \t conversion_factor = %.6g \n \n ",ECS_max,con_fac_joule,con_fac_force,conversion_factor);
					if(logfile) fprintf(logfile , "ECS_max = %f \t conversion_factor_energy = %f J \t conversion_factor_force = %f N \t conversion_factor = %.6g \n \n ",ECS_max,con_fac_joule,con_fac_force,conversion_factor);
					
					if(screen) fprintf(screen, "x_min = %f \t x_max = %f \t y_min = %f \t y_max = %f \t z_min = %f \t z_max = %f \n \n", x_min,x_max,y_min,y_max,z_min,z_max);
					if(logfile) fprintf(logfile, "x_min = %f \t x_max = %f \t y_min = %f \t y_max = %f \t z_min = %f \t z_max = %f \n \n", x_min,x_max,y_min,y_max,z_min,z_max);						
						
					while(sieves_series[index_lower] < min_daughter_size)			//counts number of sieves smaller than the min daughter size
						{
							index_lower++;			//Nos. lesser than min_daughter_size
						}
						
					while(max_feed_particle_size < sieves_series[sieves_series_length - 1 - index_upper])		//counts number of sieves bigger than the max particle size
						{
							index_upper++;			//Nos. bigger than max_feed_particle_size
						}
					//index_upper--;					//for safe side

					if(screen) fprintf(screen ,"\n sieves_series_length = %d, index_lower = %d, index_upper = %d \n",sieves_series_length,index_lower,index_upper);
					if(logfile) fprintf(logfile ,"\n sieves_series_length = %d, index_lower = %d, index_upper = %d \n",sieves_series_length,index_lower,index_upper);
												
					sieves = NULL;					
					if(sieves) memory->destroy(sieves);	
					if(sieves) memory->sfree(sieves);	
					if(sieves) delete[] sieves;		
					sieves = new double[(sieves_series_length - index_lower - index_upper)];		//an array to contain only selected sieves	

						//if(screen) fprintf(screen ,"\n Sieves \n");
						//if(logfile) fprintf(logfile ,"\n Sieves \n");
						//for(int j = 0; j < (sieves_series_length - index_lower - index_upper); j++)		//record sieves that are actually in use
						//	{
						//		sieves[j] = sieves_series[index_lower + j];
								//if(screen) fprintf(screen ,"%f \t", sieves[j]);
								//if(logfile) fprintf(logfile ,"%f \t", sieves[j]);
						//	}			
					
						//if(screen) fprintf(screen ,"\n \n \n");
						//if(logfile) fprintf(logfile ,"\n \n \n");
						/*********************************************************************************************/

					/*--------------------------------------------------------------------------------------------------*/
					fstream SievesSelected;   ///////
					
					SievesSelected.open("Sieves_selected", std::ios::out);	
					SievesSelected<<"Sieves selected \n";
						
					for(int j = 0; j < (sieves_series_length - index_lower - index_upper); j++)		//prints the selected sieves into the file
						{
							sieves[j] = sieves_series[index_lower + j];
							SievesSelected<<sieves[j];
							SievesSelected<<"\t";
						}
					SievesSelected<<"\n";
						
					SievesSelected.close();	
					/*----------------------------------------------------------------------------------------------------*/

					mass_distribution_truncated = NULL;
					if(mass_distribution_truncated) memory->destroy(mass_distribution_truncated);
					mass_distribution_truncated = new double[sieves_series_length - index_lower - index_upper];	
					for(int i = 0; i < (sieves_series_length - index_lower - index_upper); i++) 	mass_distribution_truncated[i] = 0.0; 	//creating the zeros array of the size of the size of selected size	

					/*----------------------Obatining previous values if restart was performed----------------------------*/
			/*		if(update->ntimestep > 1000)
						{
							fstream MassDisTruncated_RetriveFromRestart;
						
							MassDisTruncated_RetriveFromRestart.open("mass_distribution_truncated", std::ios::in);

								for(int i = 0; i < (sieves_series_length - index_lower - index_upper); i++) 	mass_distribution_truncated[i]<<MassDisTruncated_RetriveFromRestart;

							MassDisTruncated_RetriveFromRestart.close();
						}		*/
					/*----------------------------------------------------------------------------------------------------*/

					/*---------------------------Initiating values in truncated mass file --------------------------------*/
					/*	if(passing_truncation_first_time_switch == true)
						{
							fstream MassDisTruncated;
							
							MassDisTruncated.open("Truncated_Mass_Distribution", std::ios::out);	
													
							for(int j = 0; j < (sieves_series_length - index_lower - index_upper); j++)		
								{
										MassDisTruncated<<0.0;
										MassDisTruncated<<"\t";
								}
							MassDisTruncated<<"\n";
											
							MassDisTruncated.close();	

							passing_truncation_first_time_switch = false;		
						}	*/
					/*----------------------------------------------------------------------------------------------------*/
					
					fill_bounds(x_min, x_max, y_min, y_max, z_min, z_max, bin_x_min, bin_x_max, bin_y_min, bin_y_max, bin_z_min, bin_z_max);		//function defined in breakage_coordinates.h
					
					/*--------------------------------------------------------------------------------------------------------------------*/					
					assert(JK_parameter_A > 1.0 && JK_parameter_A < 100.0);
					assert(JK_parameter_b > 0.0 && JK_parameter_b < 10.0);
					assert(WI > 0.0 && WI < 100.0);
					assert(ECS_max > ECS_break && ECS_max < 10.0);
					assert(min_daughter_size > sieves_series[0] && min_daughter_size < min_parent_size_to_break && min_parent_size_to_break < max_feed_particle_size);
					assert(con_fac_joule > 0.0);
					assert(con_fac_force > 0.0);
					assert(x_max > x_min);
					assert(y_max > y_min);
					assert(z_max > z_min);
					//	assert(mill_axis == 1 || mill_axis == 2 || mill_axis == 3 || mill_axis == -1 || mill_axis == -2 || mill_axis == -3);	
					assert(bin_x_max > bin_x_min);
					assert(bin_y_max > bin_y_min);
					assert(bin_z_max > bin_z_min);
					assert(passing_truncation_switch == 0 || passing_truncation_switch == 1);
					/*--------------------------------------------------------------------------------------------------------------------*/
			}
		
		if(screen) fprintf(screen ,"\n \n ");
		if(logfile) fprintf(logfile ,"\n \n ");
		
		// no fixed insertion target (such as # particles) since insertion triggered by break-ups//
		ninsert_exists = 0;   
			//changed from default value 1 to 0//

		// turn off overlap check and turn off start stats since both not required//
		check_ol_flag = 0;	
			//changed from default value 1 to 0//
		print_stats_start_flag = 0;

		breakdata = NULL; 
		ECS = NULL;
		size = NULL;	
			//No breakdata yet, initiation of variable//
			
		x_sphere = NULL;
		r_sphere = NULL;
			
		fix_break = NULL;		
			//No atom properties yet for new particles, default initiation//

		nevery = 1;	
			//action after every 1 timestep//

		n_break = 0;	
			//number of particles. initiation//
		mass_break = 0.;	
			//mass to be break. initiation//

		if(maxrad >= 1.) error->all(FLERR,"Fix breakparticle/force: Particle distribution must be relative, max radius must be < 1");
		
		if(screen) fprintf(screen ,"\n Exiting function FixBreakparticleForce(). \n");
		if(logfile) fprintf(logfile ,"\n Exiting function FixBreakparticleForce(). \n");
	}	//this means new radius is given in relative terms//
	

/* ----------------------------------Constructors------------------------------------ */

void FixBreakparticleForce::post_create() 		 
	//Setting-up property/atom args & default values//
	{
//		if(screen) fprintf(screen ,"\n Entering function post_create(). \n");
//		if(logfile) fprintf(logfile ,"\n Entering function post_create(). \n");
		
		FixInsert::post_create();

		if(!fix_break)		
			//if fix_break is not defined then it is being defined//
			{		//search for fix property/atom command in manual//This loop is for defining that//
					//Syntax for this fix	(see at: http://nf.nci.org.au/facilities/software/LIGGGHTS/doc/fix_property.html)
					//fix fix-id group-id property/atom variablename style(possible attributes: scalar,vector,atomtype,matrix,atomtypepair) restartvalue(possible attributes: yes/no) comm_ghost_values(possible attributes: yes/no) comm_reverse_value(possible attributes: yes/no) defaultvalue(s)//
				char *breakvar_name = new char[10+strlen(id)];
				char* fixarg[9];
				fixarg[0] = new char[10+strlen(id)];	//fix-id//
				
				sprintf(fixarg[0],"break_%s",id);
				
				fixarg[1]="all";	//group-id//
				fixarg[2]="property/atom";		//action style//reserves properties per atom to be accessed later//
				fixarg[3] = new char[10+strlen(id)];	//variable name//
				sprintf(breakvar_name,"break_%s",id);		
				strcpy(fixarg[3],breakvar_name);	
				fixarg[4]="scalar"; 	//scalar style//reserves one value per atom//
				fixarg[5]="yes";    //restartvalue//
				fixarg[6]="yes";    //comm_ghost_value//
				fixarg[7]="no";    //comm_reverse_ghost//
				fixarg[8]="0.";		//defaultvalues//
									//it will be fixarg[break_id all property/atom break_id scalar yes yes no 0]//
							
				fix_break = modify->add_fix_property_atom(9,fixarg,style);		
					//class FixPropertyAtom* add_fix_property_atom(int narg,char **arg,const char *)//declared in modify.h.//
					//creating and initiatin fix_break//
					//See fix_property atoms for actions and comments//
					//basically grows memory appropriately and assign default values at appropriate places//
				
				delete []fixarg[0];  //delete break_id//
				delete []fixarg[3];	 //delete variable name break_id//
				delete []breakvar_name;		//delete breakvar_name break_id//
			}
//		if(screen) fprintf(screen ,"\n Exiting function post_create(). \n");
//		if(logfile) fprintf(logfile ,"\n Exiting function post_create(). \n");
	}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::pre_delete() 		
	//pre_delete declared in fix_breakparticle_force.h//
	{
//		if(screen) fprintf(screen ,"\n Entering function pre_delete(). \n");
//		if(logfile) fprintf(logfile ,"\n Entering function pre_delete(). \n");
		
		modify->delete_fix(fix_break->id);	
			//deleteing id of fix_break//
//		if(screen) fprintf(screen ,"\n Exiting function pre_delete(). \n");
//		if(logfile) fprintf(logfile ,"\n Exiting function pre_delete(). \n");
	}

/* ---------------------------------------------------------------------- */

FixBreakparticleForce::~FixBreakparticleForce()		//destructor//
	{
//		if(screen) fprintf(screen ,"\n Entering function ~FixBreakparticleForce(). \n");
//		if(logfile) fprintf(logfile ,"\n Entering function ~FixBreakparticleForce(). \n");
		
		if(breakdata) memory->sfree(breakdata);     
		if(ECS) memory->sfree(ECS);
		if(size) memory->sfree(size);
		if(x_sphere) memory->sfree(x_sphere);
		if(r_sphere) memory->sfree(r_sphere);
		if(sieves) memory->sfree(sieves);
		if(mass_distribution) memory->sfree(mass_distribution);
		if(mass_distribution_truncated) memory->sfree(mass_distribution_truncated);
		if(number_per_break) memory->sfree(number_per_break);
//		if(volume_sphere) memory->sfree(volume_sphere);
//		if(mass_sphere) memory->sfree(mass_sphere);
		
		//deleting any available breakadata//
//		if(screen) fprintf(screen ,"\n Exiting function ~FixBreakparticleForce(). \n");
//		if(logfile) fprintf(logfile ,"\n Exiting function ~FixBreakparticleForce(). \n");
	}

/* ---------------------------------------------------------------------- */
void FixBreakparticleForce::init_defaults()   //default constructor//
	{
//		if(screen) fprintf(screen ,"\n Entering function init_defaults(). \n");
//		if(logfile) fprintf(logfile ,"\n Entering function init_defaults(). \n");
		
		f_break = 0.0;	
				
		/**************Initiating value for ECS_break****************/
		ECS_flag = 0;
		ECS_break = 0.0;		
		
		n_fragments = 0;
		size_iparticle = 0.;
		ECS_iparticle = 0.;
		nspheres = 0;
		/************************************************************/
		
			//default for thresold force and thresold ECS//
		volumefraction = 0.5;	
		fix_fragments = NULL;  
			//fix_fragments is a class pointer of FixTemplateMultispheres
//		if(screen) fprintf(screen ,"\n Exiting function init_defaults(). \n");
//		if(logfile) fprintf(logfile ,"\n Exiting function init_defaults(). \n");	
			
	}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::init()
	{
//		if(screen) fprintf(screen ,"\n Enetering function init(). \n");
//		if(logfile) fprintf(logfile ,"\n Entering function init(). \n");
		
		FixInsert::init();              
			//init() is defined in fix_insert.cpp//
			//assigns appropriate data to keyword multisphere//
			//perfomrs some error checks and assigns data fix multisphere or creates fix multisphere with NULL data//
			
//		if(screen) fprintf(screen ,"\n Exiting function init(). \n");
//		if(logfile) fprintf(logfile ,"\n Exiting function init(). \n");
	
	}

/* ----------------------------------------------------------------------
   calculate ninsert, insert_every, ninsert_per, massinsert, flowrates etc
   also perform error checks
------------------------------------------------------------------------- */

void FixBreakparticleForce::calc_insertion_properties()
	//Getting args from commands in script & executes accordingly//
	{
//		if(screen) fprintf(screen ,"\n Entering function calc_insertion_properties(). \n");
//		if(logfile) fprintf(logfile ,"\n Entering function calc_insertion_properties(). \n");
		
		// error checks
		/****************Added error check for ECS_break value*****************/
		if(f_break == 0. && ECS_break == 0.)  	
			//False//f_break was assigned force_thresold in line 82 of this file//
			error->all(FLERR,"Illegal fix breakparticle/force command, you have to specify 'force_treshold' or 'ECS_tresold'");
		if(nflowrate > 0. || massflowrate > 0.)  
		 	//False//Values were turned 0 with void FixInsert::init_defaults() in this file//
			//cannot assign number flowrate or massflowrate//default values are 0//
			error->all(FLERR,"Illegal fix breakparticle/force command, specifying 'nflowrate' or 'massflowrate' is not allowed");
		if(ninsert > 0 || massinsert > 0.)		
			//False//values were turned 0 in this file//
			//connot assign number of particles or mass to be inserted//defaults are 0//
			error->all(FLERR,"Illegal fix breakparticle/force command, specifying 'nparticles' or 'mass' is not allowed");
		if(insert_every <= 0)   
			//it cannot be negative//
			error->all(FLERR,"Illegal fix breakparticle/force command, specifying 'every' must be > 0");

		// fix holding particle fragments
		if(fix_distribution->n_particletemplates() != 1)   
			//class FixParticledistributionDiscrete *fix_distribution; declared in fix_insert.h//
			//class is defined in fix_particledistribution_discrete.h//
			//n_particletemplates() returns ntemplates. ntemplates is no of templates.//
			//See at line no 81 in fix_particledistribution.h//
			//Particles with only one type of atom_type can be considered//
			//See in fix_particledistribution_discrete.cpp. ntemplate is 4th argument in breakage command example//
			//"fix		pdfragments all particledistribution/discrete 1.  1 fragments 1.", therefore it is 1//
			//First 1 tell the total no of particle profiles, second one tells first profile, third 1 is for profile weightage//
			
			//Getting no of templates//
			//Next step will be getting the style in which this template will be executed; which is multisphere//
			error->all(FLERR,"Illegal fix breakparticle/force command, fix of type particledistribution/discrete must hold exactly one template");
			
		if(strcmp(fix_distribution->particletemplates()[0]->style,"particletemplate/multiplespheres")) 
			//inline class FixTemplateSphere** particletemplates(). declared in fix_particledistribution_discrete.h//class defined in fix_templete_sphere.h//
			//particletemplates() returns templates. Therefore, inline class FixTemplateSphere** templates()//declaring a 2d array//
			//templates = new FixTemplateSphere*[ntemplates];. ntemplates is 1//Action in all rows and first column//
			//only multisphere style is allowed//
			
			//Previous step has got the no of template//
			//This step gets the style in which this will be executed//
			error->all(FLERR,"Illegal fix breakparticle/force command, fix of type particledistribution/discrete must hold exactly one template of type fix particletemplate/multiplespheres");
			
	
		fix_fragments = static_cast<FixTemplateMultiplespheres*>(fix_distribution->particletemplates()[0]); 	
			//forcing to cast according to class FixTemplateMultispheres
			//casts particletemplates()[0] according to class FixTemplateMultispheres//
			//then transfers it to FixTemplateMultispheres *fix_fragments//
			//FixTemplateMultispheres *fix_fragments//
			//assigning particletemplates to fix_fragments which is class pointer of FixTemplateMultisphere//
			//So now class FixTemplateMultispheres know that its has to execute as per multisphere style and number of templates is only 1//

			// get number of fragments from fix
		n_fragments = fix_fragments->number_spheres();   
						
		density_particle = fix_fragments->particle_density;
			//if(screen) fprintf(screen, "FBPF: density_particle = %f \n",density_particle);
			//if(logfile) fprintf(logfile, "FBPF: density_particle = %f \n",density_particle);
			//density_particle = class_multisphere->particle_density;
						
			////if(screen) fprintf(screen ,"\n n_fragments assigned as %d in function calc_insertion_properties. \n",n_fragments);
			////if(logfile) fprintf(logfile ,"\n n_fragments assigned as %d in function calc_insertion_properties. \n",n_fragments);
		
			//This is equivalent to FixTemplateMultispheres -> number_spheres //
			//Therefore takes the value of number of spheres from FixTEmplateMultispheres class//
			
			//fix_fragments is a class pointer of FixTemplateMultiplespheres in fix_breakparticle_force.h//
			//FixTempleteMultisphere is defined in fix_templete_multispheres.h//
			//In this file itself int number_spheres(); is declared//
			//function in fix_template_multispheres.cpp returns nspheres which is no of spheres//
			//nspheres in fix_template_multispheres.cpp is no of spheres 10 in command in breakage example 
			//"fix		fragments all particletemplate/multiplespheres 1 atom_type 1 density constant 2500 nspheres 10 ntry 1000000 & spheres file fragmentfile scale 1.0"// 
			//assigning number of spheres from class FixTempleteMultisphere as number of fragments//
			//Data for 10 spheres writes in fragment file with scale value 1.0. Means actaul values will be printing//
		
		// do not need ninsert_per
//		if(screen) fprintf(screen ,"\n Exiting function calc_insertion_properties(). \n \n");
//		if(logfile) fprintf(logfile ,"\n Exiting function calc_insertion_properties(). \n \n");
	}

/* ---------------------------------------------------------------------- */
//This is needed to set the end-of-step function into the pipeline,
//END_OF_STEP is defined in namespace FixConst, so I added it above !

int FixBreakparticleForce::setmask()
	{   
//		if(screen) fprintf(screen ,"\n Entering function setmask(). \n");
//		if(logfile) fprintf(logfile ,"\n Entering function setmask(). \n");
		
		int mask = FixInsert::setmask();
		mask |= END_OF_STEP;
		
//		if(screen) fprintf(screen ,"\n Exiting function setmask(). \n");
//		if(logfile) fprintf(logfile ,"\n Exiting function setmask(). \n");
		
		return mask;		
	}

/* ---------------------------------------------------------------------- */

inline int FixBreakparticleForce::is_nearby(int i)
	{
//		if(screen) fprintf(screen ,"\n Entering function is_nearby(). \n");
//		if(logfile) fprintf(logfile ,"\n Entering function is_nearby(). \n");
	
//		if(screen) fprintf(screen ,"\n Exiting function is_nearby(). \n");
//		if(logfile) fprintf(logfile ,"\n Exiting function is_nearby(). \n");
		
		// need not check overlap with existing particles since we
		// know space originally taken by deleted particles is free
		return 0;		
	}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::end_of_step()  	
	//assigning properties//
	{
		if(screen) fprintf(screen ,"\n Entering function end_of_step(). \n \n ");
		if(logfile) fprintf(logfile ,"\n Entering function end_of_step(). \n \n");
		
		int nlocal = atom->nlocal;	
			//nlocal is number of atoms in this processor//
		double *flag = fix_break->vector_atom;		
			//vector_atom contains the default values for atoms//
			//Deafult values are taken from FixBreakparticleForce::post_create() in this file//
			//See fix_propoerty_atom.cpp//
			//assigning vector_atom from class FixPropertyAtom//
			//action taken in fix_property_atom//
		double **x = atom->x;	
			//these are the coordinates of atoms//
		double **f = atom->f;
			//force values//
		double *radius = atom->radius;		 
			//particles radius assignments//
		
		/******************Added particle mass*******************/
		double *rmass = atom->rmass;		
			//particles mass assignment//	
		/********************************************************/			
		int *mask = atom->mask;
		/****************Adding features to handle breakage by ECS_tresold********************/
		double f_sqr,f_break_sqr, ECS_particle;
		double EI_particle;
			//f_break_sqr = f_break * f_break;    
			//converting allotted thresold force as force acting on the particle//
			// checking breakage criterion for all local particles
		int flag_count = 0;	
		
		rmass_break = 0.0;
		
		for(int i = 0; i < nlocal; i++)			//checks the particles which are going to be breakon//
			{
				if (mask[i] & groupbit) 	
					//groupbit is the check whether particle belongs to current group or not//mask is the group number to which particle belongs//
					{
						if(ECS_flag == 1) 
							{			
								EI_particle = conversion_factor * sqrt(vectorMag3DSquared(f[i]));
				
								ECS_particle = EI_particle / rmass[i] * 1000; 
									//Calulating ECS acting on the particle//In kWh/ton//
														
									////		if(screen) fprintf(screen ,"end_of_step(): i = %d, f = %f N, radius = %f m, rmass = %f kg, EI_partcile = %f kWh, ECS_particle = %f kWh/ton \n",i,sqrt(vectorMag3DSquared(f[i])),radius[i],rmass[i],EI_particle,ECS_particle);
									////		if(logfile) fprintf(logfile ,"end_of_step(): i = %d, f = %f N, radius = %f m, rmass = %f kg, EI_partcile = %f kWh, ECS_particle = %f kWh/ton \n",i,sqrt(vectorMag3DSquared(f[i])),radius[i],rmass[i],EI_particle,ECS_particle);
							
								double temp_density = rmass[i] * 3.0 / (4.0 * pi * radius[i] * radius[i] * radius[i]);
								
									//			fprintf(screen ,"\n Extracted density_particle = %f from fix_template_multisphere.h \n", density_particle);
									//			fprintf(logfile ,"\n Extracted density_particle = %f from fix_template_multisphere.h \n", density_particle);										
								if(ECS_break <= ECS_particle && ECS_particle < ECS_max && radius[i] >= (min_parent_size_to_break / 2.0) && (0.9 * density_particle) <= temp_density && (temp_density <= 1.1 * density_particle))	//particle_density is coming from fix_template_sphere.h
									{
										flag[i] = 1;	
										flag_count++;
										rmass_break += rmass[i];
									////		if(screen) fprintf(screen ,"Flag [%d] turned to 1 with ECS_break (%f) <= ECS (%f) < ECS_max (%f) \n", i,ECS_break,ECS_particle,ECS_max);
									////		if(logfile) fprintf(logfile ,"Flag [%d] turned to 1 with ECS_break (%f) <= ECS (%f) < ECS_max (%f) \n", i,ECS_break,ECS_particle,ECS_max);
									}
									//current force acting on particle is greater than thresold force//
									//flag i=1 that will make that particle is eligible for breaking//
							}	
						else if(ECS_flag == 0)
							{
								f_break_sqr = f_break * f_break;   
								f_sqr = vectorMag3DSquared(f[i]); 
								//current force acting on particle//
								if(f_sqr >= f_break_sqr)	
								//current force acting on particle is greater than thresold force//
								flag[i] = 1;						
								//flag i=1 that will make that particle is eligible for breaking//
							}else 
							{
								error->all(FLERR,"Illegal ECS_flag value");
							}					
					}
			}
		
		if(ECS_flag == 1)
			{
				previous_particles_count = 0;
				rmass_total = 0.0;	
									
				double temp_dia = 0.0;
				for(int i = 0; i < nlocal; i++)
					{
						double temp_density = rmass[i] * 3.0 / (4.0 * pi * radius[i] * radius[i] * radius[i]);
								
						if(((0.9 * density_particle) < temp_density) && (temp_density < (1.1 * density_particle)))	
							{
								rmass_total += rmass[i];		//counts total mass in simulation (of particular density)
								previous_particles_count++;			//counts total number of particles in simulation (of particular density)							
							}		
					}
			}					
		
				
		if((ECS_flag == 1) && ((update->ntimestep == massdis_previous_ntimestep + massdis_nevery) || massdis_previous_ntimestep_flag == false))
			//for printing the mass distribution with mass distribution frequency (number of steps)
			{
				if(massdis_previous_ntimestep_flag == false)		//if this is being done for the first time
					{
						massdis_previous_ntimestep = update->ntimestep;
						massdis_previous_ntimestep_flag = true;
					}	
				
				mass_distribution = NULL;
				if(mass_distribution) memory->destroy(mass_distribution);
				mass_distribution = new double[sieves_series_length - index_lower - index_upper];	
				for(int i = 0; i < (sieves_series_length - index_lower - index_upper); i++) 	mass_distribution[i] = 0.0; 	//creating the zeros array of the size of the size of selected size
													
				double temp_dia = 0.0;
				for(int i = 0; i < nlocal; i++)		//adds mass to appropriate size for each particle
					{
						double temp_density = rmass[i] * 3.0 / (4.0 * pi * radius[i] * radius[i] * radius[i]);
								
						if(((0.9 * density_particle) < temp_density) && (temp_density < (1.1 * density_particle)))	
							{											
								temp_dia = 2.0 * radius[i];
								for(int j = 0; j < (sieves_series_length - index_lower - index_upper); j++)
									{
										if((0.95 * sieves[j] < temp_dia) && (temp_dia < 1.05 * sieves[j])) 
											{
												mass_distribution[j] += rmass[i];
												break;
											}
										if((j == (sieves_series_length - index_lower - index_upper) -1) && (1.05 * sieves[j] < temp_dia)) 	mass_distribution[j] += rmass[i];
									}																
							}		
					}
				
				for(int i = 1; i < (sieves_series_length - index_lower - index_upper); i++)	mass_distribution[i] += mass_distribution_truncated[i];			//adds the truncated mass	

				/*--------------------------Prints differential mass distribution----------------------------*/										
				fstream diff_massdistribution;
				diff_massdistribution.open("Diff_mass_distribution_file",std::ios::out | std::ios::app);			
			 
				if(massdis_flag_1 == false)
					{
						diff_massdistribution<<"\n \n \n\*--------------------------------------------------------------------------------------------------------------------------*\ \n";
						diff_massdistribution<<"\nntimestep \t Particles \t Differential mass distribution \n";
						massdis_flag_1 = true;
					}	
								
				writefile(diff_massdistribution, previous_particles_count, update->ntimestep, mass_distribution, sieves_series_length - index_lower - index_upper, 1); 
				diff_massdistribution.close();
				/*--------------------------------------------------------------------------------------------*/
								
				for(int i = 1; i < (sieves_series_length - index_lower - index_upper); i++)	mass_distribution[i] += mass_distribution[i-1];			//converts to cumulative mass distribution				
						
				/*-------------------------Prints cumulative mass distribution-------------------------------*/										
				fstream cum_massdistribution;
				cum_massdistribution.open("Cum_mass_distribution_file",std::ios::out | std::ios::app);			
			 
				if(massdis_flag_2 == false)
					{
						cum_massdistribution<<"\n \n \n\*--------------------------------------------------------------------------------------------------------------------------*\ \n";
						cum_massdistribution<<"\nntimestep \t Particles \t Cumulative mass distribution \n";
						massdis_flag_2 = true;
					}	
								
				writefile(cum_massdistribution, previous_particles_count, update->ntimestep, mass_distribution, sieves_series_length - index_lower - index_upper, 1);
				cum_massdistribution.close();	
				/*--------------------------------------------------------------------------------------------*/

				/*--------------Printing truncated mass to a file for use while restarting simulation----------------*/
				fstream MassDisTruncated;
					
				MassDisTruncated.open("mass_distribution_truncated", std::ios::out);	
											
				for(int j = 0; j < (sieves_series_length - index_lower - index_upper); j++)		
					{
						MassDisTruncated<<mass_distribution_truncated[j];
						MassDisTruncated<<"\t";
					}
									
				MassDisTruncated.close();	
				/*----------------------------------------------------------------------------------------------------*/

				massdis_previous_ntimestep = update->ntimestep;
			}			
					
		if(flag_count > 0)
			{
				if(screen) fprintf(screen, "\n \n \n \n \n\*-------------------------------------------------------------------------------------------------------------------------*\ \n \n ");
				if(logfile) fprintf(logfile, "\n \n \n \n \n\*-------------------------------------------------------------------------------------------------------------------------*\ \n \n");
				if(update->dt != previous_timestep)
					{
						if(screen) fprintf(screen, "timestep size = %.6g \n", update->dt);
						if(logfile) fprintf(logfile, "timestep size = %.6g \n", update->dt);
						previous_timestep = update->dt;	
					}
				if(screen) fprintf(screen, "\n%d flags (mass = %f kg) out of %d  turned to 1 at ntimestep = %d \n", flag_count,rmass_break,nlocal,update->ntimestep);
				if(logfile) fprintf(logfile, "\n%d flags (mass = %f kg) out of %d turned to 1 at ntimestep = %d \n", flag_count,rmass_break,nlocal,update->ntimestep);
			}
		if(screen) fprintf(screen ,"\n Exiting function end_of_step(). \n \n");
		if(logfile) fprintf(logfile ,"\n Exiting function end_of_step(). \n \n");	
	}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::pre_insert()		
	//things to be done before particle insertion//
	{
		if(screen) fprintf(screen ,"\n Entering function pre_insert(). \n");
		if(logfile) fprintf(logfile ,"\n Entering function pre_insert(). \n");
		
		int i,ibreak;
		int nlocal = atom->nlocal;	
			//nlocal is number of particles in current processor//
		int *mask = atom->mask;		
		double **x = atom->x;	
			//particles coordinates assignment//
		double **v = atom->v;		
			//particles forces assignments//
		double *radius = atom->radius;		 
			//particles radius assignments//
		double *rmass = atom->rmass;		
			//particles mass assignment//	
		double *flag = fix_break->vector_atom;		
			//assigning vector_atom from class FixPropertyAtom//
		AtomVec *avec = atom->avec;		
			//avec is a class pointer to AtomVec//also, declared in atom.h//
		double **f = atom->f;
			//force values//
		int *particle_type = atom->type;

		n_break_this_local = 0;
		mass_break_this_local = 0.;
		
	//	if(screen) fprintf(screen ,"\n");	
	//	if(logfile) fprintf(logfile ,"\n");
	//	if(screen) fprintf(screen ,"Real-time counter for n_break_this_local & mass break_this_local \t \n");	
	//	if(logfile) fprintf(logfile ,"Real-time counter for n_break_this_local & mass_break_this_local \t \n");
		
		for( i = 0; i < nlocal; i++)
			{
				if (mask[i] & groupbit && flag[i] == 1.) 
					//we have made flag 1 earlier for particles that need to be broken////mask is group number//
					{											
						//groupbit is the check whether particle belongs to current group or not//
						n_break_this_local++;
							//total mass of particle to be broken//rmass[i] is mass of particle based on r//
						mass_break_this_local += rmass[i];
							
			////		if(screen) fprintf(screen ,"i = %d (cum_n_break_this_local = %d, cum_mass_break_this_local = %f) \t \n",i,n_break_this_local,mass_break_this_local);	
			////		if(logfile) fprintf(logfile ,"i = %d (cum_n_break_this_local = %d, cum_mass_break_this_local = %f) \t \n",i,n_break_this_local,mass_break_this_local);							
					}
			}
	/*		
		if(n_break_this_local > 0)
			{
				if(screen) fprintf(screen ,"n_break_this_local = %d, mass_break_this_local = %f \t \n",n_break_this_local,mass_break_this_local);	
				if(logfile) fprintf(logfile ,"n_break_this_local = %d, mass_break_this_local = %f \t \n",n_break_this_local,mass_break_this_local);
			}	*/
	//	if(screen) fprintf(screen ,"\n");	
	//	if(logfile) fprintf(logfile ,"\n");
		
	//	if(screen) fprintf(screen ,"n_break_this_local = %d, mass_break_this_local = %f. \n \n",n_break_this_local,mass_break_this_local);	
	//	if(logfile) fprintf(logfile ,"n_break_this_local = %d, mass_break_this_local = %f. \n \n",n_break_this_local,mass_break_this_local);	
	

	   // tally stats

	   LAMMPS_NS::MPI_Sum_Scalar(n_break_this_local,n_break_this,world); 
			//defined in mpi.h//this function actually performs copying and adding the memeory into the destination///copies and add memory equivalent to int in n_break_this_local to n_break_this//
	   n_break += n_break_this;   
			//Getting total number of particles to be broken//
	  			
	   LAMMPS_NS::MPI_Sum_Scalar(mass_break_this_local,mass_break_this,world);    
			//copies memory equivalent to doubles in mass_break_this_local to mass_break_this//
	   mass_break += mass_break_this;
	   
	   if(n_break_this_local > 0)
			{
				if(screen) fprintf(screen ,"n_break_this_local = %d, mass_break_this_local = %f, n_break_this = %d, mass_break_this = %f, n_break = %d, mass_break = %f. \n \n",n_break_this_local,mass_break_this_local,n_break_this,mass_break_this,n_break,mass_break);	
				if(logfile) fprintf(logfile ,"n_break_this_local = %d, mass_break_this_local = %f, n_break_this = %d, mass_break_this = %f, n_break = %d, mass_break = %f. \n \n",n_break_this_local,mass_break_this_local,n_break_this,mass_break_this,n_break,mass_break);	
			}	

	   // allocate breakage data

	   //if(breakdata) memory->destroy_2d_double_array(breakdata);										
	   if(breakdata) memory->destroy(breakdata);
	   if(ECS) memory->destroy(ECS);
	   if(size) memory->destroy(size);
	   if(r_sphere) memory->destroy(r_sphere);
	   if(x_sphere) memory->destroy(x_sphere);
//	   if(volume_sphere) memory->destroy(volume_sphere);
//	   if(mass_sphere) memory->destroy(mass_sphere);
	   	
		//Please note that the appropriate function for this isn't available in memory.h
		//breakdata = memory->create_2d_double_array(n_break_this_local,7,"FixBreakparticleForce::breakdata");
	   //breakdata = memory->create(FixBreakparticleForce::breakdata,n_break_this_local,7,"FixBreakparticleForce::breakdata");
	   
	   breakdata = memory->create(FixBreakparticleForce::breakdata,n_break_this_local,9,"FixBreakparticleForce::breakdata");
	   
	   ECS = new double[n_break_this_local];
	   size = new double[n_break_this_local];
	 
	   	     //defined in memory.h//Its all about safe memory allocation//Its a 2d array (where each element can multiple values) with no of rows=n_break_this_local, no of columns=7, for 7 properties that are defined later in this file//
			//Please note that the appropriate function for this wasn't available in memory.h in the form of create(int,int,double **)//FixBreakparticleForce::breakdata was added which results 2d array as it can be seen at fix_breakparticle_force.h

	   // fill breakage data and remove particles
	   int nlocal_index = 0;
	   nlocal_index = nlocal;

	   i = ibreak = 0;
	   double EI_particle = 0.0; double ECS_particle = 0.0;
	   
	   while (i < nlocal) 	
		//nlocal is numbewr of particles in current processor//
		   {
				  if (mask[i] & groupbit && flag[i] == 1.)  
					//if turns true then data copy, ibreak++, nlocal-- otherwise i++//copying data of original particles//
						  {
										//copy data needed for insertion
								  vectorCopy3D(x[i],&breakdata[ibreak][0]);	
										//function defined in vector_liggghts.h//copying original parent coordinates from left to right//in a 2d array ibreak is row no and 0 is column no//copying coordinates in 1st column//
								  vectorCopy3D(v[i],&breakdata[ibreak][3]);	
										//copying from left to right//copying original parent vecoties in 4th column//
								  breakdata[ibreak][6] = radius[i];	
										//putting radius value in 7th column//
								  breakdata[ibreak][7] = rmass[i] * 3.0 / ( 4.0 * pi * radius[i] * radius[i] * radius[i]);	//density
							
									////	  fprintf(screen, "density = %f \n",rmass[i] * 3.0 / ( 4.0 * pi * radius[i] * radius[i] * radius[i]));
									////	  fprintf(screen, " breakdata[ibreak = %d][7] = %f \n",  ibreak, breakdata[ibreak][7]);
									////	  fprintf(logfile, "density = %f \n",rmass[i] * 3.0 / ( 4.0 * pi * radius[i] * radius[i] * radius[i]));
									////	  fprintf(logfile, " breakdata[ibreak = %d][7] = %f \n",  ibreak, breakdata[ibreak][7]);
										  
									////	  breakdata[ibreak][8] = particle_type[i];  //type
										  
									////	  breakdata[ibreak][9] = rmass[i];
															  
								  EI_particle = conversion_factor * sqrt(vectorMag3DSquared(f[i]));
									
								  ECS_particle = EI_particle / rmass[i] * 1000; 
										//Calulating ECS acting on the particle//In kWh/ton//	
													
								  ECS[ibreak] = ECS_particle;
									////	fprintf(screen, "ECS_particle = %f \n", ECS_particle);
									////	fprintf(logfile, "ECS_particle = %f \n", ECS_particle);
								  breakdata[ibreak][8] = ECS_particle;
								  
									////	  if(screen) fprintf(screen ,"deleting particles: i = %d, ibreak = %d, Size = %.8g m, EI = %.8g kWh, ECS = %.8g kWh/ton, f = %f N, rmass = %.8g kg,atom_type = %d, density = %f. \n",i,ibreak,size[ibreak],EI_particle,ECS[ibreak],sqrt(vectorMag3DSquared(f[ibreak])),rmass[ibreak],parent_particle_atom_type[ibreak], parent_particle_density[ibreak]  );
									////	  if(logfile) fprintf(logfile ,"deleting particles: i = %d, ibreak = %d, Size = %.8g m, EI = %.8g kWh, ECS = %.8g kWh/ton, f = %f N, rmass = %.8g kg,atom_type = %d, density = %f. \n",i,ibreak,size[ibreak],EI_particle,ECS[ibreak],sqrt(vectorMag3DSquared(f[ibreak])),rmass[ibreak],parent_particle_atom_type[ibreak], parent_particle_density[ibreak] );
																														
									////	if(screen) fprintf(screen ,"del par: i = %d, ibreak = %d, Size = %f m, density = %f, ECS = %f kWh/ton, f = %f N, rmass = %f kg \n",i,ibreak,2.0*breakdata[ibreak][6],breakdata[ibreak][7], breakdata[ibreak][8],sqrt(vectorMag3DSquared(f[i])),rmass[i]);
									////	if(logfile) fprintf(logfile ,"del par: i = %d, ibreak = %d, Size = %f m, density = %f, ECS = %f kWh/ton, f = %f N, rmass = %f kg \n",i,ibreak,2.0*breakdata[ibreak][6],breakdata[ibreak][7], breakdata[ibreak][8],sqrt(vectorMag3DSquared(f[i])),rmass[i]);																				
										 
								    // delete existing particle that is going to be broken//
								  avec->copy(nlocal-1,i,0);  
								  
								  /*
								   Please note that i is increasing from 0 to something and if it is found for deletion then it will be deleted.
								   So particle which deleted is not the last one in the list nlocal rather it is the current particle. On consequence the nlocal is reduced by 1. 
								   So don't connect i with the earlier sequence and therefore other values will be changed.
								  */
								  
										//command through atomvec.h//avec is a class pointer to AtomVec//after copying coordinates, velocities and radius of parent particle is going to be deleted.
								  nlocal--;					
										//after deleteing one particle new no of nlocal is obviously nlocal--//
										//i will remain 0 until if condition turns true//if all cases if is true then nlocal will automatically tends to 0//
										//if in any case if turns false then no data copy rather i++. so, i is the count of number of particles noy eligible for break//
										//So, i increases if particle found with ineligibility of breaking and nlocal reduces if particle with eligibility to be broken. Ultimately, both are tending to become equal at some point of time//
										//That means at that point i=nlocal and therefore false for while loop//
										
								  ibreak++;	
						  }
				  else i++;
		   }
		   
	   particles_count = previous_particles_count - ibreak;

	   if (nlocal != nlocal_index)
			{
		//		if(screen) fprintf(screen ,"\nprevious_particles_count = %d, ibreak = %d, particles_Count = %d \n",previous_particles_count,ibreak,particles_count);	
		//		if(logfile) fprintf(logfile ,"\nprevious_particles_count = %d, ibreak = %d, particles_Count = %d \n",previous_particles_count,ibreak,particles_count);	
				if(screen) fprintf(screen ,"\nnlocal reduced from %d (%d particles, mass = %f) to %d (%d particles, mass = %f). \n",nlocal_index,previous_particles_count,rmass_total,nlocal,particles_count,rmass_break);	
				if(logfile) fprintf(logfile ,"\nnlocal reduced from %d (%d particles, mass = %f) to %d (%d particles, mass = %f). \n",nlocal_index,previous_particles_count,rmass_total,nlocal,particles_count,rmass_break);	
			}
	   // update local and global # particles
	   
	   atom->nlocal = nlocal;	
			//updating no of particles in current processor//After several excutions of nlocal--//Particles that are going to be break have been deleted//
	   double rlocal = static_cast<double>(atom->nlocal);
	   MPI_Allreduce(&rlocal,&atom->natoms,1,MPI_DOUBLE,MPI_SUM,world);   
			//reducing memory in atom->natoms according to new nlocal//

	   // print stats
	   print_stats_breakage_during();     //print stats in terminal//
	   
	   if(screen) fprintf(screen ,"\n Exiting function pre_insert(). \n");
	   if(logfile) fprintf(logfile ,"\n Exiting function pre_insert(). \n");
	}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::print_stats_breakage_during()
	{
	  if(screen) fprintf(screen ,"\n Entering function print_stats_breakage_during(). \n");
	  if(logfile) fprintf(logfile ,"\n Entering function print_stats_breakage_during(). \n");
	  
	  int step = update->ntimestep;

	  if (me == 0 && n_break_this > 0)
		  {
				if (screen)
				  fprintf(screen ,"Particle breakage: broke %d particles (mass %f) at step %d\n - a total of %d particles (mass %f) broken so far \n \n",
					  n_break_this,mass_break_this,step,n_break,mass_break);	//n_break_this is no of particles broken currently whereas n_break is the total no of particles so far//

				if (logfile)
				  fprintf(logfile,"Particle breakage: broke %d particles (mass %f) at step %d\n - a total of %d particles (mass %f) broken so far \n \n",
					  n_break_this,mass_break_this,step,n_break,mass_break);
		  }
		  
	  if(screen) fprintf(screen ,"\n Exiting function print_stats_breakage_during(). \n");
	  if(logfile) fprintf(logfile ,"\n Exiting function print_stats_breakage_during(). \n");  
		  
	}

/* ---------------------------------------------------------------------- */

int FixBreakparticleForce::calc_ninsert_this()		
	{
		if(screen) fprintf(screen ,"\n Entering function calc_ninsert_this(). \n");
		if(logfile) fprintf(logfile ,"\n Entering function calc_ninsert_this(). \n");

/****************************************************************************************************************************************/		
		if(ECS_flag == 1)
			{			
				ninsert_daughter = 0;		//Number of total daughter particles to be inserted
				
				int iparticle = 0;		//index for the particle to be broken
				int cum_ninsert_daughter = 0;		//incrementing the daughters to-be inserted based on iparticle increment
				
				r_sphere = NULL;		//array of radius of the particles to-be inserted (one value for one daughter. only array will contain al daughter so repetition of sizes is quite obvious)
				if(r_sphere) memory->sfree(r_sphere);

				size_r_sphere = 0;		//size/length of r_sphere
				
				number_per_break = NULL;		//number of daughter particles to-be generated for each particle breakage
				if(number_per_break) memory->destroy(number_per_break);
				number_per_break = new int[n_break_this_local];			//assigning size to the number of particles to-be broken
				
				if(n_break_this_local > 0)
					{
						if(screen) fprintf(screen ,"\niparticle \t size \t \t mass \t \t ECS \t \t alpha \t \t t10%% \t max_daughter_size \t daughters_count \t daughters_mass \n");
						if(logfile) fprintf(logfile ,"\niparticle \t size \t \t mass \t \t ECS \t \t alpha \t \t t10%% \t max_daughter_size \t daughters_count \t daughters_mass \n");
					}
				
				while(iparticle < n_break_this_local)		//incrementing for each particle to break
						{
									int cum_ninsert_daughter_iparticle = 0;			//number of total fragements for ith particle
								
									size_iparticle = 2.0 * breakdata[iparticle][6];			//diameter of the particle to break
									ECS_iparticle = breakdata[iparticle][8];		//comminution energy currently working on it
									double mass_iparticle = 4.0 * pi * pow(size_iparticle,3) / 3.0 / 8.0 * density_particle;		//mass of the ith particle
																							
									alpha = C_GM * size_iparticle + alpha_intercept;		//alpha in the King model (alpha=C_ECS*ECS+C_GM*GM+C_WI*WI+C)
															
									t10 = (JK_parameter_A / 100.0) * (1.0-pow(e,((-JK_parameter_b) * ECS_iparticle)));		//t10 in King model (t10=A(1-exp(-b*ECS)))
															
									double t10_in_percent = t10 * 100.0;	
															
									int ii = 0; double Diff_tn = 0.0, cum_Diff_tn = 0.0; double cum_mass = 0.0;//ii is for incrementing sizes of particles//Diff_tn is the differential weightage of material in size class ii//
									int number_particles = 0; //n is the number of particle in i size class//r_n is radius of particle at that size class//
										
									//	if(screen) fprintf(screen ,"\n Sieves selected \t");
									//	if(logfile) fprintf(logfile ,"\n Sieves selected \t");
															
									int number_of_sieves_smaller = 0;
									while(sieves[ii] <= size_iparticle && ii < (sieves_series_length - index_lower - index_upper))		
										{																
										//	if(screen) fprintf(screen ,"%f \t",sieves[ii]);	
										//	if(logfile) fprintf(logfile ,"%f \t",sieves[ii]);	
											number_of_sieves_smaller++;			//counting the number of sieves sizes smaller than the current ith particle					
											ii++;
										}

									if(screen) fprintf(screen ,"\n Number of sieves selected = %d \n",number_of_sieves_smaller);
									if(logfile) fprintf(logfile ,"\n Number of sieves selected = %d \n",number_of_sieves_smaller);
															
									int number_particles_per_size_class[number_of_sieves_smaller];		//to record the number of particles per size smaller than the ith particle size
										
									//initializing values
									for(int count = 0; count < number_of_sieves_smaller; count++)
										{
											number_particles_per_size_class[count] = 0;
										}
									
									if(screen) fprintf(screen ,"\n size_class \t num_par \t Diff_tn \t  \t cum_Diff_tn \t cum_ninsert_daughter_iparticle \t cum_mass \n");
									if(logfile) fprintf(logfile ,"\n size_class \t num_par \t Diff_tn \t  \t cum_Diff_tn \t cum_ninsert_daughter_iparticle \t cum_mass \n");
									
									bool exact_match_max_size = true;		//true if the top size of the number_particles_per_size_class array equals to the ith particle size
									if(size_iparticle > (1.01 * sieves[number_of_sieves_smaller - 1])) exact_match_max_size = false;
																						
									ii = 0;			
																									
									while(sieves[ii] <= size_iparticle && ii < (sieves_series_length - index_lower - index_upper))     //fragments counting for each size based on t10 model
										{
												double n_current = 0.0, n_previous = 0.0, n_next = 0.0, power_current = 0.0, power_previous = 0.0, power_next = 0.0;
																															  
												if(ii == 0)
													{
														n_current = size_iparticle / sieves[ii];
														assert(n_current != 1.0);
																						
														power_current = pow((9.0/(n_current-1.0)),alpha);	
																													
														Diff_tn = 1.0-pow((1.0-t10),power_current);
														
														assert(Diff_tn >= 0.0);
													}else if(ii > 0 && ii < number_of_sieves_smaller)
													{
														n_previous = size_iparticle / sieves[ii-1];
														n_current = size_iparticle / sieves[ii];
																						
														power_previous = pow((9.0/(n_previous-1.0)),alpha);
														power_current = pow((9.0/(n_current-1.0)),alpha);
											
														Diff_tn =  (1.0-pow((1.0-t10),power_current)) - (1.0-pow((1.0-t10),power_previous));
														
														if(ii == (number_of_sieves_smaller - 1) && exact_match_max_size == false)
															{
																n_previous = size_iparticle / sieves[ii];
																power_previous = pow((9.0/(n_previous-1.0)),alpha);
																Diff_tn += 1 - (1.0-pow((1.0-t10),power_previous));
															}														
														assert(Diff_tn >= 0.0);																																					
													}else 
													{
														error->all(FLERR,"Illegal ECS_flag value");
													}
																				
												cum_Diff_tn += Diff_tn;

												number_particles = (int)round((Diff_tn * pow((size_iparticle / sieves[ii]),3)));	//Differential tn value is multiplies by (R/r)^3 to get the number of particles in that size. Combination of (int)round() will give nearest integer
																				
												number_particles_per_size_class[ii] = number_particles;			//assigning number particles to-be generated to the appropriate size in a pre-defined array
												
												cum_mass += number_particles * 4.0 * pi * pow(sieves[ii],3) / 3.0 / 8.0 * density_particle;
																				
												cum_ninsert_daughter_iparticle += number_particles;		//incrementing the number of particles to-be generated for ith particle
																																
												if(screen) fprintf(screen ," %f \t \t %d \t \t %f \t \t %f \t \t \t %d \t \t \t \t %f \n",sieves[ii], number_particles_per_size_class[ii], Diff_tn,cum_Diff_tn,cum_ninsert_daughter_iparticle, cum_mass);
												if(logfile) fprintf(logfile ," %f \t \t %d \t \t %f \t \t %f \t \t \t %d \t \t \t \t %f \n",sieves[ii], number_particles_per_size_class[ii], Diff_tn,cum_Diff_tn,cum_ninsert_daughter_iparticle, cum_mass);
																																
												ii++;
										}

									/*-----------Adjustments for the mass loss due to round-off calculation errors----------*/
									if(cum_Diff_tn < 1.0)	
									{
										double adjustment_particles = (int)round(((1-cum_Diff_tn) * pow((size_iparticle / sieves[0]),3)));		//number particles last being assumed as belonging to smallest size
										number_particles_per_size_class[0] += adjustment_particles;			//lost mast is being adjusted to smallest size or fines
									}
									/*--------------------------------------------------------------------------------------*/
																														
									/**************Resizing the r_sphere array by cum_ninsert_daughter_iparticle*************/
									if(!r_sphere) 
										{
											r_sphere = new double[cum_ninsert_daughter_iparticle];
										}else
										{											
											double *resize_array = new double[size_r_sphere + cum_ninsert_daughter_iparticle];
		
											for(int i = 0; i < size_r_sphere; i++)
												{
													resize_array[i] = r_sphere[i];
												}	
											
											r_sphere = resize_array;
											resize_array = NULL;
											if(resize_array) memory->destroy(resize_array);
											if(resize_array) memory->sfree(resize_array);
											delete[] resize_array;	
										}
									/****************************************************************************************/	
									
									int jj=0; int kk=0; int mm=0; int count=0; 
									
									double max_daughter_size = sieves[number_of_sieves_smaller - 1];									
																								
									for(jj=0; jj < number_of_sieves_smaller; jj++)			//putting radius values to the r_sphere array. One value for one particle. So, obvious repetition of the size values.
										{							
											for(kk=0; kk < number_particles_per_size_class[jj]; kk++)
												{																				
													r_sphere[size_r_sphere + count + kk] = sieves[jj] / 2.0;
												}																			
											count = count + number_particles_per_size_class[jj];
										}
									
									size_r_sphere += cum_ninsert_daughter_iparticle;		//incrementing the index to show size of the r_sphere array to match the current size of r_sphere
									number_per_break[iparticle] = cum_ninsert_daughter_iparticle;		//assigning fragments number to appropriate size of number_per_break array
									cum_ninsert_daughter += cum_ninsert_daughter_iparticle;			//incrementing the total fragments to-be generated

									if(screen) fprintf(screen ,"%d \t \t %f \t %f \t %f \t %f \t %f \t %f \t \t %d \t \t \t %f \n \n", iparticle, size_iparticle,mass_iparticle,ECS_iparticle,alpha,t10_in_percent,max_daughter_size, cum_ninsert_daughter_iparticle, cum_mass);
									if(logfile) fprintf(logfile ,"%d \t \t %f \t %f \t %f \t %f \t %f \t %f \t \t %d \t \t \t %f \n \n", iparticle, size_iparticle,mass_iparticle,ECS_iparticle,alpha,t10_in_percent,max_daughter_size,cum_ninsert_daughter_iparticle, cum_mass);

									for(int mm=0; mm<cum_ninsert_daughter; mm++)
										{
											if(screen) fprintf(screen ,"r_sphere[%d] = %f \t", mm, r_sphere[mm]);
											if(logfile) fprintf(logfile ,"r_sphere[%d] = %f \t", mm, r_sphere[mm]);
										}
									
									if(screen) fprintf(screen ,"Total number of fragments to insert in this loop = %d", cum_ninsert_daughter);
									if(logfile) fprintf(logfile ,"Total number of fragments to insert in this loop = %d", cum_ninsert_daughter);
										
									iparticle++;	//while loop increment//
						}
				
/********************************************************************************************************************************************************/
				ninsert_daughter = cum_ninsert_daughter;		//ninsert_daughter is the total number of fragments to-be generated for breakage of all eligible parent particles
			}

		if(screen) fprintf(screen ,"\nn_break_this = %d \n",n_break_this);
		if(logfile) fprintf(logfile ,"\nn_break_this = %d \n",n_break_this);
	
		if(screen) fprintf(screen ,"\n Exiting function calc_ninsert_this(). \n");
		if(logfile) fprintf(logfile ,"\n Exiting function calc_ninsert_this(). \n");

		// number of ptis to insert this timestep
		// will effectively insert n_break_this * n_fragments spheres
	
		return n_break_this;

	}

/* ----------------------------------------------------------------------
   generate new particles at positions where old particles were deleted
   function is executed locally on each process as opposed to
   FixInsertPack::x_v_omega()

   overlap check is not needed since space around broken particles is empty

   returns # bodies and # spheres that could actually be inserted
------------------------------------------------------------------------- */

void FixBreakparticleForce::x_v_omega(int ninsert_this,int &ninserted_this, int &ninserted_spheres_this, double &mass_inserted_this)
	{
		if(screen) fprintf(screen ,"\n \n \n \n \n \n Entering function x_v_omega() \n");
		if(logfile) fprintf(logfile ,"\n \n \n \n \n \n Entering function x_v_omega() \n");
		
		double pos_ins[3],v_ins[3],omega_ins[3],quat_ins[4],rad_broken;
		int iparticle, nins;
		ParticleToInsert *pti;		
			//ParticleToInsert is a class declared in particleToInsert.h//

		vectorZeroize3D(omega_ins);	
			//turning all values to zero//
		vectorZeroize4D(quat_ins);
			//turning all values to zero//

		double mass_inserted_this_local = 0.;		// local insertion for local processor
		mass_inserted_this = 0.;					// global insertion will be added by mass_inserted_this_local. variable is in function parameter
				
		int ninserted_spheres_this_local = 0;		// local insertion for local processor
		ninserted_spheres_this = 0;					// global insertion will be added by ninserted_spheres_this_local. variable is in function parameter
		
		int ninserted_this_local = 0;				// local insertion
		ninserted_this = 0;							// global insertion will be added by ninserted_this_local. variable is in function parameter
				
		iparticle = 0; 		//ith particle index

		if(ECS_flag == 0)
				{
						//	if(screen) fprintf(screen ,"\n In while loop of function x_v_omega() of ECS_flag = %d \n \n",ECS_flag);
						//	if(logfile) fprintf(logfile ,"\n In while loop of function x_v_omega() of ECS_flag = %d \n \n",ECS_flag);
							
						//	if(screen) fprintf(screen ,"\n A \n");
						//	if(logfile) fprintf(logfile ,"\n A \n");	
							
							while(iparticle < n_break_this_local)
									{										
													//		if(screen) fprintf(screen ,"\n B \n");
													//		if(logfile) fprintf(logfile ,"\n B \n");	
											vectorCopy3D(&breakdata[iparticle][0],pos_ins);	
													//copies positions from left to right//
													//		if(screen) fprintf(screen ,"\n B- \n");
													//		if(logfile) fprintf(logfile ,"\n B- \n");	
											vectorCopy3D(&breakdata[iparticle][3],v_ins);	
													//copies velocities from left to right//
													//		if(screen) fprintf(screen ,"\n B-- \n");
													//		if(logfile) fprintf(logfile ,"\n B-- \n");	
											rad_broken = breakdata[iparticle][6];			
													//copies original particle radius to rad_broken//
													// get pti and scale it down with radius of broken particle
													//		if(screen) fprintf(screen ,"\n B--- \n");
													//		if(logfile) fprintf(logfile ,"\n B--- \n");	
											pti = fix_distribution->pti_list[iparticle];  
													//	if(screen) fprintf(screen ,"\n B---- \n");
													//	if(logfile) fprintf(logfile ,"\n B---- \n");	

													//	if(screen) fprintf(screen ,"\n C \n");
													//	if(logfile) fprintf(logfile ,"\n C \n");	 
													//fixdistribution is a class pointer to FixParticledistributionDiscrete (declared in fix_insert.h, class defined in fix_particledistribution_discrete.h)//
													//pti_list is double array of class ParticleToInsert (defined in fix_particledistribution_discrete.h, class defined in particletoinsert.h)//
													//See in particleToInsert.h where class ParticleToInsert is defined//
													//pti_list contains details of particle properties of one particle, pti is class pointer of particletoinsert//
													//so, particle properties from fix_distribution or indirectly from fix_particledistribution_discrete.h is being transferred to ParticleToInsert class//																	
											pti->scale_pti(rad_broken);  
													//scaling the fragments obtained from fragment file in the simulation folder to the actual parent size (fragments given in file are for parent size of 1 unit)
													//	if(screen) fprintf(screen ,"\n D \n");
													//	if(logfile) fprintf(logfile ,"\n D \n");	
													//rad_broken should be a relative values such 0.3 or 0.5//
													//scale_pti is declared in particletoinsert.h//function defined in particletoinsert.cpp//
													//scale_pti(r_scale) multiplies poistions,velocities and r_bound_ins with r_scale, and volume and mass with r_scale^3//
													//These are relative operations and actual values will be set up in next line//
																		
											nins = pti->set_x_v_omega(pos_ins,v_insert,omega_ins,quat_insert); 	//function defined in particletoinsert.cpp//
											
													//	if(screen) fprintf(screen ,"\n E \n");
													//	if(logfile) fprintf(logfile ,"\n E \n");	
													//Copies velocities and omega. Returns values nspheres to nins. Add relative values to position//
						
											// tally stats
											ninserted_spheres_this_local += nins;  
													//Please note it was 0 earlier so now  ninserted_spheres_this_local=nins//
											mass_inserted_this_local += pti->mass_ins;
													//Again mass_inserted_this_local=mass_ins
											ninserted_this_local++;	
													//increment for number of particles inserted////
											iparticle++;
											
										//	if(screen) fprintf(screen ,"\n F \n");
										//	if(logfile) fprintf(logfile ,"\n F \n");	
//											if(screen) fprintf(screen ,"\n while loop ending here \t \n \n ");		
//											if(logfile) fprintf(logfile ,"\n while loop ending here \t \n \n ");
									}	
					
				
				}else if(ECS_flag == 1)
				{					
						  if(screen) fprintf(screen ,"\nn_break_this_local = %d \n \n \n",n_break_this_local);
						  if(logfile) fprintf(logfile ,"\nn_break_this_local = %d \n \n \n",n_break_this_local);	
						  
						  if(screen) fprintf(screen ,"\nNumber of daughters per break \n");
						  if(logfile) fprintf(logfile ,"\nNumber of daughters per break \n");
						  					
						  int total_break = 0;
						  
						  while(iparticle < n_break_this_local)		//calculates the total number of fragments to-be generated
						    	{
									if(screen) fprintf(screen ,"%d \t",number_per_break[iparticle]);
									if(logfile) fprintf(logfile ,"%d \t",number_per_break[iparticle]);
									total_break += number_per_break[iparticle];
									iparticle++;
							    }

						  iparticle = 0;			  
						  if(screen) fprintf(screen ,"\nTotal = %d \n  \n",total_break);
						  if(logfile) fprintf(logfile ,"\nTotal = %d \n  \n",total_break);
						  										
						  int daughter = 0;
						  int previous_sum = 0;
						  					  						  
						  while(iparticle < n_break_this_local)				
								{
										previous_sum = daughter;
										
										for(daughter; daughter < (previous_sum + number_per_break[iparticle]); daughter++)
											{
												if(r_sphere[daughter] > (min_parent_size_to_break / 2.0))
													{
														for(int iii=0; iii<3; iii++)
															{
																pos_ins[iii] = random_coordinates_generator(iii, update->ntimestep * random * 554);
																random++;
															}
																			//		if(screen) fprintf(screen ,"\n daughter + 1 = %d, number_per_break[%d] = %d, previous_sum + 1 + number_per_break[%d]) = %d, coordinates(%f,%f,%f) \n",daughter+1,iparticle,number_per_break[iparticle],iparticle,previous_sum + number_per_break[iparticle],pos_ins[0],pos_ins[1],pos_ins[2]);
																			//		if(logfile) fprintf(logfile ,"\n daughter + 1 = %d, number_per_break[%d] = %d, previous_sum + 1 + number_per_break[%d]) = %d, coordinates(%f,%f,%f) \n",daughter+1,iparticle,number_per_break[iparticle],iparticle,previous_sum + number_per_break[iparticle],pos_ins[0],pos_ins[1],pos_ins[2]);
											  			vectorCopy3D(&breakdata[iparticle][3],v_ins);	
																			//copies velocities from left to right//
																			//rad_broken = breakdata[iparticle][6];			
																			//copies original particle radius to rad_broken//
																																				
														/********************Most important part that controls the creation of new particles***********************/
																			// get pti and scale it down with radius of broken particle
														pti = fix_distribution->pti_list[daughter];   
																			//fixdistribution is a class pointer to FixParticledistributionDiscrete (declared in fix_insert.h, class defined in fix_particledistribution_discrete.h)//
																			//pti_list is double array of class ParticleToInsert (defined in fix_particledistribution_discrete.h, class defined in particletoinsert.h)//
																			//See in particleToInsert.h where class ParticleToInsert is defined//
																			//pti_list contains details of particle properties of one particle, pti is class pointer of particletoinsert//
																			//so, particle properties from fix_distribution or indirectly from fix_particledistribution_discrete.h is being transferred to ParticleToInsert class//
																			
														if(screen) fprintf(screen ,"r_sphere[%d] = %f \n", daughter, r_sphere[daughter]);		
														if(logfile) fprintf(logfile ,"r_sphere[%d] = %f \n", daughter, r_sphere[daughter]);

														pti->set_r_mass_vol_rboundins(r_sphere[daughter]);
																								
														pti->scale_pti(1.0);  
																			//rad_broken should be a relative values such 0.3 or 0.5//
																			//scale_pti is declared in particletoinsert.h//function defined in particletoinsert.cpp//
																			//scale_pti(r_scale) multiplies poistions,velocities and r_bound_ins with r_scale, and volume and mass with r_scale^3//
																			//These are relative operations and actual values will be set up in next line//
																													
														nins = pti->set_x_v_omega(pos_ins,v_insert,omega_ins,quat_insert); 	//function defined in particletoinsert.cpp//
																			//Copies velocities and omega. Returns values nspheres to nins. Add relative values to position//
														/********************************************************************************************************/
																							
																			// tally stats
														ninserted_spheres_this_local += nins;  
																			//Please note it was 0 earlier so now  ninserted_spheres_this_local=nins//
														mass_inserted_this_local += pti->mass_ins;
																			//Again mass_inserted_this_local=mass_ins	
													}else
													{
														if(passing_truncation_switch == true)
															{
																double fragment_size = 2.0 * r_sphere[daughter];
																int jj = 0;
																while((sieves[jj] < fragment_size) && (jj < (sieves_series_length - index_lower - index_upper)))
																	{
																		jj++;
																	}
																mass_distribution_truncated[jj] += density_particle * 4.0 * pi * r_sphere[daughter] * r_sphere[daughter] * r_sphere[daughter] / 3.0;
															}else
															{
																if(regenerate_in_bin_switch == false)
																	{
																		for(int iii=0; iii<3; iii++)
																			{
																				pos_ins[iii] = random_coordinates_generator(iii, update->ntimestep * random * 2554);
																				random++;
																			}
																							//		if(screen) fprintf(screen ,"\n daughter + 1 = %d, number_per_break[%d] = %d, previous_sum + 1 + number_per_break[%d]) = %d, coordinates(%f,%f,%f) \n",daughter+1,iparticle,number_per_break[iparticle],iparticle,previous_sum + number_per_break[iparticle],pos_ins[0],pos_ins[1],pos_ins[2]);
																							//		if(logfile) fprintf(logfile ,"\n daughter + 1 = %d, number_per_break[%d] = %d, previous_sum + 1 + number_per_break[%d]) = %d, coordinates(%f,%f,%f) \n",daughter+1,iparticle,number_per_break[iparticle],iparticle,previous_sum + number_per_break[iparticle],pos_ins[0],pos_ins[1],pos_ins[2]);
															  			vectorCopy3D(&breakdata[iparticle][3],v_ins);	
																							//copies velocities from left to right//
																							//rad_broken = breakdata[iparticle][6];			
																							//copies original particle radius to rad_broken//
																																								
																		/********************Most impornat part that controls the creation of new particles***********************/
																							// get pti and scale it down with radius of broken particle
																		pti = fix_distribution->pti_list[daughter];   
																							//fixdistribution is a class pointer to FixParticledistributionDiscrete (declared in fix_insert.h, class defined in fix_particledistribution_discrete.h)//
																							//pti_list is double array of class ParticleToInsert (defined in fix_particledistribution_discrete.h, class defined in particletoinsert.h)//
																							//See in particleToInsert.h where class ParticleToInsert is defined//
																							//pt+i_list contains details of particle properties of one particle, pti is class pointer of particletoinsert//
																							//so, particle properties from fix_distribution or indirectly from fix_particledistribution_discrete.h is being transferred to ParticleToInsert class//
																							
																		pti->set_r_mass_vol_rboundins(r_sphere[daughter]);
																												
																		pti->scale_pti(1.0);  
																							//rad_broken should be a relative values such 0.3 or 0.5//
																							//scale_pti is declared in particletoinsert.h//function defined in particletoinsert.cpp//
																							//scale_pti(r_scale) multiplies poistions,velocities and r_bound_ins with r_scale, and volume and mass with r_scale^3//
																							//These are relative operations and actual values will be set up in next line//
																																	
																		nins = pti->set_x_v_omega(pos_ins,v_insert,omega_ins,quat_insert); 	//function defined in particletoinsert.cpp//
																							//Copies velocities and omega. Returns values nspheres to nins. Add relative values to position//
																		/********************************************************************************************************/
																											
																							// tally stats
																		ninserted_spheres_this_local += nins;  
																							//Please note it was 0 earlier so now  ninserted_spheres_this_local=nins//
																		mass_inserted_this_local += pti->mass_ins;
																							//Again mass_inserted_this_local=mass_ins	
																	}else
																	{
																		for(int iii=0; iii<3; iii++)
																			{
																				pos_ins[iii] = random_coordinates_generator_bin(iii, update->ntimestep * random * 5554);
																				random++;
																			}
																							//		if(screen) fprintf(screen ,"\n daughter + 1 = %d, number_per_break[%d] = %d, previous_sum + 1 + number_per_break[%d]) = %d, coordinates(%f,%f,%f) \n",daughter+1,iparticle,number_per_break[iparticle],iparticle,previous_sum + number_per_break[iparticle],pos_ins[0],pos_ins[1],pos_ins[2]);
																							//		if(logfile) fprintf(logfile ,"\n daughter + 1 = %d, number_per_break[%d] = %d, previous_sum + 1 + number_per_break[%d]) = %d, coordinates(%f,%f,%f) \n",daughter+1,iparticle,number_per_break[iparticle],iparticle,previous_sum + number_per_break[iparticle],pos_ins[0],pos_ins[1],pos_ins[2]);
															  			vectorCopy3D(&breakdata[iparticle][3],v_ins);	
																							//copies velocities from left to right//
																							//rad_broken = breakdata[iparticle][6];			
																							//copies original particle radius to rad_broken//
																																								
																		/********************Most impornat part that controls the creation of new particles***********************/
																							// get pti and scale it down with radius of broken particle
																		pti = fix_distribution->pti_list[daughter];   
																							//fixdistribution is a class pointer to FixParticledistributionDiscrete (declared in fix_insert.h, class defined in fix_particledistribution_discrete.h)//
																							//pti_list is double array of class ParticleToInsert (defined in fix_particledistribution_discrete.h, class defined in particletoinsert.h)//
																							//See in particleToInsert.h where class ParticleToInsert is defined//
																							//pti_list contains details of particle properties of one particle, pti is class pointer of particletoinsert//
																							//so, particle properties from fix_distribution or indirectly from fix_particledistribution_discrete.h is being transferred to ParticleToInsert class//
																							
																		pti->set_r_mass_vol_rboundins(r_sphere[daughter]);
																												
																		pti->scale_pti(1.0);  
																							//rad_broken should be a relative values such 0.3 or 0.5//
																							//scale_pti is declared in particletoinsert.h//function defined in particletoinsert.cpp//
																							//scale_pti(r_scale) multiplies poistions,velocities and r_bound_ins with r_scale, and volume and mass with r_scale^3//
																							//These are relative operations and actual values will be set up in next line//
																																	
																		nins = pti->set_x_v_omega(pos_ins,v_insert,omega_ins,quat_insert); 	//function defined in particletoinsert.cpp//
																							//Copies velocities and omega. Returns values nspheres to nins. Add relative values to position//
																		/********************************************************************************************************/
																											
																							// tally stats
																		ninserted_spheres_this_local += nins;  
																							//Please note it was 0 earlier so now  ninserted_spheres_this_local=nins//
																		mass_inserted_this_local += pti->mass_ins;
																							//Again mass_inserted_this_local=mass_ins		
																	}
															}
													}		
											}
																											
									ninserted_this_local++;	
										//increment for number of particles inserted////
									iparticle++;
											
									//if(screen) fprintf(screen ,"\n while loop ending here \t \n \n ");		
									//if(logfile) fprintf(logfile ,"\n while loop ending here \t \n \n ");										
								}																									
																																		
				} else
				{  
					error->all(FLERR,"Illegal ECS_flag value");
				}												
		
		// tally stats, have to do this since operation is locally on each process
		// as opposed to e.g. FixInsertPack::x_v_omega()
		
	//	if(screen) fprintf(screen ,"\n G \n");
     //		  if(logfile) fprintf(logfile ,"\n G \n");	

		LAMMPS_NS::MPI_Sum_Scalar(ninserted_spheres_this_local,ninserted_spheres_this,world);	//copies memory//
	//	if(screen) fprintf(screen ,"\n H \n");
     //		  if(logfile) fprintf(logfile ,"\n H \n");
		LAMMPS_NS::MPI_Sum_Scalar(ninserted_this_local,ninserted_this,world);
	//	if(screen) fprintf(screen ,"\n I \n");
    // 		  if(logfile) fprintf(logfile ,"\n I \n");
		LAMMPS_NS::MPI_Sum_Scalar(mass_inserted_this_local,mass_inserted_this,world);																														
	//		if(screen) fprintf(screen ,"\n J \n");
     //		  if(logfile) fprintf(logfile ,"\n J \n");													
														
				//increment for number of particles already broken////
//		if(screen) fprintf(screen ,"\n ninserted_spheres_this_local = %d, mass_inserted_this_local= %f \t \n ", ninserted_spheres_this_local, mass_inserted_this_local);		
//		if(logfile) fprintf(logfile ,"\n ninserted_spheres_this_local = %d, mass_inserted_this_local= %f \t \n ", ninserted_spheres_this_local, mass_inserted_this_local);	
							
	    if(screen) fprintf(screen ,"\n Exiting function x_v_omega(). \n ");
	    if(logfile) fprintf(logfile ,"\n Exiting function x_v_omega(). \n "); 
	  
//	    if(screen) fprintf(screen ,"\n *********************************************************************** \n \n \n \n \n \n \n \n \n \n");
//	    if(logfile) fprintf(logfile ,"\n ***********************************************************************. \n \n \n \n \n \n \n \n \n \n"); 
	}

//not sure if we realy need this function, however it needs to be a function with this name here for the parent class fix-instert
double FixBreakparticleForce::insertion_fraction()
	{
		if(screen) fprintf(screen ,"\n Entering function insertion_fraction(). \n");
		if(logfile) fprintf(logfile ,"\n Entering function insertion_fraction(). \n");	
		
		if(screen) fprintf(screen ,"\n Exiting function insertion_fraction(). \n");
		if(logfile) fprintf(logfile ,"\n Exiting function insertion_fraction(). \n");
		
		return volumefraction;
	}
