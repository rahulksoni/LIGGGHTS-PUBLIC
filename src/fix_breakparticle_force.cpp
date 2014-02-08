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
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"										/*************/
#include "atom_vec.h"									/************/
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "domain.h"
#include "random_park.h"								/***********/
#include "memory.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_particledistribution_discrete.h"			/**************/
#include "fix_template_multiplespheres.h"				/************/
#include "particleToInsert.h"							/************/


/***************************************************************/
#include "fix_breakparticle_force.h"
#include "breakage_coordinates.h"
#include "assert.h"
/***************************************************************/

/*************************For reading data from files*************************************/
#include <iostream>
#include "stdio.h"
#include <fstream>
#include <limits>
#include <cmath>
#include <cstdlib>		//for accessing rand number generator function
#include <ctime>		//for accessing current time to create a seed value for the random number generator
/*****************************************************************************************/

using namespace LAMMPS_NS;
using namespace FixConst; //needed for END_OF_STEP

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
/*---------------------------------------------------------------------------------------------*/

FixBreakparticleForce::FixBreakparticleForce(LAMMPS *lmp, int narg, char **arg) :
  FixInsert(lmp, narg, arg)   
				//FixBreakparticleForce construction, inclusion of FixInsert, narg is no of arguments in fixbreakparticle/force command, **arg is corresponding keywords//
	{
//		if(screen) fprintf(screen ,"\n \n \n \n Entering function FixBreakparticleForce(). \n");
//		if(logfile) fprintf(logfile ,"\n \n \n \n Entering function FixBreakparticleForce(). \n");
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
					int sieves_series_length = 40;
					
					fstream parameters;
					
					/**********************************************************************************************/
					parameters.open("parametersfile");

					/*--------------------------------------------------------------------------------------------*/
					GotoLine(parameters,5);

					parameters>>JK_parameter_A;
					parameters>>JK_parameter_b;
					parameters>>WI;
					parameters>>C_GM;
					parameters>>C_WI;
					parameters>>C_intercept;
					/*--------------------------------------------------------------------------------------------*/
					
					/*--------------------------------------------------------------------------------------------*/					
					GotoLine(parameters,9);
					
					parameters>>min_parent_size_to_break;
					parameters>>min_daughter_size;		
					/*--------------------------------------------------------------------------------------------*/
					
					/*--------------------------------------------------------------------------------------------*/								
					GotoLine(parameters,13);
					
					parameters>>ECS_max;
					parameters>>con_fac_joule;
					parameters>>con_fac_force;
					/*--------------------------------------------------------------------------------------------*/
					GotoLine(parameters,17);
					
					parameters>>x_min;
					parameters>>x_max;
					parameters>>y_min;
					parameters>>y_max;
					parameters>>z_min;
					parameters>>z_max;
					parameters>>mill_axis;
					/*-------------------------------------------------------------------------------------------*/
					
					parameters.close();
					/*********************************************************************************************/
					
					conversion_factor = (con_fac_joule / con_fac_force) * 2.77777778 / 10000000.0;	// = 7.3491e-10	
						//Estimating EI on particle based on force comparison with figure 6 of Nikhil's paper//In KWh//
						//1757.08 N is force back-calculated for EI = 4.65 Joule from Nikhil's paper//
						//"E.T. Tuzcu, N. Dhawan and R.K. Rajamani. Coarse particle fracture with the ultrafast load cell. Minerals & Metallurgical Processing, 2011, Vol. 28, No. 4, pp. 176-186"//
						
					alpha_intercept = C_WI * WI + C_intercept;	
									
					if(screen) fprintf(screen ,"JK_parameter_A = %f \t JK_parameter_b = %f \t WI = %f \t C_GM = %f \t C_WI = %f \t C_intercept = %f \t alpha_intercept = %f \n \n", JK_parameter_A,JK_parameter_b,WI,C_GM,C_WI,C_intercept,alpha_intercept);
					if(logfile) fprintf(logfile ,"JK_parameter_A = %f \t JK_parameter_b = %f \t WI = %f \t C_GM = %f \t C_WI = %f \t C_intercept = %f \t alpha_intercept = %f \n \n", JK_parameter_A,JK_parameter_b,WI,C_GM,C_WI,C_intercept,alpha_intercept);

					if(screen) fprintf(screen , "min_parent_size_to_break = %f \t min_daughter_size = %f \n \n ",min_parent_size_to_break,min_daughter_size);
					if(logfile) fprintf(logfile , "min_parent_size_to_break = %f \t min_daughter_size = %f \n \n ",min_parent_size_to_break,min_daughter_size);
					
					if(screen) fprintf(screen , "ECS_max = %f \t conversion_factor_energy = %f J \t conversion_factor_force = %f N \t conversion_factor = %.6g \n \n ",ECS_max,con_fac_joule,con_fac_force,conversion_factor);
					if(logfile) fprintf(logfile , "ECS_max = %f \t conversion_factor_energy = %f J \t conversion_factor_force = %f N \t conversion_factor = %.6g \n \n ",ECS_max,con_fac_joule,con_fac_force,conversion_factor);
					
					if(screen) fprintf(screen, "x_min = %f \t x_max = %f \t y_min = %f \t y_max = %f \t z_min = %f \t z_max = %f \t mill_axis = %d \n \n", x_min,x_max,y_min,y_max,z_min,z_max,mill_axis);
					if(logfile) fprintf(logfile, "x_min = %f \t x_max = %f \t y_min = %f \t y_max = %f \t z_min = %f \t z_max = %f \t mill_axis = %d \n \n", x_min,x_max,y_min,y_max,z_min,z_max,mill_axis);

					int index = 0;
					while(sieves_series[index] < min_daughter_size)
						{
							index++;
						}

					sieves = NULL;
					if(sieves) memory->destroy(sieves);
					sieves = new double[sieves_series_length - index];

					if(screen) fprintf(screen ,"\n sieves \n");
					if(logfile) fprintf(logfile ,"\n sieves \n");
					for(int j=0; j < (sieves_series_length - index); j++)
						{
							sieves[j] = sieves_series[index + j];
							if(screen) fprintf(screen ,"%f \t", sieves[j]);
							if(logfile) fprintf(logfile ,"%f \t", sieves[j]);
						}		
						
					if(screen) fprintf(screen ,"\n \n ");
					if(logfile) fprintf(logfile ,"\n \n ");
						
					fill_bounds(mill_axis, x_min, x_max, y_min, y_max, z_min, z_max);		
					
					/*--------------------------------------------------------------------------------------------------------------------*/					
					assert(JK_parameter_A > 1.0 && JK_parameter_A < 100.0);
					assert(JK_parameter_b > 0.0 && JK_parameter_b < 10.0);
					assert(WI > 0.0 && WI < 100.0);
					assert(ECS_max > ECS_break && ECS_max < 10.0);
					assert(min_daughter_size > sieves_series[0] && min_daughter_size < min_parent_size_to_break);
					assert(con_fac_joule > 0.0);
					assert(con_fac_force > 0.0);
					assert(x_max > x_min);
					assert(y_max > y_min);
					assert(z_max > z_min);
					assert(mill_axis == 1 || mill_axis == 2 || mill_axis == 3 || mill_axis == -1 || mill_axis == -2 || mill_axis == -3);	
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
		
//		if(screen) fprintf(screen ,"\n Exiting function FixBreakparticleForce(). \n");
//		if(logfile) fprintf(logfile ,"\n Exiting function FixBreakparticleForce(). \n");
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
		if(x_sphere) memory->destroy(x_sphere);
		if(r_sphere) memory->sfree(r_sphere);
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
			
		//		if(ECS_flag == 0)
		//			{
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
			//			if(screen) fprintf(screen, "FBPF: density_particle = %f \n",density_particle);
			//			if(logfile) fprintf(logfile, "FBPF: density_particle = %f \n",density_particle);
//						density_particle = class_multisphere->particle_density;
						
		////				if(screen) fprintf(screen ,"\n n_fragments assigned as %d in function calc_insertion_properties. \n",n_fragments);
		////				if(logfile) fprintf(logfile ,"\n n_fragments assigned as %d in function calc_insertion_properties. \n",n_fragments);
		//			}
		
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
//		if(screen) fprintf(screen ,"\n Entering function end_of_step(). \n \n ");
//		if(logfile) fprintf(logfile ,"\n Entering function end_of_step(). \n \n");
		
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
		
		
		double rmass_total = 0.0;
		for(int i = 0; i < nlocal; i++)
			{
				rmass_total += rmass[i];
			}
			
		if(update->dt != previous_timestep)
			{
				if(screen) fprintf(screen, "timestep = %f \n", update->dt);
				if(logfile) fprintf(logfile, "timestep = %f \n", update->dt);	
				previous_timestep = update->dt;	
			}	
		
		double rmass_break = 0.0;
		
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
		if(flag_count != 0)
			{
				if(screen) fprintf(screen, "\n %d flags (mass = %f) out of %d (mass = %f) turned to 1 at timestep = %d \n", flag_count,rmass_break,nlocal,rmass_total,update->ntimestep);
				if(logfile) fprintf(logfile, "\n %d flags (mass = %f) out of %d (mass = %f) turned to 1 at timestep = %d \n", flag_count,rmass_break,nlocal,rmass_total,update->ntimestep);
			}
			
//		if(screen) fprintf(screen ,"\n Exiting function end_of_step(). \n \n");
//		if(logfile) fprintf(logfile ,"\n Exiting function end_of_step(). \n \n");	
	}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::pre_insert()		
	//things to be done before particle insertion//
	{
//		if(screen) fprintf(screen ,"\n Entering function pre_insert(). \n");
//		if(logfile) fprintf(logfile ,"\n Entering function pre_insert(). \n");
		
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
						//total number of particles to be broken//
						n_break_this_local++;
							//total mass of particle to be broken//rmass[i] is mass of particle based on r//
						mass_break_this_local += rmass[i];
							
			////			if(screen) fprintf(screen ,"i = %d (cum_n_break_this_local = %d, cum_mass_break_this_local = %f) \t \n",i,n_break_this_local,mass_break_this_local);	
			////			if(logfile) fprintf(logfile ,"i = %d (cum_n_break_this_local = %d, cum_mass_break_this_local = %f) \t \n",i,n_break_this_local,mass_break_this_local);							
					}
			}
			
		if(n_break_this_local > 0)
			{
				if(screen) fprintf(screen ,"n_break_this_local = %d, mass_break_this_local = %f) \t \n",n_break_this_local,mass_break_this_local);	
				if(logfile) fprintf(logfile ,"n_break_this_local = %d, mass_break_this_local = %f) \t \n",n_break_this_local,mass_break_this_local);
			}
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
	   
	   if(n_break_this > 0)
			{
				if(screen) fprintf(screen ,"n_break_this = %d, mass_break_this = %f, n_break = %d, mass_break = %f. \n \n",n_break_this,mass_break_this,n_break,mass_break);	
				if(logfile) fprintf(logfile ,"n_break_this = %d, mass_break_this = %f, n_break = %d, mass_break = %f. \n \n",n_break_this,mass_break_this,n_break,mass_break);	
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
	   double EI_particle; double ECS_particle;
	   
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
					////			  fprintf(screen, "ECS_particle = %f \n", ECS_particle);
					////			  fprintf(logfile, "ECS_particle = %f \n", ECS_particle);
								  breakdata[ibreak][8] = ECS_particle;
								  
							////	  if(screen) fprintf(screen ,"deleting particles: i = %d, ibreak = %d, Size = %.8g m, EI = %.8g kWh, ECS = %.8g kWh/ton, f = %f N, rmass = %.8g kg,atom_type = %d, density = %f. \n",i,ibreak,size[ibreak],EI_particle,ECS[ibreak],sqrt(vectorMag3DSquared(f[ibreak])),rmass[ibreak],parent_particle_atom_type[ibreak], parent_particle_density[ibreak]  );
							////	  if(logfile) fprintf(logfile ,"deleting particles: i = %d, ibreak = %d, Size = %.8g m, EI = %.8g kWh, ECS = %.8g kWh/ton, f = %f N, rmass = %.8g kg,atom_type = %d, density = %f. \n",i,ibreak,size[ibreak],EI_particle,ECS[ibreak],sqrt(vectorMag3DSquared(f[ibreak])),rmass[ibreak],parent_particle_atom_type[ibreak], parent_particle_density[ibreak] );
																												
					////		  if(screen) fprintf(screen ,"del par: i = %d, ibreak = %d, Size = %f m, density = %f, ECS = %f kWh/ton, f = %f N, rmass = %f kg \n",i,ibreak,2.0*breakdata[ibreak][6],breakdata[ibreak][7], breakdata[ibreak][8],sqrt(vectorMag3DSquared(f[i])),rmass[i]);
					////		  if(logfile) fprintf(logfile ,"del par: i = %d, ibreak = %d, Size = %f m, density = %f, ECS = %f kWh/ton, f = %f N, rmass = %f kg \n",i,ibreak,2.0*breakdata[ibreak][6],breakdata[ibreak][7], breakdata[ibreak][8],sqrt(vectorMag3DSquared(f[i])),rmass[i]);																				
								 
								 
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

	   if (nlocal != nlocal_index)
			{
				if(screen) fprintf(screen ,"\n \n nlocal reduced from %d to %d. \n",nlocal_index,nlocal);	
				if(logfile) fprintf(logfile ,"\n \n nlocal reduced from %d to %d. \n",nlocal_index,nlocal);	
			}
	   // update local and global # particles
	   
	   atom->nlocal = nlocal;	
			//updating no of particles in current processor//After several excutions of nlocal--//Particles that are going to be break have been deleted//
	   double rlocal = static_cast<double>(atom->nlocal);
	   MPI_Allreduce(&rlocal,&atom->natoms,1,MPI_DOUBLE,MPI_SUM,world);   
			//reducing memory in atom->natoms according to new nlocal//

	   // print stats
	   print_stats_breakage_during();     //print stats in terminal//
	   
//	   if(screen) fprintf(screen ,"\n Exiting function pre_insert(). \n");
//	   if(logfile) fprintf(logfile ,"\n Exiting function pre_insert(). \n");
	}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::print_stats_breakage_during()
	{
//	  if(screen) fprintf(screen ,"\n Entering function print_stats_breakage_during(). \n");
//	  if(logfile) fprintf(logfile ,"\n Entering function print_stats_breakage_during(). \n");
	  
	  int step = update->ntimestep;

	  if (me == 0 && n_break_this > 0)
		  {
				if (screen)
				  fprintf(screen ,"Particle breakage: broke %d particles (mass %f) at step %d\n - a total of %d particles (mass %f) broken so far \n \n \n \n \n",
					  n_break_this,mass_break_this,step,n_break,mass_break);	//n_break_this is no of particles broken currently whereas n_break is the total no of particles so far//

				if (logfile)
				  fprintf(logfile,"Particle breakage: broke %d particles (mass %f) at step %d\n - a total of %d particles (mass %f) broken so far \n \n \n \n \n",
					  n_break_this,mass_break_this,step,n_break,mass_break);
		  }
		  
//	  if(screen) fprintf(screen ,"\n Exiting function print_stats_breakage_during(). \n");
//	  if(logfile) fprintf(logfile ,"\n Exiting function print_stats_breakage_during(). \n");  
		  
	}

/* ---------------------------------------------------------------------- */

int FixBreakparticleForce::calc_ninsert_this()		
	{
//		if(screen) fprintf(screen ,"\n Entering function calc_ninsert_this(). \n");
//		if(logfile) fprintf(logfile ,"\n Entering function calc_ninsert_this(). \n");

/****************************************************************************************************************************************/		
		if(ECS_flag == 1)
			{
				ninsert_daughter = 0;
				
				int iparticle = 0;
				int cum_ninsert_daughter = 0;
				
				r_sphere = NULL;
				if(r_sphere) memory->sfree(r_sphere);

				size_r_sphere = 0;
				
				number_per_break = NULL;
				if(number_per_break) memory->destroy(number_per_break);
				number_per_break = new int[n_break_this_local];
				
				while(iparticle < n_break_this_local)
						{
									int cum_ninsert_daughter_iparticle = 0;
								
									size_iparticle = 2.0 * breakdata[iparticle][6];
									ECS_iparticle = breakdata[iparticle][8];		
																							
									alpha = C_GM * size_iparticle + alpha_intercept;
															
									t10 = (JK_parameter_A / 100.0) * (1.0-pow(e,((-JK_parameter_b) * ECS_iparticle)));
															
									double t10_in_percent = t10 * 100.0;	double size_in_mm = 1000.0 * size_iparticle;
															
									if(screen) fprintf(screen,"\n \n \n Size (%f m), ECS (%f kWh/ton), JK_A (%f), JK_b (%f), alpha = %f, t10_in_percent = %f. \n", size_iparticle,ECS_iparticle,JK_parameter_A,JK_parameter_b,alpha,t10_in_percent);
									if(logfile) fprintf(logfile,"\n \n \n Size (%f m), ECS (%f kWh/ton), JK_A (%f), JK_b (%f), alpha = %f, t10_in_percent = %f. \n", size_iparticle,ECS_iparticle,JK_parameter_A,JK_parameter_b,alpha,t10_in_percent);
															
									int ii=0; double Diff_tn=0.0, cum_Diff_tn = 0.0;  //ii is for incrementing sizes of particles//Diff_tn is the differential weightage of material in size class ii//
									int number_particles = 0; //n is the number of particle in i size class//r_n is radius of particle at that size class//
										
									if(screen) fprintf(screen ,"\n Sieves selected \t");
									if(logfile) fprintf(logfile ,"\n Sieves selected \t");
															
									int number_of_sieves_smaller = 0;
									while(sieves[ii] < size_iparticle)
										{
											number_of_sieves_smaller++;
																														
											if(screen) fprintf(screen ,"%f  \t",sieves[ii]);	
											if(logfile) fprintf(logfile ,"%f  \t",sieves[ii]);	
																			
											ii++;
										}
														
						//			if(screen) fprintf(screen ,"\n Number of sieves selected = %d \n",number_of_sieves_smaller);
						//			if(logfile) fprintf(logfile ,"\n Number of sieves selected = %d \n",number_of_sieves_smaller);
															
									int number_particles_per_size_class[number_of_sieves_smaller];
										
									//initializing values
									for(int count=0; count < number_of_sieves_smaller; count++)
										{
												number_particles_per_size_class[count] = 0;
										}
																						
									ii = 0;		
									
									if(screen) fprintf(screen ,"\n size_class \t num_par \t Diff_tn \t  \t \t cum_Diff_tn \t cum_ninsert_daughter_iparticle \t \n");
									if(logfile) fprintf(logfile ,"\n size_class \t num_par \t Diff_tn \t  \t \t cum_Diff_tn \t cum_ninsert_daughter_iparticle \t \n",sieves[ii],number_particles_per_size_class[ii],  Diff_tn,100.0*t10,alpha,cum_Diff_tn,cum_ninsert_daughter_iparticle);
																
									while(sieves[ii] < size_iparticle)
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
														
														assert(Diff_tn >= 0.0);
																																					
													}else 
													{
														error->all(FLERR,"Illegal ECS_flag value");
													}
																				
												cum_Diff_tn += Diff_tn;

												number_particles = (int)round((Diff_tn * pow((size_iparticle / sieves[ii]),3)) + 0.5);	//Differential tn value is multiplies by (R/r)^3 to get the number of particles in that szie. 0.5 is added so that round commnad can work properly. combination of <int>(round) will give nearest integer
																				
												number_particles_per_size_class[ii] = number_particles;
																				
												cum_ninsert_daughter_iparticle += number_particles;
																																
												if(screen) fprintf(screen ," %f \t %d \t \t %f \t  \t \t %f \t \t \t %d \n",sieves[ii], number_particles_per_size_class[ii], Diff_tn,cum_Diff_tn,cum_ninsert_daughter_iparticle);
												if(logfile) fprintf(logfile ," %f \t %d \t \t %f \t  \t \t %f \t \t \t %d \n",sieves[ii], number_particles_per_size_class[ii], Diff_tn,cum_Diff_tn,cum_ninsert_daughter_iparticle);
																																
												ii++;
										}
																														
									
									if(!r_sphere) 
										{
											r_sphere = new double[cum_ninsert_daughter_iparticle];
										}else
										{
											/****************************Resizing the r_sphere array*********************************/
											double *resize_array = new double[size_r_sphere + cum_ninsert_daughter_iparticle];
		
											for(int i = 0; i < size_r_sphere; i++)
												{
													resize_array[i] = r_sphere[i];
												}	
											
											r_sphere = resize_array;
											resize_array = NULL;
											if(resize_array) memory->sfree(resize_array);
											delete[] resize_array;	
																				
											//	resize(cum_ninsert_daughter_iparticle, size_r_sphere, r_sphere);
											/****************************************************************************************/
										}
									
									int jj=0; int kk=0; int mm=0; int count=0; int cum_num_par = 0;
									
									if(screen) fprintf(screen ,"\n iparticle \t size_r_sphere \t \t r_sphere \t daughter_iparticle \t \n");
									if(logfile) fprintf(logfile ,"\n iparticle \t size_r_sphere \t \t r_sphere \t daughter_iparticle \t \n");
																								
									for(jj=0; jj < number_of_sieves_smaller; jj++)
										{
												cum_num_par += number_particles_per_size_class[jj];
																			
												for(kk=0; kk < number_particles_per_size_class[jj]; kk++)
													{																				
															r_sphere[size_r_sphere + count + kk] = sieves[jj] / 2.0;
													}
																			
												count = count + number_particles_per_size_class[jj];
										}
									
									size_r_sphere += cum_ninsert_daughter_iparticle;	
									number_per_break[iparticle] = cum_ninsert_daughter_iparticle;
									cum_ninsert_daughter += cum_ninsert_daughter_iparticle;
										
									if(screen) fprintf(screen ,"\t %d \t \t %d \t  \t \t %f \t \t %d  \n", iparticle, size_r_sphere, r_sphere[size_r_sphere-1], cum_num_par);
									if(logfile) fprintf(logfile ,"\t %d \t \t %d \t  \t \t %f \t \t %d  \n", iparticle, size_r_sphere, r_sphere[size_r_sphere-1], cum_num_par);	
										
									iparticle++;	//while loop increment//
						}
				
/********************************************************************************************************************************************************/
				ninsert_daughter = cum_ninsert_daughter;		
			}
	
//		if(screen) fprintf(screen ,"\n Exiting function calc_ninsert_this(). \n");
//		if(logfile) fprintf(logfile ,"\n Exiting function calc_ninsert_this(). \n");
		
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
//		if(screen) fprintf(screen ,"\n \n \n \n \n \n Entering function x_v_omega() \n");
//		if(logfile) fprintf(logfile ,"\n \n \n \n \n \n Entering function x_v_omega() \n");
		
		double pos_ins[3],v_ins[3],omega_ins[3],quat_ins[4],rad_broken;
		int iparticle, nins;
		ParticleToInsert *pti;		
			//ParticleToInsert is a class declared in particleToInsert.h//

		vectorZeroize3D(omega_ins);	
			//turning all values to zero//
		vectorZeroize4D(quat_ins);
			//turning all values to zero//

		double mass_inserted_this_local = 0.;		// local insertion
		mass_inserted_this = 0.;					// global insertion will be added by mass_inserted_this_local. variable is in function parameter
		
		int ninserted_this_local = 0;				// local insertion
		ninserted_this = 0;							// global insertion will be added by ninserted_this_local. variable is in function parameter
			
		int ninserted_spheres_this_local = 0;		// local insertion
		ninserted_spheres_this = 0;					// global insertion will be added by ninserted_spheres_this_local. variable is in function parameter
				
		iparticle = 0; 

		if(ECS_flag == 0)
				{
//							if(screen) fprintf(screen ,"\n In while loop of function x_v_omega() of ECS_flag = %d \n \n",ECS_flag);
//							if(logfile) fprintf(logfile ,"\n In while loop of function x_v_omega() of ECS_flag = %d \n \n",ECS_flag);
							
							while(iparticle < n_break_this_local)
									{
											vectorCopy3D(&breakdata[iparticle][0],pos_ins);	
													//copies positions from left to right//
											vectorCopy3D(&breakdata[iparticle][3],v_ins);	
													//copies velocities from left to right//
											rad_broken = breakdata[iparticle][6];			
													//copies original particle radius to rad_broken//
											// get pti and scale it down with radius of broken particle
											pti = fix_distribution->pti_list[iparticle];   
													//fixdistribution is a class pointer to FixParticledistributionDiscrete (declared in fix_insert.h, class defined in fix_particledistribution_discrete.h)//
													//pti_list is double array of class ParticleToInsert (defined in fix_particledistribution_discrete.h, class defined in particletoinsert.h)//
													//See in particleToInsert.h where class ParticleToInsert is defined//
													//pti_list contains details of particle properties of one particle, pti is class pointer of particletoinsert//
													//so, particle properties from fix_distribution or indirectly from fix_particledistribution_discrete.h is being transferred to ParticleToInsert class//
																	
											pti->scale_pti(rad_broken);  
													//rad_broken should be a relative values such 0.3 or 0.5//
													//scale_pti is declared in particletoinsert.h//function defined in particletoinsert.cpp//
													//scale_pti(r_scale) multiplies poistions,velocities and r_bound_ins with r_scale, and volume and mass with r_scale^3//
													//These are relative operations and actual values will be set up in next line//
																		
											nins = pti->set_x_v_omega(pos_ins,v_insert,omega_ins,quat_insert); 	//function defined in particletoinsert.cpp//
													//Copies velocities and omega. Returns values nspheres to nins. Add relative values to position//
																	
																	
											// tally stats
											ninserted_spheres_this_local += nins;  
													//Please note it was 0 earlier so now  ninserted_spheres_this_local=nins//
											mass_inserted_this_local += pti->mass_ins;
													//Again mass_inserted_this_local=mass_ins
											ninserted_this_local++;	
													//increment for number of particles inserted////
											iparticle++;
											
//											if(screen) fprintf(screen ,"\n while loop ending here \t \n \n ");		
//											if(logfile) fprintf(logfile ,"\n while loop ending here \t \n \n ");
									}	
					
				
				}else if(ECS_flag == 1)
				{					
						  if(screen) fprintf(screen ,"\n n_break_this_local = %d \n \n \n",n_break_this_local);
						  if(logfile) fprintf(logfile ,"\n n_break_this_local = %d \n \n \n",n_break_this_local);	
						  
						  if(screen) fprintf(screen ,"\n number_per_break[iparticle] \n");
						  if(logfile) fprintf(logfile ,"\n number_per_break[iparticle] \n");
						  					
						  while(iparticle < n_break_this_local)
						    	{
									if(screen) fprintf(screen ,"%d[%d] \t",number_per_break[iparticle],iparticle);
									if(logfile) fprintf(logfile ,"%d[%d] \t",number_per_break[iparticle],iparticle);
									iparticle++;
							    }
						  iparticle = 0;			  
						  if(screen) fprintf(screen ,"\n  \n");
						  if(logfile) fprintf(logfile ,"\n  \n");
						  
					//	  if(screen) fprintf(screen ,"r_sphere[size_r_sphere] \n");
					//	  if(logfile) fprintf(logfile ,"r_sphere[size_r_sphere] \n");
											
					//	  for(int i = 0; i < size_r_sphere; i++)
					//			{
					//				if(screen) fprintf(screen ,"%f[%d] \t",r_sphere[i],i);
					//				if(logfile) fprintf(logfile ,"%f[%d] \t",r_sphere[i],i);
					//			}
									  
					//    if(screen) fprintf(screen ,"\n  \n");
					//	  if(logfile) fprintf(logfile ,"\n  \n");
											
						  int daughter = 0;
						  int previous_sum = 0;
						  					  						  
						  while(iparticle < n_break_this_local)				
								  {
											previous_sum = daughter;
										
											for(daughter; daughter < (previous_sum + number_per_break[iparticle]); daughter++)
													{
															for(int iii=0; iii<3; iii++)
																{
																	pos_ins[iii] = random_coordinates_generator(iii);
																}
															
															if(screen) fprintf(screen ,"\n daughter + 1 = %d, number_per_break[%d] = %d, previous_sum + 1 + number_per_break[%d]) = %d, coordinates(%f,%f,%f) \n",daughter+1,iparticle,number_per_break[iparticle],iparticle,previous_sum + number_per_break[iparticle],pos_ins[0],pos_ins[1],pos_ins[2]);
															if(logfile) fprintf(logfile ,"\n daughter + 1 = %d, number_per_break[%d] = %d, previous_sum + 1 + number_per_break[%d]) = %d, coordinates(%f,%f,%f) \n",daughter+1,iparticle,number_per_break[iparticle],iparticle,previous_sum + number_per_break[iparticle],pos_ins[0],pos_ins[1],pos_ins[2]);
						  					
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
																											
											ninserted_this_local++;	
													//increment for number of particles inserted////
											iparticle++;
											
//											if(screen) fprintf(screen ,"\n while loop ending here \t \n \n ");		
//											if(logfile) fprintf(logfile ,"\n while loop ending here \t \n \n ");	
										
								}																									
																																		
				} else
				{  
						{
								error->all(FLERR,"Illegal ECS_flag value");
						}
			    }												
		
		// tally stats, have to do this since operation is locally on each process
		// as opposed to e.g. FixInsertPack::x_v_omega()

		LAMMPS_NS::MPI_Sum_Scalar(ninserted_spheres_this_local,ninserted_spheres_this,world);	//copies memory//
		LAMMPS_NS::MPI_Sum_Scalar(ninserted_this_local,ninserted_this,world);
		LAMMPS_NS::MPI_Sum_Scalar(mass_inserted_this_local,mass_inserted_this,world);																														
																		
														
				//increment for number of particles already broken////
//		if(screen) fprintf(screen ,"\n ninserted_spheres_this_local = %d, mass_inserted_this_local= %f \t \n ", ninserted_spheres_this_local, mass_inserted_this_local);		
//		if(logfile) fprintf(logfile ,"\n ninserted_spheres_this_local = %d, mass_inserted_this_local= %f \t \n ", ninserted_spheres_this_local, mass_inserted_this_local);	
							
//	    if(screen) fprintf(screen ,"\n Exiting function x_v_omega(). \n ");
//	    if(logfile) fprintf(logfile ,"\n Exiting function x_v_omega(). \n "); 
	  
//	    if(screen) fprintf(screen ,"\n *********************************************************************** \n \n \n \n \n \n \n \n \n \n");
//	    if(logfile) fprintf(logfile ,"\n ***********************************************************************. \n \n \n \n \n \n \n \n \n \n"); 
	}

//not sure if we realy need this function, however it needs to be a function with this name here for the parent class fix-instert
double FixBreakparticleForce::insertion_fraction()
	{
//		if(screen) fprintf(screen ,"\n Entering function insertion_fraction(). \n");
//		if(logfile) fprintf(logfile ,"\n Entering function insertion_fraction(). \n");
		
		if(volumefraction != 0.5)
			{
				if(screen) fprintf(screen ,"volumefraction = %f. \n", volumefraction);		
				if(logfile) fprintf(logfile ,"volumefraction = %f. \n", volumefraction);	
			}
		
		
		
//		if(screen) fprintf(screen ,"\n Exiting function insertion_fraction(). \n");
//		if(logfile) fprintf(logfile ,"\n Exiting function insertion_fraction(). \n");
		
		return volumefraction;
	}
