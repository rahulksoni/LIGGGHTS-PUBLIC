/* ----------------------------------------------------------------------
LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
Transfer Simulations

www.liggghts.com | www.cfdem.com
Christoph Kloss, christoph.kloss@cfdem.com

LIGGGHTS is based on LAMMPS
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov

Copyright (2003) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporatioin, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.

See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(breakparticle/force,FixBreakparticleForce)

#else

#ifndef LMP_FIX_BREAKPARTICLE_FORCE_H
#define LMP_FIX_BREAKPARTICLE_FORCE_H

#include <iostream>
#include "fix_insert.h"						//Inclusion of fix_insert//

using namespace std;

namespace LAMMPS_NS 				//I think modification of namespace LAMMPS_NS//
	{

		class FixBreakparticleForce : public FixInsert 			//inclusion of public part of fix_insert//
			{
				 public:
					
					  FixBreakparticleForce(class LAMMPS *, int, char **);  //Extract values from script command
					  ~FixBreakparticleForce();		
						
									/******************************************************************************************/
									double e = 2.71828;
									double pi = 3.14159265359;
									double JK_parameter_A = 0.0;
									double JK_parameter_b = 0.0;
									double WI = 0.0; 			//Limestone
									double C_GM = 0.0; 
									double C_WI = 0.0; 
									double C_intercept = 0.0;
									double ECS_max = 0.0;	//kwh/ton
									double min_parent_size_to_break = 0.0;	//Particle below this size won't be eligible for breakage									
									double min_daughter_size = 0.0;		
									double max_feed_particle_size = 0.0;	//feed size of particles (if two or more sizes then max of it)
									double con_fac_joule = 0.0; //Energy in joule in coversion factor
									double con_fac_force = 0.0; //force in Newton in conversion factor
									double conversion_factor = 0.0;
									
									double previous_timestep;	
									bool timestep_resetflag = false;
									bool massdis_flag_1 = false; 		//flag that table first column is printed or not in the differential mass distribution file
									bool massdis_flag_2 = false;		//flag that table first column is printed or not in the cumulative mass distribution file
									int massdis_nevery = 0;
									int massdis_previous_ntimestep = 0; 	bool massdis_previous_ntimestep_flag = false;
									int previous_particles_count = 0;	double rmass_total;	double rmass_break = 0.0;
									int particles_count = 0;
									int random = 0;
									/******************************************************************************************/
									
									/***************************************King Model*****************************************/
									double t10;
									//t10 = JK_parameter_A * (1-pow(2.718,((-JK_parameter_b) * ECS)))//
									 
									double alpha;
									//alpha = C_ECS * ECS + C_GM * GM + C_WI_0 * WI_0 + C_intercept //
									//alpha ~= C_GM * GM + C_WI * WI + C_intercept //
									//alpha = C_GM * GM + alpha_intercept //
									 
									//double alpha_intercept;
									double alpha_intercept = 0.0;
									/******************************************************************************************/
									
									/**********************************Geometry parameters*************************************/
									double x_min = 0.0;		
									double x_max = 0.0;
									double y_min = 0.0;
									double y_max = 0.0;
									double z_min = 0.0;
									double z_max = 0.0;
									int mill_axis = 0;						//(x=1, y=2,z=3)
									double bin_x_min = 0.0;					// bin values are the dimension where undersize will be generated
									double bin_x_max = 0.0;
									double bin_y_min = 0.0; 
									double bin_y_max = 0.0;
									double bin_z_min = 0.0;
									double bin_z_max = 0.0;
									bool regenerate_in_bin_switch = false; 						// whether to re-generate smaller (than min_particle size to break) in the bin or not
									bool passing_truncation_switch = false; 		// if true undersize won't generate in the bin to reduce computational cost. Their mass will be noted down on size basis for future reference of results analysis.
									//bool passing_truncation_first_time_switch = true;	//true means the simulation is running for the first time, else at least one simulation restart
									/******************************************************************************************/								
									


					  void post_create();		//FixInsert::post_create(); , char* fixarg[9];
					  void pre_delete();		
		
					  virtual int setmask();																		
					  virtual void init();		//FixInsert::init(); , updating timestep
					  void init_defaults();			//Sets defaults
					  int calc_ninsert_this();		//return n_break_this;
					  
					  virtual void end_of_step();		//Marks particles for breakage
					  
					   // force threshold and number of fragments per particle
					  double ECS_break; //kWh/ton//tresold//
					  
					  double f_break;	
					  
					  double volumefraction;
					  int n_fragments;
					  
					  /***********This two arrays will be used in JK model and developed model for alpha to calculate particle distribution********************/
					  double *ECS;			//An array that will contain the ECS values for different particles//
					  double *size;  		//An array that will contain the diameters of the particles//
					  double *sieves;		//Will contain sieves equal and larger than the minimum prescribed daughter size
					  double *mass_distribution;	//Will contain the size based mass-distribution of particles present in the simulation
					  double *mass_distribution_truncated;	//Will contain the size based mass-distribution of particles truncated (fragments not re-generated)
					 
					  int *number_per_break;	//array containing number of daughters per particle//
					 				
					  double size_iparticle;		//Current size of particle, based on which alpha will be calculated//
					  double ECS_iparticle;		//Current ECS for particle, based on which alpha will be calculated//
					  
					  int size_r_sphere = 0;		//Array containing the radius values of particles to be inserted//
					  
					  double density_particle = 0.0;    //Will retrieve density of particle value from fix_template_multisphere.h//
					  
					  int sieves_series_length = 0;	int index_lower = 0; int index_upper = 0;		//index will be used to defined smallest and biggest particle in the region//
					  /****************************************************************************************************************************************/

				 protected:

					  // inherited functions
					  virtual void calc_insertion_properties(); 		//fix_fragments = static_cast<FixTemplateMultiplespheres*>(fix_distribution->particletemplates()[0]); 	
					  virtual void pre_insert();		//Count particles to break, copy their data and delete particles
					  int is_nearby(int);		
					  void x_v_omega(int ninsert_this,int &ninserted_this, int &ninserted_spheres_this, double &mass_inserted_this);		
											//generate new particles

					  // functions declared in this class
					  inline void generate_random(double *pos, double rad_broken,double rad_insert);		
					  void print_stats_breakage_during();		

					  // per breakage flag
					  class FixPropertyAtom *fix_break;				  
					  
					  // template holding data of the fragments
					  class FixTemplateMultiplespheres *fix_fragments;			

					  class FixCheckTimestepGran *check_timestep;	//will be used to reset timestep flag in fix_check_timestep_gran.h
			//		  class FixCheckTimestepGran *randomize_single(); 

					  // number of spheres in template
					  int nspheres;

					  // coords of each sphere with respect to center of mass
					  double **x_sphere;    //values of positions will be copied from fragment file in these//

					  // radius of each sphere
					  double *r_sphere;		//radius values will be copied from fragment file in these//In ECS case values will be calculated and copied//
			//		  double *volume_sphere;
			//		  double *mass_sphere;
					
					  // stats for breakage
					  int n_break,n_break_this,n_break_this_local;
					  double mass_break,mass_break_this,mass_break_this_local;

					  // data for particles that to be broken
					  double **breakdata;					  

					  //need non-virtual function for parent class
					  double insertion_fraction();		//returns volume fraction
			};

	}

#endif
#endif
