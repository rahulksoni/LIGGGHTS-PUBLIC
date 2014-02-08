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
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.

See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(breakparticle/force,FixBreakparticleForce)

#else

#ifndef LMP_FIX_BREAKPARTICLE_FORCE_H
#define LMP_FIX_BREAKPARTICLE_FORCE_H

#include "fix_insert.h"						//Inclusion of fix_insert//
#include "fix_check_timestep_gran.h"



namespace LAMMPS_NS 				//I think modification of namespace LAMMPS_NS//
	{

		class FixBreakparticleForce : public FixInsert 			//inclusion of public part of fix_insert//
			{
				 public:

											
											  
					  FixBreakparticleForce(class LAMMPS *, int, char **);  //Extract values from script command
					  ~FixBreakparticleForce();		
						
					  double previous_timestep;	
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
									double *sieves;
									double con_fac_joule = 0.0; //Energy in joule in coversion factor
									double con_fac_force = 0.0; //force in Newton in conversion factor
									double conversion_factor = 0.0;
									
									bool timestep_resetflag = false;
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
					 
					  int *number_per_break;	//array containing number of daughters per particle//
					 				
					  double size_iparticle;		//Current size of particle, based on which alpha will be calculated//
					  double ECS_iparticle;		//Current ECS for particle, based on which alpha will be calculated//
					  
					  int size_r_sphere = 0;
					  
					  double density_particle = 0.0;    //Will retrieve value from fix_template_multisphere.h
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
