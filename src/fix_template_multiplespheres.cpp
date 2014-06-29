/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
Thanks to Chris Stoltz (P&G) for providing
a Fortran version of the MC integrator
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_template_multiplespheres.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "vector_liggghts.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "output.h"
#include "memory.h"
#include "error.h"
#include "random_mars.h"
#include "random_park.h"
#include "fix_rigid.h"
#include "particleToInsert.h"
#include "input_multisphere.h"


#include "assert.h"
/***************************************************************/

using namespace LAMMPS_NS;
using namespace LMP_PROBABILITY_NS;

#define LARGE 1e8
#define EPSILON 1.0e-7
#define N_SHUFFLE_BOUND 200

/* ---------------------------------------------------------------------- */

FixTemplateMultiplespheres::FixTemplateMultiplespheres(LAMMPS *lmp, int narg, char **arg) :
  FixTemplateSphere(lmp, narg, arg)
	{
		//	if(screen) fprintf(screen ,"\n ===>>> TMS:  FixTemplateMultiplespheres(LAMMPS *lmp, int narg, char **arg) . \n");
		//	if(logfile) fprintf(logfile ,"\n ===>>> TMS:  function FixTemplateMultiplespheres(LAMMPS *lmp, int narg, char **arg) . \n");
			
			  if(pdf_density->rand_style() != RANDOM_CONSTANT) error->all(FLERR,"Fix particletemplate/multiplespheres currently only supports constant density");
			  if(pdf_radius) error->fix_error(FLERR,this,"currently does not support keyword 'radius'");
			  if(domain->dimension != 3) error->fix_error(FLERR,this,"only supports 3D simulations");

			  particle_density = expectancy(pdf_density);   //Will be used in fix_breakparticle_force.cpp for marking the particles for breakage
			//  if(screen) fprintf(screen, "FTM: particle_density = %f \n",particle_density);
			//  if(logfile) fprintf(logfile, "FTM: particle_density = %f \n",particle_density);

			  // parse number of spheres
			  if (strcmp(arg[iarg++],"nspheres") != 0) error->fix_error(FLERR,this,"expecting argument 'nspheres'");
			  
			  /**********Important: In our case nspheres need to be controlled**********/
			  nspheres = atoi(arg[iarg++]);
	
			  if (narg < iarg+4) error->fix_error(FLERR,this,"not enough arguments");
			  if (strcmp(arg[iarg++],"ntry") != 0) error->fix_error(FLERR,this,"need 'ntry' to be defined");
			  ntry = static_cast<int>(atoi(arg[iarg++]));

			  memory->create(x_sphere,nspheres,3,"FixTemplateMultiplespheres:x_sphere");
						  
			  ///////////////radius array for no of spheres////////////////
			  r_sphere = new double[nspheres];

			  // re-create pti with correct nspheres
			  delete pti;
			  pti = new ParticleToInsert(lmp,nspheres);
			
			  for (int i = 0; i < 3; i++) 
					{
						x_min[i] = LARGE;
						x_max[i] = -LARGE;
					}

			  bool spheres_read = false;

			  bool hasargs = true;
			  while (iarg < narg && hasargs)
				  {
						hasargs = false;

						if (strcmp(arg[iarg],"spheres") == 0)
							{
								  hasargs = true;
								  spheres_read = true;
								  iarg++;

								  if (strcmp(arg[iarg],"file") == 0)		//file to read data for fragments//
									  {
										  iarg++;
										  if (narg < iarg+3) error->fix_error(FLERR,this,"not enough arguments");

										  char *clmp_filename = arg[iarg++];

										  if (strcmp(arg[iarg++],"scale") != 0) error->fix_error(FLERR,this,"you have to specify a scale factor");
										  scale_fact = atof(arg[iarg++]);	
										  if(scale_fact<=0.) error->fix_error(FLERR,this,"scale factor must be >0");

										  // allocate input class, try to open file, read data from file
										  InputMultisphere *myclmp_input = new InputMultisphere(lmp,0,NULL);
																  
										  /*********************Important: This part is responsible for grabbing data from the fragment file*******************/
										  myclmp_input->clmpfile(clmp_filename,x_sphere,r_sphere,nspheres);
																		/*void clmpfile(const char *,double **,double*,int); of input_multisphere.h is executed
																		  If we see the function definition in input_multisphere.cpp then it copies the content of clmp_filename (actually our fragment file)
																		  to x_sphere and r_sphere for nspheres no of times. positions are copied in x_sphere (in 3D) and radius values are copied in r_sphere.
																		  x_sphere is a 2D array and r_sphere is a 1D array defined in fix_template_multispheres.h 
																		*/
										  /********************************************************************************************************************/
																														
										  delete myclmp_input;
												
										  if(nspheres == 1)
												for(int i = 0; i < 3; i++)
													{
														x_sphere[0][i] = 0.0;
													}
													  
										  for(int i = 0; i < nspheres; i++)
												{
													  if(r_sphere[i] <= 0.) error->fix_error(FLERR,this,"radius must be > 0");
													  r_sphere[i] *= (scale_fact*force->cg());		
															//cg is coarsegraining which possess value 1.0 in force.cpp//
															//multiplying radius values with scale value which is 1.0 in breakage example//
													  vectorScalarMult3D(x_sphere[i],scale_fact*force->cg());
															//multiplying coordinates with scale value//
												}

										  // calculate bounding box
										  /*************If particle is out of already created bounding then boundingg box is being modified************/
										  for(int i = 0; i < nspheres; i++)
											    {
													  for(int j=0;j<3;j++)
														  {
																if (x_sphere[i][j]-r_sphere[i]<x_min[j]) x_min[j] = x_sphere[i][j]-r_sphere[i];
																if (x_sphere[i][j]+r_sphere[i]>x_max[j]) x_max[j] = x_sphere[i][j]+r_sphere[i];
														  }
												}
										   /************************************************************************************************************/		  
									  }
								  else 			//Not of our interest, just an error check//
									  {
											if (narg < iarg + 4*nspheres) error->fix_error(FLERR,this,"not enough arguments");

											  //read sphere r and coos, determine min and max
											for(int i = 0; i < nspheres; i++)
												{
													  r_sphere[i] = atof(arg[iarg+3])*force->cg();
													  if(r_sphere[i] <= 0.) error->fix_error(FLERR,this,"radius must be >0");
													  for(int j = 0; j < 3; j++)
															{
																	x_sphere[i][j] = atof(arg[iarg+j])*force->cg();
																	if (x_sphere[i][j]-r_sphere[i]<x_min[j]) x_min[j]=x_sphere[i][j]-r_sphere[i];
																	if (x_sphere[i][j]+r_sphere[i]>x_max[j]) x_max[j]=x_sphere[i][j]+r_sphere[i];
															}
											          iarg+=4;
												}
									  }
							}
					   else if(strcmp(style,"particletemplate/multiplespheres") == 0)		//Not of interest, just an error check//
							error->fix_error(FLERR,this,"unknown keyword");
				  }

				  if(!spheres_read) error->fix_error(FLERR,this,"need to define spheres for the template");		//error check only; if sphere properties reading was unsuccessful//

				  if(comm->me == 0 && screen) fprintf(screen,"Calculating the properties of the given template.\n   Depending on ntry, this may take a while...\n");

				  if(ntry < 1e3) error->fix_error(FLERR,this,"ntry is too low");
				  
				  if(comm->me == 0 && ntry < 1e5) error->warning(FLERR,"fix particletemplate/multisphere: ntry is very low");
				  
	//	if(screen) fprintf(screen ,"\n <<<=== TMS:  FixTemplateMultiplespheres(LAMMPS *lmp, int narg, char **arg) . \n");
	//	if(logfile) fprintf(logfile ,"\n <<<=== TMS:  function FixTemplateMultiplespheres(LAMMPS *lmp, int narg, char **arg) . \n");

	}

/* ---------------------------------------------------------------------- */

FixTemplateMultiplespheres::~FixTemplateMultiplespheres()
	{
	//	if(screen) fprintf(screen ,"\n ===>>> TMS:  ~FixTemplateMultiplespheres(). \n");
	//	if(logfile) fprintf(logfile ,"\n ===>>> TMS:  ~FixTemplateMultiplespheres(). \n");
		
		memory->destroy(x_sphere);
		delete []r_sphere;
	//	delete []mass_sphere;
	//	if(screen) fprintf(screen ,"\n <<<=== TMS:   ~FixTemplateMultiplespheres(). \n");
	//	if(logfile) fprintf(logfile ,"\n <<<=== TMS:   ~FixTemplateMultiplespheres(). \n");
	
	}

/* ---------------------------------------------------------------------- */

void FixTemplateMultiplespheres::post_create()
	{
	//	if(screen) fprintf(screen ,"\n ===>>> TMS:  post_create(). \n");
	//	if(logfile) fprintf(logfile ,"\n ===>>> TMS:  post_create(). \n");
		// calculate bounding sphere and center of mass
		// also transforms sphere coordinates so that com = 0/0/0
		
			calc_bounding_sphere();
			calc_center_of_mass();
			
	//	if(screen) fprintf(screen ,"\n <<<=== TMS:   post_create(). \n");
	//	if(logfile) fprintf(logfile ,"\n <<<=== TMS:   post_create(). \n");	
	
	}


/* ----------------------------------------------------------------------
   calc bounding sphere with iterative procedure

   do this multiple times, randomizing point order every time, see
   http://hacksoflife.blogspot.com/2009/01/randomized-bounding-spheres.html
   choose optimal result at the end - this gives linear run-time
------------------------------------------------------------------------- */

void FixTemplateMultiplespheres::calc_bounding_sphere()
	{
	//	if(screen) fprintf(screen ,"\n ===>>> TMS:  calc_bounding_sphere(). \n");
	//	if(logfile) fprintf(logfile ,"\n ===>>> TMS:  calc_bounding_sphere(). \n");
		
		  r_bound = LARGE;			//A set limit for the r_bound
		  int *visited = new int[nspheres];		//An empty array of the nsphere length
		  double d[3],delta,dist;		//d is an empty vector of length 3

		  for(int shuffle = 0; shuffle < N_SHUFFLE_BOUND; shuffle ++)		//N_SHUFFLE_BOUND = 200
		  {
			  for(int i = 0; i < nspheres; i++) visited[i] = 0;			//Filling all elements of visited with 0

			  int isphere = -1;			//initializing with negative value that will be able to start upcoming loop
			  int nvisited = 0;			//initializing nvisited. it is a counter to count number of particles that have been considered
			  double x_bound_temp[3],rbound_temp;		//An empty vector of length 3. temporary radius

			  while(isphere < 0 || visited[isphere] || isphere >= nspheres )
				  isphere = static_cast<int>(random->uniform()*nspheres);		//Randomly choosing a value of 0 < isphere < nspheres

			  nvisited++;		//increment since 1 particle visited as isphere
			  visited[isphere] = 1;		//

			  vectorCopy3D(x_sphere[isphere],x_bound_temp);
			  rbound_temp = r_sphere[isphere];

			  while(nvisited < nspheres)
				  {
					  while(isphere < 0 || visited[isphere] || isphere >= nspheres )
						   isphere = static_cast<int>(random->uniform()*nspheres);

					  nvisited++;
					  visited[isphere] = 1;

					  vectorSubtract3D(x_sphere[isphere],x_bound_temp,d);
					  dist = vectorMag3D(d);

					  // do nothing if sphere is completely contained in bounding sphere
					  // if not contained, shift and extend bounding sphere
					  if(dist + r_sphere[isphere] > rbound_temp)
						  {
							  double fact = (dist + r_sphere[isphere] - rbound_temp) / (2. * dist);
							  vectorScalarMult3D(d,fact);
							  vectorAdd3D(x_bound_temp,d,x_bound_temp);
							  rbound_temp += vectorMag3D(d);
						  }
					  
				  }
			  if(rbound_temp < r_bound)
				  {
					  r_bound = rbound_temp;
					  vectorCopy3D(x_bound_temp,x_bound);
				  }
			  
		  }
		  delete []visited;

		  // do a coarse check on the validity of the bounding sphere calc
		  for(int i = 0; i < nspheres; i++)
			  {
				  double temp[3];
				  vectorSubtract3D(x_bound,x_sphere[i],temp);
				  if(vectorMag3D(temp) > r_bound) error->fix_error(FLERR,this,"Bounding sphere calculation for template failed");
			  }
			  
	//	if(screen) fprintf(screen ,"\n <<<=== TMS:   calc_bounding_sphere(). \n");
	//	if(logfile) fprintf(logfile ,"\n <<<=== TMS:   calc_bounding_sphere(). \n");	  

	}


/* ----------------------------------------------------------------------
   calc center of mass
------------------------------------------------------------------------- */

void FixTemplateMultiplespheres::calc_center_of_mass()		//For particles of same sizes it calculates center of mass otherwise it is centroid
	{
		
	//	if(screen) fprintf(screen ,"\n ===>>> TMS:  calc_center_of_mass(). \n");
	//	if(logfile) fprintf(logfile ,"\n ===>>> TMS:  calc_center_of_mass(). \n");
		  
	//	if(ECS_flag_break_h == 0)
	//		{
				  // mc integration, calc volume and com, mass weight
				  int nsuccess = 0;

				  double x_try[3],xcm[3],dist_j_sqr;
				  int n_found = 0;

				  vectorZeroize3D(xcm);

				  bool alreadyChecked = false;
				  for(int i = 0; i < ntry; i++)
					  {
						  generate_xtry(x_try);

						  alreadyChecked = false;
						  for(int j = 0; j < nspheres; j++)
							  {
								  dist_j_sqr = dist_sqr(j,x_try);

								  // only count once if contained in multiple spheres
								  if (alreadyChecked) break;
								  if(dist_j_sqr < r_sphere[j]*r_sphere[j])
									  {
										  xcm[0] = (xcm[0]*static_cast<double>(nsuccess)+x_try[0])/static_cast<double>(nsuccess+1);
										  xcm[1] = (xcm[1]*static_cast<double>(nsuccess)+x_try[1])/static_cast<double>(nsuccess+1);
										  xcm[2] = (xcm[2]*static_cast<double>(nsuccess)+x_try[2])/static_cast<double>(nsuccess+1);
										  nsuccess++;
										  alreadyChecked = true;
									  }
							  }
					  }

				  // expectancy values
				  volume_expect = static_cast<double>(nsuccess)/static_cast<double>(ntry)*(x_max[0]-x_min[0])*(x_max[1]-x_min[1])*(x_max[2]-x_min[2]);
				  mass_expect = volume_expect*expectancy(pdf_density);
				  
				 
				  
				  r_equiv = pow(6.*mass_expect/(8.*expectancy(pdf_density)*M_PI),1./3.);

				  // transform into a system with center of mass=0/0/0
				  
				 
						  for(int i = 0; i < nspheres; i++)
						  vectorSubtract3D(x_sphere[i],xcm,x_sphere[i]);
					

				  vectorSubtract3D(x_min,xcm,x_min);
				  vectorSubtract3D(x_max,xcm,x_max);
				  vectorSubtract3D(x_bound,xcm,x_bound);
		//	}
		
	//	if(screen) fprintf(screen ,"\n <<<=== TMS:   calc_center_of_mass(). \n");
	//	if(logfile) fprintf(logfile ,"\n <<<=== TMS:   calc_center_of_mass(). \n");
	}
	
	
/* ----------------------------------------------------------------------*/

double FixTemplateMultiplespheres::max_r_bound()
	{
	//	if(screen) fprintf(screen ,"\n ===>>> TMS:  max_r_bound(). \n");
	//	if(logfile) fprintf(logfile ,"\n ===>>> TMS:  max_r_bound(). \n");
		
	//	if(screen) fprintf(screen ,"\n <<<=== TMS:   max_r_bound(). \n");
	//	if(logfile) fprintf(logfile ,"\n <<<=== TMS:   max_r_bound(). \n");
		
		return r_bound;
	}
	
/* ----------------------------------------------------------------------*/

double FixTemplateMultiplespheres::max_rad()
	{
		
	//	if(screen) fprintf(screen ,"\n ===>>> TMS:  max_rad(). \n");
	//	if(logfile) fprintf(logfile ,"\n ===>>> TMS:  max_rad(). \n");
		
		double rmax = 0.;

		for(int j = 0; j < nspheres; j++)
		  if(rmax < r_sphere[j]) rmax = r_sphere[j];

	//	if(screen) fprintf(screen ,"\n <<<=== TMS:   max_rad(). \n");
	//	if(logfile) fprintf(logfile ,"\n <<<=== TMS:   max_rad(). \n");

		return rmax;
	}


/* ----------------------------------------------------------------------*/

double FixTemplateMultiplespheres::min_rad()
	{
		
	//	if(screen) fprintf(screen ,"\n ===>>> TMS:  min_rad(). \n");
	//	if(logfile) fprintf(logfile ,"\n ===>>> TMS:  min_rad(). \n");
		
		double rmin = 0.;

		for(int j = 0; j < nspheres; j++)
		  if(rmin > r_sphere[j]) rmin = r_sphere[j];

	//	if(screen) fprintf(screen ,"\n <<<=== TMS:   min_rad(). \n");
	//	if(logfile) fprintf(logfile ,"\n <<<=== TMS:   min_rad(). \n");

		return rmin;
	}

/* ----------------------------------------------------------------------*/

int FixTemplateMultiplespheres::number_spheres()
	{
	//	if(screen) fprintf(screen ,"\n ===>>> TMS: number_spheres(). \n");
	//	if(logfile) fprintf(logfile ,"\n ===>>> TMS: number_spheres(). \n");
		
	//	if(screen) fprintf(screen ,"\n <<<=== TMS:  number_spheres(). \n");
	//	if(logfile) fprintf(logfile ,"\n <<<=== TMS:  number_spheres(). \n");
		
		return nspheres;
	}
	
	

/* ----------------------------------------------------------------------*/

void FixTemplateMultiplespheres::randomize_single()
	{
	//	if(screen) fprintf(screen ,"\n ===>>> TMS: randomize_single(). \n");
	//	if(logfile) fprintf(logfile ,"\n ===>>> TMS: randomize_single(). \n");
	  
		  pti->nspheres = nspheres;
		  pti->density_ins = expectancy(pdf_density);
	//	  if(ECS_flag_break_h == 0)
	//		{
				pti->volume_ins = volume_expect;
				pti->mass_ins = mass_expect;
	//		}else
	//		{
	//			pti->volume_ins = volume_particles_inserted;
	//			pti->mass_ins = mass_particles_inserted;
	//		}
		  pti->r_bound_ins = r_bound;
		  pti->atom_type = atom_type;

		  for(int j = 0; j < nspheres; j++)
			  {
				  pti->radius_ins[j] = r_sphere[j];
				  vectorCopy3D(x_sphere[j],pti->x_ins[j]);
			  }

		  vectorZeroize3D(pti->v_ins);
		  vectorZeroize3D(pti->omega_ins);

		  pti->groupbit = groupbit;
		  
	//    if(screen) fprintf(screen ,"\n <<<=== TMS: randomize_single(). \n");
	//	if(logfile) fprintf(logfile ,"\n <<<=== TMS: randomize_single(). \n");
	}

/* ----------------------------------------------------------------------*/

void FixTemplateMultiplespheres::init_ptilist(int n_random_max)
	{
		if(screen) fprintf(screen ,"\n ===>>> TMS: init_ptilist(int n_random_max). \n");
		if(logfile) fprintf(logfile ,"\n ===>>> TMS: init_ptilist(int n_random_max). \n");

		if(screen) fprintf(screen ,"\n TMS: n_random_max = %d \n", n_random_max);
		if(logfile) fprintf(logfile ,"\n TMS: n_random_max = %d \n", n_random_max);

		if(screen) fprintf(screen ,"\n TMS: sizeof(ParticleToInsert*) = %d \n", sizeof(ParticleToInsert*));
		if(logfile) fprintf(logfile ,"\n TMS: sizeof(ParticleToInsert*) = %d \n", sizeof(ParticleToInsert*));
		
		if(pti_list) error->all(FLERR,"invalid FixTemplateSphere::init_list()");
		n_pti_max = n_random_max;
		pti_list = (ParticleToInsert**) memory->smalloc(n_pti_max*sizeof(ParticleToInsert*),"pti_list");
		for(int i = 0; i < n_pti_max; i++)
		   pti_list[i] = new ParticleToInsert(lmp,nspheres);

		if(screen) fprintf(screen ,"\n TMS: nspheres = %d \n", nspheres);
		if(logfile) fprintf(logfile ,"\n TMS: nspheres = %d \n", nspheres);
		   
		if(screen) fprintf(screen ,"\n <<<=== TMS: init_ptilist(int n_random_max). \n");
		if(logfile) fprintf(logfile ,"\n <<<=== TMS: init_ptilist(int n_random_max). \n");   
	}


/* ----------------------------------------------------------------------*/

void FixTemplateMultiplespheres::randomize_ptilist(int n_random,int distribution_groupbit)
	{
		
	if(screen) fprintf(screen ,"\n ===>>> TMS: randomize_ptilist(int n_random,int distribution_groupbit). \n");
	if(logfile) fprintf(logfile ,"\n ===>>> TMS: randomize_ptilist(int n_random,int distribution_groupbit). \n");
		
		if(screen) fprintf(screen ,"\n TMS: n_random = %d \n", n_random);
		if(logfile) fprintf(logfile ,"\n TMS: n_random = %d \n", n_random);

		for(int i = 0; i < n_random; i++)
			{
	//			if(screen) fprintf(screen ,"\n ===>>> TMS: randomize_ptilist(int n_random,int distribution_groupbit): 1 n_random = %d, i = %d \n",n_random,i);
	//			if(logfile) fprintf(logfile ,"\n ===>>> TMS: randomize_ptilist(int n_random,int distribution_groupbit): 1 n_random = %d, i = %d \n",n_random,i);
				ParticleToInsert *pti = pti_list[i];
						/* An example that how can pointers be defined:
							  	Pointers can be initialized either to the address of a variable (such as in the case above), or to the value of another pointer (or array):
										1		int myvar;
										2		int *foo = &myvar;
										3		int *bar = foo;
							*/

						/*This function updates the the values in pti_list[i] for all i values*/

				pti->density_ins = expectancy(pdf_density);
			
		//		  if(!ECS_flag_break_h)
		//			{
						pti->volume_ins = volume_expect;
				
						pti->mass_ins = mass_expect;
				
		//			}else
		//			{
		//				pti->volume_ins = volume_particles_inserted;
		//				pti->mass_ins = mass_particles_inserted;
		//			}
				  pti->r_bound_ins = r_bound;
			
				  pti->atom_type = atom_type;
				

				  for(int j = 0; j < nspheres; j++)
					  {
			
						  pti->radius_ins[j] = r_sphere[j];
	
						  vectorCopy3D(x_sphere[j],pti->x_ins[j]);
						
					  }

				  vectorZeroize3D(pti->v_ins);
				
				  vectorZeroize3D(pti->omega_ins);
			

				  pti->groupbit = groupbit | distribution_groupbit; 
			
			}
			
			
		if(screen) fprintf(screen ,"\n <<<=== TMS: randomize_ptilist(int n_random,int distribution_groupbit). \n");
		if(logfile) fprintf(logfile ,"\n <<<=== TMS: randomize_ptilist(int n_random,int distribution_groupbit). \n");  
	}



/********************Protected parts of class******************************/



	
/* ----------------------------------------------------------------------
   sqr distance from x_sphere[j] to xtest
------------------------------------------------------------------------- */

double FixTemplateMultiplespheres::dist_sqr(int j, double *xtest)
	{
	//	if(screen) fprintf(screen ,"\n ===>>> TMS: dist_sqr(int j, double *xtest). \n");
	//	if(logfile) fprintf(logfile ,"\n ===>>> TMS: dist_sqr(int j, double *xtest). \n");
		
		double dSqr = 0.;
		dSqr += (xtest[0]-x_sphere[j][0])*(xtest[0]-x_sphere[j][0]);
		dSqr += (xtest[1]-x_sphere[j][1])*(xtest[1]-x_sphere[j][1]);
		dSqr += (xtest[2]-x_sphere[j][2])*(xtest[2]-x_sphere[j][2]);
		
	//	if(screen) fprintf(screen ,"\n <<<=== TMS: dist_sqr(int j, double *xtest). \n");
	//	if(logfile) fprintf(logfile ,"\n <<<=== TMS: dist_sqr(int j, double *xtest). \n");  
		return dSqr;
	}

/* ----------------------------------------------------------------------
   generate random point in bbox
------------------------------------------------------------------------- */

void FixTemplateMultiplespheres::generate_xtry(double *x_try)
	{
	//	if(screen) fprintf(screen ,"\n ===>>> TMS: generate_xtry(double *x_try). \n");
	//	if(logfile) fprintf(logfile ,"\n ===>>> TMS: generate_xtry(double *x_try). \n");
		
		x_try[0] = x_min[0]+(x_max[0]-x_min[0])*random->uniform();
		x_try[1] = x_min[1]+(x_max[1]-x_min[1])*random->uniform();
		x_try[2] = x_min[2]+(x_max[2]-x_min[2])*random->uniform();
		
	//	if(screen) fprintf(screen ,"\n <<<=== TMS: generate_xtry(double *x_try). \n");
	//	if(logfile) fprintf(logfile ,"\n <<<=== TMS: generate_xtry(double *x_try). \n");  
	}




