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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_template_sphere.h"
#include "fix_particledistribution_discrete.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "region.h"
#include "update.h"
#include "modify.h"
#include "output.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"
#include "particleToInsert.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define MIN(A,B) (((A) < (B)) ? (A) : (B))
#define MAX(A,B) (((A) > (B)) ? (A) : (B))

#define LMP_DEBUGMODE_SPHERE false

/* ---------------------------------------------------------------------- */

FixParticledistributionDiscrete::FixParticledistributionDiscrete(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
	{
//		  if(screen) fprintf(screen ,"\n ===>>> PDD:  function FixParticledistributionDiscrete(LAMMPS *lmp, int narg, char **arg). \n");
//		  if(logfile) fprintf(logfile ,"\n ===>>> PDD:  function FixParticledistributionDiscrete(LAMMPS *lmp, int narg, char **arg). \n");
		
		  restart_global = 1;

		  // random number generator, same for all procs

		  if (narg < 7)
				error->all(FLERR,"Illegal fix particledistribution/discrete command, not enough arguments");
		  seed = atoi(arg[3]) + comm->me;
		  random = new RanPark(lmp,seed);
		  ntemplates = atoi(arg[4]);
		  if(ntemplates < 1)
				error->all(FLERR,"Illegal fix particledistribution/discrete command, illegal number of templates");

		  templates = new FixTemplateSphere*[ntemplates];
		  distweight = new double[ntemplates];			//distribution based weigthage for each template
		  cumweight = new double[ntemplates];		//Cumulative wirght distribution
		  parttogen = new int[ntemplates];		
		  distorder = new int[ntemplates];		//distribution reorganized in order in ascending or descending order

		  iarg = 5;

		  int itemp=0;

		  if(narg != iarg+2*ntemplates)			//2 values template id and its weightage for each template
				error->all(FLERR,"Illegal fix particledistribution/discrete command, # of templates does not match # of arguments");

		  // parse further args
		  do 
			  {
					if(itemp == ntemplates) break;
					if(narg < iarg+1)
						error->all(FLERR,"Illegal fix particledistribution/discrete command, not enough arguments");
					int ifix = modify->find_fix(arg[iarg]);

					if(ifix < 0)
						error->all(FLERR,"Illegal fix particledistribution/discrete command, invalid ID for fix particletemplate provided");

					if(strncmp(modify->fix[ifix]->style,"particletemplate/",16))
						error->all(FLERR,"Illegal fix particledistribution/discrete command, fix is not of type particletemplate");

					templates[itemp] = static_cast<FixTemplateSphere*>(modify->fix[ifix]);
					distweight[itemp] = atof(arg[iarg+1]);
					if (distweight[itemp] < 0) error->all(FLERR,"Illegal fix particledistribution/discrete command, invalid weight");
					itemp++;
					iarg += 2;
			  } while (iarg < narg);

		  // check for double use of template which is not allowed
		  for(int i = 0; i < ntemplates; i++)
			  for(int j = 0; j < i; j++)
					if(templates[i] == templates[j])
							error->all(FLERR,"Illegal fix particledistribution/discrete command, cannot use the same template twice");

		  // normalize distribution
		  double weightsum = 0;
		  for(int i = 0; i < ntemplates; i++)
				weightsum += distweight[i];		//calculating total weightsum

		  if(comm->me == 0 && fabs(weightsum-1.) > 0.00001)		//total wieghtsum should be equal to 1
				error->warning(FLERR,"particledistribution/discrete: sum of distribution weights != 1, normalizing distribution");

		  for(int i = 0; i < ntemplates; i++)		
				distweight[i]/=weightsum;

		  if(comm->me == 0 && screen)
			  {
				  fprintf(screen,"Fix particledistribution/discrete (id %s): distribution based on mass%%:\n",this->id);
				  for(int i = 0; i < ntemplates; i++)
						fprintf(screen,"    %s: d=%e (max. bounding sphere) mass%%=%f%%\n",templates[i]->id,2.*templates[i]->max_r_bound(),100.*distweight[i]);
			  }

		  // convert distribution from mass% to number%  ///////////////////////////////////////////////////////////
		  for(int i=0;i<ntemplates; i++)
				distweight[i]=distweight[i]/templates[i]->massexpect();

		  weightsum=0;
		  for(int i=0;i<ntemplates; i++)
				weightsum+=distweight[i];

		  for(int i=0;i<ntemplates; i++)
				distweight[i]/=weightsum;

		  if(comm->me == 0 && screen)
			  {
				  fprintf(screen,"Fix particledistribution/discrete (id %s): distribution based on number%%:\n",this->id);
				  for(int i = 0; i < ntemplates; i++)
						fprintf(screen,"    %s: d=%e (max. bounding sphere) number%%=%f%%\n",templates[i]->id,2.*templates[i]->max_r_bound(),100.*distweight[i]);
			  }
			
		  cumweight[0] = distweight[0];
		  for(int i = 1; i < ntemplates; i++)
				cumweight[i] = distweight[i]+cumweight[i-1];

		  volexpect=0.; massexpect=0.;

		  for(int i = 0; i < ntemplates; i++)
			  {
				  volexpect  += templates[i]->volexpect()  * distweight[i];		//getting the desired volume according to its weightage
				  massexpect += templates[i]->massexpect() * distweight[i];		//getting the desired mass according to its weightage
			  }

		  //get min/maxtype
		  maxtype = 0;
		  mintype = 10000;
		  for(int i = 0; i < ntemplates; i++)
			  {
					if(templates[i]->type() > maxtype)		//max number of atom_type (defined in fix_template_sphere.cpp)
						maxtype = templates[i]->type();
					if(templates[i]->type() < mintype)		//min number of atom_type
						mintype = templates[i]->type();
			  }

		  // check which template has the most spheres /////////////////////////
		  maxnspheres = 0;
		  for(int i = 0; i < ntemplates;i++)
				if(templates[i]->number_spheres() > maxnspheres)
						maxnspheres=templates[i]->number_spheres();

		  // sort the distributions by insertion volume (in descending order)
		  // use bubble sort
		  for(int i = 0; i < ntemplates; i++)
				distorder[i]=i;

		  bool swaped;
		  int n = ntemplates;
		  do
			  {
				  swaped = false;
				  for(int i = 0; i < ntemplates-1; i++)
					  {
						  if(templates[distorder[i]]->volexpect() < templates[distorder[i+1]]->volexpect())
							  {
									//swap
									int tmp = distorder[i];
									distorder[i] = distorder[i+1];
									distorder[i+1] = tmp;
									swaped = true;
							  }
					  }
				  n--;
			  } while(swaped && n > 0);

		  pti = templates[distorder[0]]->pti;		//getting pti from template

		  pti_list = NULL;    //Nulified
		  n_pti = n_pti_max = 0;	

		  //calc max radius and bounding sphere radius

		  maxrad = maxrbound = 0.;
		  minrad = 1000.;

		  for(int i = 0; i < ntemplates;i++)		//getting maxrbound
				if(templates[i]->max_r_bound() > maxrbound)
					maxrbound = templates[i]->max_r_bound();

		  for(int i = 0; i < ntemplates;i++)		//getting maxrad
				if(templates[i]->max_rad() > maxrad)
					maxrad = templates[i]->max_rad();

		  for(int i = 0; i < ntemplates;i++)		//getting minrad
				if(templates[i]->min_rad() < minrad)
					minrad = templates[i]->min_rad();
					
//		  if(screen) fprintf(screen ,"\n <<<=== PDD:  function FixParticledistributionDiscrete(LAMMPS *lmp, int narg, char **arg). \n");
//		  if(logfile) fprintf(logfile ,"\n <<<=== PDD:  function FixParticledistributionDiscrete(LAMMPS *lmp, int narg, char **arg). \n");

	}

/* ---------------------------------------------------------------------- */

FixParticledistributionDiscrete::~FixParticledistributionDiscrete()
	{
//		if(screen) fprintf(screen ,"\n ===>>> PDD:  function ~FixParticledistributionDiscrete() \n");
//		if(logfile) fprintf(logfile ,"\n ===>>> PDD:  function ~FixParticledistributionDiscrete() \n");
		
		delete []templates;
		delete []distweight;
		delete []cumweight;
		delete []parttogen;
		delete []distorder;
		if(pti_list) delete []pti_list;		
		delete random;
		
//		if(screen) fprintf(screen ,"\n <<<=== PDD:  function ~FixParticledistributionDiscrete(). \n");
//		if(logfile) fprintf(logfile ,"\n <<<=== PDD:  function ~FixParticledistributionDiscrete(). \n");
	}

/* ----------------------------------------------------------------------*/

int FixParticledistributionDiscrete::setmask()
	{
//		if(screen) fprintf(screen ,"\n ===>>> PDD:  function setmask() \n");
//		if(logfile) fprintf(logfile ,"\n ===>>> PDD:  function setmask() \n");
		
		
		
//		if(screen) fprintf(screen ,"\n <<<=== PDD:  function setmask() \n");
//		if(logfile) fprintf(logfile ,"\n <<<=== PDD:  function setmask() \n");
		
		int mask = 0;
		return mask;
	}

/* ----------------------------------------------------------------------
   prepares the fix for a series of randomize_single() commands
   typically called once per insertion step
------------------------------------------------------------------------- */

int FixParticledistributionDiscrete::random_init_single(int ntotal)
	{
//		if(screen) fprintf(screen ,"\n ===>>> PDD:  random_init_single(int ntotal) \n");
//		if(logfile) fprintf(logfile ,"\n ===>>> PDD:  random_init_single(int ntotal) \n");
		
		ninsert = ntotal;		//ninsert is the total number of particles to be inserted
		ninserted = 0;

		for(int i = 0; i < ntemplates; i++)
				parttogen[i] = static_cast<int>(static_cast<double>(ninsert) * distweight[i] + random->uniform());	//calculating number of particles for a particular template

		ninsert = 0;
		for(int i = 0; i < ntemplates; i++)
				ninsert += parttogen[i];   	//once again calculating the ninsert (total number of particles to insert)
			
		
//		if(screen) fprintf(screen ,"\n <<<=== PDD:  random_init_single(int ntotal) \n");
//		if(logfile) fprintf(logfile ,"\n <<<=== PDD:  random_init_single(int ntotal) \n");
		
		return ninsert;
		
	}

/* ----------------------------------------------------------------------
   request one template to generate one pti
------------------------------------------------------------------------- */

Region* FixParticledistributionDiscrete::randomize_single()
	{
//		if(screen) fprintf(screen ,"\n ===>>> PDD:  randomize_single() \n");
//		if(logfile) fprintf(logfile ,"\n ===>>> PDD:  randomize_single() \n");
		
		if(ntemplates == 1)			//In our case this loop will itself terminate the loop
			{
				 templates[0]->randomize_single();
				 
				 return templates[0]->region(); 
			}

		//choose a template from the discrete distribution, beginning from large to small particles
		int chosen = 0;
		int chosendist = distorder[chosen];
		int ntoinsert = parttogen[chosendist];
		while(ninserted >= ntoinsert && chosen < ntemplates-1)
			{
				chosen++;
				chosendist = distorder[chosen];
				ntoinsert += parttogen[chosendist];			//calculating total number of particles to generate
			}

		templates[chosendist]->randomize_single();			//randomizing template (executing pti)

		pti = templates[chosendist]->pti;		//transferring pti to pti

		ninserted++;	

		return templates[chosendist]->region();
		
//		if(screen) fprintf(screen ,"\n <<<=== PDD:  randomize_single() \n");
//		if(logfile) fprintf(logfile ,"\n <<<=== PDD:  randomize_single() \n");
	}

/* ----------------------------------------------------------------------
   prepares the fix for a series of randomize_list() command
   also prepares templates
       - deletes their old lists if present and allocates new lists
   typically only called once before first insertion step

   allocates for max # particles

   can be called by multiple fix insert commands, so check first if max #
   particles to be inserted is exceeded and only re-allocate in this case
------------------------------------------------------------------------- */

void FixParticledistributionDiscrete::random_init_list(int ntotal)
{
	if(screen) fprintf(screen ,"\n ===>>> PDD:  random_init_list(int ntotal) \n");
	if(logfile) fprintf(logfile ,"\n ===>>> PDD:  random_init_list(int ntotal) \n");
	
    int parttogen_max_i, n_pti_max_requested;
    int nprocs = comm->nprocs;

    ntotal += 2 * ntemplates;

    if(screen) fprintf(screen ,"\n PDD: ntotal = %d \n", ntotal);
	if(logfile) fprintf(logfile ,"\n PDD: ntotal = %d \n", ntotal);

    // number of requested pti
    n_pti_max_requested = 0;

    for(int i = 0; i < ntemplates; i++)		//Checking for number of particles to generate//If finds true then reinitilizes the init_ptilist
		{
			parttogen_max_i = static_cast<int>(static_cast<double>(ntotal) * distweight[i] + static_cast<double>(1.01)*(ntemplates+nprocs));
			n_pti_max_requested += parttogen_max_i;

			if(screen) fprintf(screen ,"\n PDD: parttogen_max_i = %d \n", parttogen_max_i);
			if(logfile) fprintf(logfile ,"\n PDD: parttogen_max_i = %d \n", parttogen_max_i);

			// re-allocated if need more ptis in this template than allocated so far
			if(parttogen_max_i > templates[i]->n_pti_max)
				{
					templates[i]->delete_ptilist();
					templates[i]->init_ptilist(parttogen_max_i);
				}
			
		}

    // re-allocate if need more total ptis in distribution than allocated so far
    if(n_pti_max_requested > n_pti_max)
		{
			n_pti_max = n_pti_max_requested;
			if(pti_list) delete []pti_list;
			pti_list = new ParticleToInsert*[n_pti_max];
		}
		
	if(screen) fprintf(screen ,"\n <<<=== PDD:  random_init_list(int ntotal) \n");
	if(logfile) fprintf(logfile ,"\n <<<=== PDD:  random_init_list(int ntotal) \n");	
}

/* ----------------------------------------------------------------------
   tell all templates to generate their pti_list, wire their pti_list to
   the list in this fix. returns number of particles to be inserted.
   typically called once per insertion step

   for exact_number = 1, truncate distribution so to exactly meet
                               requested # particles
   for exact_number = 0, use random gen to fulfil distribution
------------------------------------------------------------------------- */

int FixParticledistributionDiscrete::randomize_list(int ntotal,int insert_groupbit,int exact_number)
	{
		
		if(screen) fprintf(screen ,"\n ===>>> PDD:  randomize_list(int ntotal,int insert_groupbit,int exact_number) \n");
		if(logfile) fprintf(logfile ,"\n ===>>> PDD:  randomize_list(int ntotal,int insert_groupbit,int exact_number) \n");
		
	//	if(ntotal > n_pti_max)
	//		{
	//			error->one(FLERR,"Faulty implementation: FixParticledistributionDiscrete::randomize_list() called for more particles than defined in random_init_list()");
	//		}

		ninsert = ntotal;
		ninserted = 0;

		if(screen) fprintf(screen ,"\n PDD: ninsert = %d \n", ninsert);
		if(logfile) fprintf(logfile ,"\n PDD: ninsert = %d \n", ninsert);

		// use random generator so long-time average of insertion will represent distribution correctly
		if(exact_number == 0)
			{
				for(int i = 0; i < ntemplates; i++)
					{
							parttogen[i] = static_cast<int>(static_cast<double>(ninsert) * distweight[i] + random->uniform());
					}
			}
		// truncate distribution so # particles to insert is met exactly
			else
			{
				int ninsert_truncated = 0, j;
				double *remainder = new double[ntemplates], rsum, r;

				// distribute particles and calculate remainder
				for(int i = 0; i < ntemplates; i++)
					{
						   parttogen[i] = static_cast<int>(static_cast<double>(ninsert) * distweight[i]);
						   ninsert_truncated += parttogen[i];
						   remainder[i] = static_cast<double>(ninsert) * distweight[i] - static_cast<double>(parttogen[i]);
					}

				int ninsert_gap = ninsert - ninsert_truncated;

				// distribute remaining ninsert_gap particles
				for(int i = 0; i < ninsert_gap; i++)
						{
								r = random->uniform() * static_cast<double>(ninsert_gap);
								j = 0;
								rsum = remainder[0];

								while(j < (ntemplates-1) && rsum < r)
									{
										rsum += remainder[j];
										j++;
									}
								
								parttogen[j]++;
						}

				delete []remainder;
			}

		// count total particle number to be inserted, let templates generate a pti_list
		ninsert = 0;
		for(int i = 0; i < ntemplates; i++)
		{
			ninsert += parttogen[i];

			if(screen) fprintf(screen ,"\n PDD: parttogen[i] = %d \n", parttogen[i]);
			if(logfile) fprintf(logfile ,"\n PDD: parttogen[i] = %d \n", parttogen[i]);

			templates[i]->randomize_ptilist(parttogen[i],groupbit | insert_groupbit);
		}

		// wire lists, make sure in correct order (large to small particles)

		n_pti = 0;
		for(int i = 0; i < ntemplates; i++)
			{
				int chosendist = distorder[i];
				for (int j = 0; j < parttogen[chosendist]; j++)
					{
						pti_list[n_pti + j] = templates[chosendist]->pti_list[j];
					}
				n_pti += parttogen[chosendist];
			}

		if(n_pti != ninsert) error->all(FLERR,"Internal error in FixParticledistributionDiscrete::randomize_list");

		ninserted = ninsert;
		
		
		
		if(screen) fprintf(screen ,"\n <<<=== PDD:  randomize_list(int ntotal,int insert_groupbit,int exact_number) \n");
		if(logfile) fprintf(logfile ,"\n <<<=== PDD:  randomize_list(int ntotal,int insert_groupbit,int exact_number) \n");
		
		return ninsert;
	}

/* ----------------------------------------------------------------------
   set particle properties - only pti needs to know which properties to set
   loop to n, not n_pti, since not all particles may have been inserted
------------------------------------------------------------------------- */

int FixParticledistributionDiscrete::insert(int n)
{
	if(screen) fprintf(screen ,"\n ===>>> PDD:  insert(int n) \n");
	if(logfile) fprintf(logfile ,"\n ===>>> PDD:  insert(int n) \n");
	
	if(screen) fprintf(screen ,"\n PDD: n = %d \n", n);
    if(logfile) fprintf(logfile ,"\n PDD: n = %d \n", n);

    int ninserted_spheres_local = 0;
    
    if(screen) fprintf(screen, "to be inserted: n = %d \n", n);
	if(logfile) fprintf(logfile, "to be inserted: n = %d \n", n);
    
    for(int i = 0; i < n; i++)
    {
        ninserted_spheres_local += pti_list[i]->insert();
        if(screen) fprintf(screen ,"\n PDD: ninserted_spheres_local = %d \n", ninserted_spheres_local);
    	if(logfile) fprintf(logfile ,"\n PDD: ninserted_spheres_local = %d \n", ninserted_spheres_local);

    }
    
    
    if(screen) fprintf(screen ,"\n <<<=== PDD:  insert(int n) \n");
	if(logfile) fprintf(logfile ,"\n <<<=== PDD:  insert(int n) \n");
	
	return ninserted_spheres_local;
}

/* ----------------------------------------------------------------------
   wrap up insertion
------------------------------------------------------------------------- */

void FixParticledistributionDiscrete::finalize_insertion()
{
	if(screen) fprintf(screen ,"\n ===>>> PDD:  finalize_insertion() \n");
	if(logfile) fprintf(logfile ,"\n ===>>> PDD:  finalize_insertion() \n"); 
	
    for(int i = 0; i < ntemplates; i++)
        templates[i]->finalize_insertion();
        
   if(screen) fprintf(screen ,"\n <<<=== PDD:  finalize_insertion() \n");
	if(logfile) fprintf(logfile ,"\n <<<=== PDD:  finalize_insertion() \n");   
}

/* ----------------------------------------------------------------------*/

double FixParticledistributionDiscrete::vol_expect()
{
//	if(screen) fprintf(screen ,"\n ===>>> PDD:  vol_expect() \n");
//	if(logfile) fprintf(logfile ,"\n ===>>> PDD:  vol_expect() \n");
    
//    if(screen) fprintf(screen ,"\n <<<=== PDD:  vol_expect() \n");
//	if(logfile) fprintf(logfile ,"\n <<<=== PDD:  vol_expect() \n");
	
	return volexpect;
}

/* ----------------------------------------------------------------------*/

double FixParticledistributionDiscrete::mass_expect()
{
//	if(screen) fprintf(screen ,"\n ===>>> PDD:  mass_expect() \n");
//	if(logfile) fprintf(logfile ,"\n ===>>> PDD:  mass_expect() \n");
    
//    if(screen) fprintf(screen ,"\n <<<=== PDD:  mass_expect() \n");
//	if(logfile) fprintf(logfile ,"\n <<<=== PDD:  mass_expect() \n");
	
	return massexpect;
}

/* ----------------------------------------------------------------------*/

double FixParticledistributionDiscrete::min_rad(int type)
{	
//	if(screen) fprintf(screen ,"\n ===>>> PDD:  min_rad(int type) \n");
//	if(logfile) fprintf(logfile ,"\n ===>>> PDD:  min_rad(int type) \n");
    //get minrad
    double minrad_type = 1000.;
    for(int i = 0; i < ntemplates;i++)
      if( templates[i]->type() == type  && templates[i]->min_rad() < minrad_type)
        minrad_type = templates[i]->min_rad();

    
    
//    if(screen) fprintf(screen ,"\n <<<=== PDD:  min_rad(int type) \n");
//	if(logfile) fprintf(logfile ,"\n <<<=== PDD:  min_rad(int type) \n");
	
	return minrad_type;
}

/* ----------------------------------------------------------------------*/

double FixParticledistributionDiscrete::max_rad(int type)
{	
//	if(screen) fprintf(screen ,"\n ===>>> PDD:  max_rad(int type) \n");
//	if(logfile) fprintf(logfile ,"\n ===>>> PDD:  max_rad(int type) \n");
    //get maxrad
    double maxrad_type = 0.;
    for(int i = 0; i < ntemplates;i++)
		if( templates[i]->type() == type  && templates[i]->max_rad() > maxrad_type)
			maxrad_type = templates[i]->max_rad();

    
//    if(screen) fprintf(screen ,"\n <<<=== PDD:  max_rad(int type) \n");
//	if(logfile) fprintf(logfile ,"\n <<<=== PDD:  max_rad(int type) \n");
	
	return maxrad_type;
}

/* ----------------------------------------------------------------------*/

int FixParticledistributionDiscrete::max_type()
{
//	if(screen) fprintf(screen ,"\n ===>>> PDD:  max_type() \n");
//	if(logfile) fprintf(logfile ,"\n ===>>> PDD:  max_type() \n");
    
//    if(screen) fprintf(screen ,"\n <<<=== PDD:  max_type() \n");
//	if(logfile) fprintf(logfile ,"\n <<<=== PDD:  max_type() \n");
	return maxtype;
}

/* ----------------------------------------------------------------------*/

int FixParticledistributionDiscrete::min_type()
{
//	if(screen) fprintf(screen ,"\n ===>>> PDD:  min_type() \n");
//	if(logfile) fprintf(logfile ,"\n ===>>> PDD:  min_type() \n");
    
//    if(screen) fprintf(screen ,"\n <<<=== PDD:  min_type() \n");
//	if(logfile) fprintf(logfile ,"\n <<<=== PDD:  min_type() \n");
	
	return mintype;
}

/* ----------------------------------------------------------------------*/

int FixParticledistributionDiscrete::max_nspheres()

{
//	if(screen) fprintf(screen ,"\n ===>>> PDD:  min_type() \n");
//	if(logfile) fprintf(logfile ,"\n ===>>> PDD:  min_type() \n");
    
//    if(screen) fprintf(screen ,"\n <<<=== PDD:  min_type() \n");
//	if(logfile) fprintf(logfile ,"\n <<<=== PDD:  min_type() \n");
	
	return maxnspheres;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixParticledistributionDiscrete::write_restart(FILE *fp)
{
	if(screen) fprintf(screen ,"\n ===>>> PDD:  write_restart(FILE *fp) \n");
	if(logfile) fprintf(logfile ,"\n ===>>> PDD:  write_restart(FILE *fp) \n");
		  int n = 0;
		  double list[1];
		  list[n++] = static_cast<int>(random->state());

		  if (comm->me == 0) {
			int size = n * sizeof(double);
			fwrite(&size,sizeof(int),1,fp);
			fwrite(list,sizeof(double),n,fp);
		  }
	if(screen) fprintf(screen ,"\n <<<=== PDD:  write_restart(FILE *fp) \n");
	if(logfile) fprintf(logfile ,"\n <<<=== PDD:  write_restart(FILE *fp) \n");
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixParticledistributionDiscrete::restart(char *buf)
{
//	if(screen) fprintf(screen ,"\n ===>>> PDD:  restart(char *buf) \n");
//	if(logfile) fprintf(logfile ,"\n ===>>> PDD:  restart(char *buf) \n");
	  int n = 0;
	  double *list = (double *) buf;

	  seed = static_cast<int> (list[n++]) + comm->me;

	  random->reset(seed);
//	if(screen) fprintf(screen ,"\n <<<=== PDD:  restart(char *buf) \n");
//	if(logfile) fprintf(logfile ,"\n <<<=== PDD:  restart(char *buf) \n");  
}
