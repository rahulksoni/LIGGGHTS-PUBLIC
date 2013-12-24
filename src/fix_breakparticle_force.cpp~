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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_breakparticle_force.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_particledistribution_discrete.h"
#include "fix_template_multiplespheres.h"
#include "particleToInsert.h"

using namespace LAMMPS_NS;
using namespace FixConst; //needed for END_OF_STEP

/* ---------------------------------------------------------------------- */

FixBreakparticleForce::FixBreakparticleForce(LAMMPS *lmp, int narg, char **arg) :
  FixInsert(lmp, narg, arg)
{
  // set defaults first, then parse args
  init_defaults();

  bool hasargs = true;
  while(iarg < narg && hasargs)
    {
    hasargs = false;
    if (strcmp(arg[iarg],"volumefraction") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix breakparticle/force command");
      volumefraction = atof(arg[iarg+1]);
      if(volumefraction < 0. || volumefraction >= 1.)
        error->all(FLERR,"Illegal fix breakparticle/force command, 0 < 'volumefraction' < 1 required");
      iarg += 2;
      hasargs = true;
    }else if (strcmp(arg[iarg],"force_treshold") == 0) 
    {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix breakparticle/force command");
      f_break = atof(arg[iarg+1]);
      if(f_break <= 0. ) error->all(FLERR,"Illegal fix breakparticle/force command, 'force_treshold' must be > 0");
      iarg += 2;
      hasargs = true;
    } else error->all(FLERR,"Illegal fix breakparticle/force command, unknown keyword");
  }

  // no fixed insertion target (such as # particles) since insertion triggered by break-ups
  ninsert_exists = 0;

  // turn off overlap check and turn off start stats since both not required
  check_ol_flag = 0;
  print_stats_start_flag = 0;

  breakdata = NULL;
  fix_break = NULL;

  nevery = 1;

  n_break = 0;
  mass_break = 0.;

  if(maxrad >= 1.) error->all(FLERR,"Fix breakparticle/force: Particle distribution must be relative, max radius mus be < 1");
}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::post_create()
{
  FixInsert::post_create();

  if(!fix_break)
  {
        char *breakvar_name = new char[10+strlen(id)];
        char* fixarg[9];
        fixarg[0] = new char[10+strlen(id)];
        sprintf(fixarg[0],"break_%s",id);
        
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3] = new char[10+strlen(id)];
        sprintf(breakvar_name,"break_%s",id);
        strcpy(fixarg[3],breakvar_name);
        fixarg[4]="scalar"; 
        fixarg[5]="yes";    
        fixarg[6]="yes";    
        fixarg[7]="no";    
        fixarg[8]="0.";
        fix_break = modify->add_fix_property_atom(9,fixarg,style);
        delete []fixarg[0];
        delete []fixarg[3];
        delete []breakvar_name;
  }

}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::pre_delete()
{
    modify->delete_fix(fix_break->id);
}

/* ---------------------------------------------------------------------- */

FixBreakparticleForce::~FixBreakparticleForce()
{
    if(breakdata) memory->sfree(breakdata);
}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::init_defaults()
{
    f_break = 0.;
    n_fragments = 0;

    volumefraction = 0.5;

    fix_fragments = NULL;
}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::init()
{
    FixInsert::init();
}

/* ----------------------------------------------------------------------
   calculate ninsert, insert_every, ninsert_per, massinsert, flowrates etc
   also perform error checks
------------------------------------------------------------------------- */

void FixBreakparticleForce::calc_insertion_properties()
{
    // error checks
    if(f_break == 0.)
        error->all(FLERR,"Illegal fix breakparticle/force command, you have to specify 'force_treshold'");
    if(nflowrate > 0. || massflowrate > 0.)
        error->all(FLERR,"Illegal fix breakparticle/force command, specifying 'nflowrate' or 'massflowrate' is not allowed");
    if(ninsert > 0 || massinsert > 0.)
        error->all(FLERR,"Illegal fix breakparticle/force command, specifying 'nparticles' or 'mass' is not allowed");
    if(insert_every <= 0)
        error->all(FLERR,"Illegal fix breakparticle/force command, specifying 'every' must be > 0");

    // fix holding particle fragments
    if(fix_distribution->n_particletemplates() != 1)
        error->all(FLERR,"Illegal fix breakparticle/force command, fix of type particledistribution/discrete must hold exactly one template");
    if(strcmp(fix_distribution->particletemplates()[0]->style,"particletemplate/multiplespheres"))
        error->all(FLERR,"Illegal fix breakparticle/force command, fix of type particledistribution/discrete must hold exactly one template of type fix particletemplate/multiplespheres");
    fix_fragments = static_cast<FixTemplateMultiplespheres*>(fix_distribution->particletemplates()[0]);

    // get number of fragments from fix
    n_fragments = fix_fragments->number_spheres();

    // do not need ninsert_per
}

/* ---------------------------------------------------------------------- */
//this is needed to set the end-of-step function into the pipeline,
//END_OF_STEP is defined in namespace FixConst, so I added it above !
int FixBreakparticleForce::setmask()
{
    int mask = FixInsert::setmask();
    mask |= END_OF_STEP;
    return mask;
}

/* ---------------------------------------------------------------------- */

inline int FixBreakparticleForce::is_nearby(int i)
{
    // need not check overlap with existing particles since we
    // know space originally taken by deleted particles is free
    return 0;
}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::end_of_step()
{
    int nlocal = atom->nlocal;
    double *flag = fix_break->vector_atom;
    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;

    double f_sqr,f_break_sqr;

    f_break_sqr = f_break * f_break;

    // check breakage criterion for all local particles

    for(int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            f_sqr = vectorMag3DSquared(f[i]);
            if(f_sqr >= f_break_sqr)
                flag[i] = 1.;
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::pre_insert()
{
    int i,ibreak;
    int nlocal = atom->nlocal;
    int *mask = atom->mask;
    double **x = atom->x;
    double **v = atom->v;
    double *radius = atom->radius;
    double *rmass = atom->rmass;
    double *flag = fix_break->vector_atom;
    AtomVec *avec = atom->avec;

    // count # of particles to break

    n_break_this_local = 0;
    mass_break_this_local = 0.;
    for( i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit && flag[i] == 1.)
        {
            n_break_this_local++;
            mass_break_this_local += rmass[i];
        }
    }

   // tally stats

   LAMMPS_NS::MPI_Sum_Scalar(n_break_this_local,n_break_this,world);
   n_break += n_break_this;
   LAMMPS_NS::MPI_Sum_Scalar(mass_break_this_local,mass_break_this,world);
   mass_break += mass_break_this;
  

   // allocate breakage data

   //if(breakdata) memory->destroy_2d_double_array(breakdata);										
   if(breakdata) memory->destroy(breakdata);					
			//Please note that the appropriate function for this isn't available in memory.h
   //breakdata = memory->create_2d_double_array(n_break_this_local,7,"FixBreakparticleForce::breakdata");
   breakdata = memory->create(FixBreakparticleForce::breakdata,n_break_this_local,7,"FixBreakparticleForce::breakdata");		
         	//Please note that the appropriate function for this isn't available in memory.h in the form of create(int,int,double **) where FixBreakparticleForce::breakdata results 2d array as it can be seen at fix_breakparticle_force.h

   // fill breakage data and remove particles
   

   i = ibreak = 0;
   while (i < nlocal)
   {
      if (mask[i] & groupbit && flag[i] == 1.)
      {
          //copy data needed for insertion
          vectorCopy3D(x[i],&breakdata[ibreak][0]);
          vectorCopy3D(v[i],&breakdata[ibreak][3]);
          breakdata[ibreak][6] = radius[i];
          ibreak++;

          // delete particle
          avec->copy(nlocal-1,i,0);
          nlocal--;
      }
      else i++;
   }

   // update local and global # particles
   
   atom->nlocal = nlocal;
   double rlocal = static_cast<double>(atom->nlocal);
   MPI_Allreduce(&rlocal,&atom->natoms,1,MPI_DOUBLE,MPI_SUM,world);

   // print stats
   print_stats_breakage_during();
}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::print_stats_breakage_during()
{
  int step = update->ntimestep;

  if (me == 0 && n_break_this > 0)
  {
    if (screen)
      fprintf(screen ,"Particle breakage: broke %d particles (mass %f) at step %d\n - a total of %d particles (mass %f) broken so far.\n",
	      n_break_this,mass_break_this,step,n_break,mass_break);

    if (logfile)
      fprintf(logfile,"Particle breakage: broke %d particles (mass %f) at step %d\n - a total of %d particles (mass %f) broken so far.\n",
	      n_break_this,mass_break_this,step,n_break,mass_break);
  }
}

/* ---------------------------------------------------------------------- */

int FixBreakparticleForce::calc_ninsert_this()
{
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
    double pos_ins[3],v_ins[3],omega_ins[3],quat_ins[4],rad_broken;
    int iparticle, nins;
    ParticleToInsert *pti;

    vectorZeroize3D(omega_ins);
    vectorZeroize4D(quat_ins);

    // local insertion
    double mass_inserted_this_local = 0.;
    int ninserted_this_local = 0;
    int ninserted_spheres_this_local = 0;

    // global insertion
    ninserted_this = ninserted_spheres_this = 0;
    mass_inserted_this = 0.;

    // n_break_this ptis with n_fragments spheres each
    // n_break_this * n_fragments spheres to be generated

    iparticle = 0;

    while(iparticle < n_break_this_local)
    {
        
        // get position, velocity and radius of broken particle
        vectorCopy3D(&breakdata[iparticle][0],pos_ins);
        vectorCopy3D(&breakdata[iparticle][3],v_ins);
        rad_broken = breakdata[iparticle][6];

        //fprintf(screen,"rad_broken %f\n",rad_broken);

        // get pti and scale it down with radius of broken particle
        pti = fix_distribution->pti_list[iparticle];
        pti->scale_pti(rad_broken);

        nins = pti->set_x_v_omega(pos_ins,v_insert,omega_ins,quat_insert);

        // tally stats
        ninserted_spheres_this_local += nins;
        mass_inserted_this_local += pti->mass_ins;
        ninserted_this_local++;

        iparticle++;
    }

    // tally stats, have to do this since operation is locally on each process
    // as opposed to e.g. FixInsertPack::x_v_omega()

  LAMMPS_NS::MPI_Sum_Scalar(ninserted_spheres_this_local,ninserted_spheres_this,world);
  LAMMPS_NS::MPI_Sum_Scalar(ninserted_this_local,ninserted_this,world);
  LAMMPS_NS::MPI_Sum_Scalar(mass_inserted_this_local,mass_inserted_this,world);
    
   
}

//not sure if we realy need this function, however it needs to be a function with this name here for the parent class fix-instert
double FixBreakparticleForce::insertion_fraction()
{
    
    return volumefraction;
}
