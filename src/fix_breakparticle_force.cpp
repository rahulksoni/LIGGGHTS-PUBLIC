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
#include "fix_breakparticle_force.h"                   /***********/
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
<<<<<<< HEAD
<<<<<<< HEAD
#include "fix_particledistribution_discrete.h"			/**************/
#include "fix_template_multiplespheres.h"				/************/
#include "particleToInsert.h"							/************/
=======
#include "fix_particledistribution_discrete.h"
#include "fix_template_multiplespheres.h"
#include "particleToInsert.h"
>>>>>>> f672010c2ceb717c573186ac1cba9f340fa00056
=======
#include "fix_particledistribution_discrete.h"			/**************/
#include "fix_template_multiplespheres.h"				/************/
#include "particleToInsert.h"							/************/
>>>>>>> 177ff7b221ad559063dd8bb1249256f1abdf06f0

using namespace LAMMPS_NS;
using namespace FixConst; //needed for END_OF_STEP

/* ---------------------------------------------------------------------- */

FixBreakparticleForce::FixBreakparticleForce(LAMMPS *lmp, int narg, char **arg) :
  FixInsert(lmp, narg, arg)   //FixBreakparticleForce construction, inclusion of FixInsert, narg is no of arguments in fixbreakparticle/force command, **arg is corresponding keywords//
{
  // set defaults first, then parse args
  init_defaults();      	//comes from fix_breakparticle_force.h through fix_insert.cpp through fix_insert.h//
  
									//Meaning of init_defaulst()//
/*
									 void FixInsert::init_defaults()
{
  // default is that total # of particles to insert by this command is known
  ninsert_exists = 1;
  ninsert = ninserted = 0;  massinsert = massinserted = 0.;  nflowrate = massflowrate = 0.;  insert_every = -1;  ninsert_per = 0.;
  // 1st insertion on next timestep is default
  first_ins_step = update->ntimestep + 1;
  maxattempt = 50;  check_ol_flag = 1;  all_in_flag = 0;  exact_number = 1;  v_randomSetting = 0; vectorZeroize3D(v_insert);
  vectorZeroize3D(v_insertFluct);  vectorZeroize3D(omega_insert);  quatUnitize4D(quat_insert); quat_random_ = false;
  print_stats_during_flag = 1;
}   */

  bool hasargs = true;			//narg is the no of total arguments and iarg is an element in that array//
  while(iarg < narg && hasargs)
    {
		hasargs = false;		//actually there is no meaning of this, bacause it is not in use//
		if (strcmp(arg[iarg],"volumefraction") == 0)    //if comparison fails//
		{
			if (iarg+2 > narg) error->all(FLERR,"Illegal fix breakparticle/force command");   //if iarg is first element then narg cannot be more than 3//
			volumefraction = atof(arg[iarg+1]);     
			if(volumefraction < 0. || volumefraction >= 1.)		//checking if volumefraction has a legitimate value//
				error->all(FLERR,"Illegal fix breakparticle/force command, 0 < 'volumefraction' < 1 required");
			iarg += 2;		//increasing iarg value by +2//
			hasargs = true;
		}else if (strcmp(arg[iarg],"force_treshold") == 0)  //if comparison fails//
			{
			if (iarg+2 > narg) error->all(FLERR,"Illegal fix breakparticle/force command");
			f_break = atof(arg[iarg+1]);
			if(f_break <= 0. ) error->all(FLERR,"Illegal fix breakparticle/force command, 'force_treshold' must be > 0");
			iarg += 2;
			hasargs = true;
		} else error->all(FLERR,"Illegal fix breakparticle/force command, unknown keyword");
    }

  // no fixed insertion target (such as # particles) since insertion triggered by break-ups
  ninsert_exists = 0;   //changed from default value 1 to 0//

  // turn off overlap check and turn off start stats since both not required
  check_ol_flag = 0;	//changed from default value 1 to 0//
  print_stats_start_flag = 0;

  breakdata = NULL; 	//No breakdata yet, initiation of variable//
  fix_break = NULL;		//No atom properties yet for new particles, default initiation//

  nevery = 1;	//action after every 1 timestep//

  n_break = 0;	//initiation//
  mass_break = 0.;	//initiation//

  if(maxrad >= 1.) error->all(FLERR,"Fix breakparticle/force: Particle distribution must be relative, max radius mus be < 1");
}	//this means new radius is given in relative terms//

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::post_create()   //constructing function that was mentioned in fix_breakparticle_force.h//
{
  FixInsert::post_create();

  if(!fix_break)		//if fix_break is not defined then it is being defined//
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
					//it will be fixarg[break_id all property/atom break_id scalar yes yes no 0]//
					
        fix_break = modify->add_fix_property_atom(9,fixarg,style);		//class FixPropertyAtom* add_fix_property_atom(int narg,char **arg,const char *)   defined in modify.h.//
        
        delete []fixarg[0];  //delete break_id//
        delete []fixarg[3];	 //delete break_id//
        delete []breakvar_name;		//delete breakvar_name break_id//
  }

}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::pre_delete() //pre_delete declared in fix_breakparticle_force.h//
{
    modify->delete_fix(fix_break->id);	//deleteing id of fix_break//
}

/* ---------------------------------------------------------------------- */

FixBreakparticleForce::~FixBreakparticleForce()		//destructor//
{
    if(breakdata) memory->sfree(breakdata);     //deleting any available breakadata//
}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::init_defaults()   //default constructor//
{
    f_break = 0.;	//default for thresold force//
    n_fragments = 0; 
    volumefraction = 0.5;	
    fix_fragments = NULL;  //fix_fragments is a class pointer of FixTemplateMultispheres
}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::init()
{
    FixInsert::init();              //init() is defined in fix_insert.cpp//
									//assigns appropriate data to keyword multisphere//
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
    if(nflowrate > 0. || massflowrate > 0.)   //cannot assign number flowrate or massflowrate//default values are 0//
        error->all(FLERR,"Illegal fix breakparticle/force command, specifying 'nflowrate' or 'massflowrate' is not allowed");
    if(ninsert > 0 || massinsert > 0.)		//connot assign number of particles or mass to be inserted//defaults are 0//
        error->all(FLERR,"Illegal fix breakparticle/force command, specifying 'nparticles' or 'mass' is not allowed");
    if(insert_every <= 0)   //it cannot be negative//
        error->all(FLERR,"Illegal fix breakparticle/force command, specifying 'every' must be > 0");

    // fix holding particle fragments
    if(fix_distribution->n_particletemplates() != 1)    //Particles with only one type of properties can be considered//
        error->all(FLERR,"Illegal fix breakparticle/force command, fix of type particledistribution/discrete must hold exactly one template");
    if(strcmp(fix_distribution->particletemplates()[0]->style,"particletemplate/multiplespheres")) //only multisphere style is allowed//
        error->all(FLERR,"Illegal fix breakparticle/force command, fix of type particledistribution/discrete must hold exactly one template of type fix particletemplate/multiplespheres");
    
    fix_fragments = static_cast<FixTemplateMultiplespheres*>(fix_distribution->particletemplates()[0]); 	//forcing to cast accroding to multisphere template//assigning particletemplates to fix_fragments whcih is class pointer of FixTemplateMultisphere//

    // get number of fragments from fix
    n_fragments = fix_fragments->number_spheres();    //assigning number of spheres from class FixTempleteMultisphere as number of fragments//

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

void FixBreakparticleForce::end_of_step()  	//assigning properties//
{
    int nlocal = atom->nlocal;	//nlocal is number of atoms in this processor//
    double *flag = fix_break->vector_atom;		//assigning vector_atom from class FixPropertyAtom//
    double **x = atom->x;		//these are the coordinates of atoms//
    double **f = atom->f;		//force values//
    int *mask = atom->mask;

    double f_sqr,f_break_sqr;

    f_break_sqr = f_break * f_break;      //converting allotted thresold force as force acting on the particle//

    // check breakage criterion for all local particles

    for(int i = 0; i < nlocal; i++)			//checks the particles which are going to be breakon//
    {
        if (mask[i] & groupbit) 	//groupbit is the check whether particle belongs to current group or not//mask is the group number to which particle belongs//
        {
            f_sqr = vectorMag3DSquared(f[i]);		//current force acting on particle//
            if(f_sqr >= f_break_sqr)				//current force acting on particle is greater than thresold force//
                flag[i] = 1.;						//flag i=1 that will make that particle eligible for breaking//
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixBreakparticleForce::pre_insert()	//things to be done before particle insertion//
{
    int i,ibreak;
    int nlocal = atom->nlocal;		//nlocal is number of particles in current processor//
    int *mask = atom->mask;		
    double **x = atom->x;		//particles coordinates assignment//
    double **v = atom->v;		//particles forces assignments//
    double *radius = atom->radius;		//particles radius assignments//
    double *rmass = atom->rmass;		//particles mass assignment//	
    double *flag = fix_break->vector_atom;		//assigning vector_atom from class FixPropertyAtom//
    AtomVec *avec = atom->avec;		//avec is a class pointer to AtomVec//also, defined in atom.h//

    // count # of particles to break

    n_break_this_local = 0;
    mass_break_this_local = 0.;
    for( i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit && flag[i] == 1.) 	//we have made flag 1 earlier for particles that need to be breakon////mask is group number//
        {											//groupbit is the check whether particle belongs to current group or not//
            n_break_this_local++;					//total number of particles to be broken//
            mass_break_this_local += rmass[i];		//total mass of particle to be broken//rmass[i] is mass of particle based on r//
        }
    }

   // tally stats

   LAMMPS_NS::MPI_Sum_Scalar(n_break_this_local,n_break_this,world); 	//defined in mpi.h//this function actually performs copying the memeory into the destination///copies memory equivalent to int in n_break_this_local to _break_this//
   n_break += n_break_this;    //Getting total number of particles to be breakon//
   LAMMPS_NS::MPI_Sum_Scalar(mass_break_this_local,mass_break_this,world);     //copies memory equivalent to doubles in mass_break_this_local to mass_break_this//
   mass_break += mass_break_this;
  

   // allocate breakage data

   //if(breakdata) memory->destroy_2d_double_array(breakdata);										
   if(breakdata) memory->destroy(breakdata);					
			//Please note that the appropriate function for this isn't available in memory.h
   //breakdata = memory->create_2d_double_array(n_break_this_local,7,"FixBreakparticleForce::breakdata");
   breakdata = memory->create(FixBreakparticleForce::breakdata,n_break_this_local,7,"FixBreakparticleForce::breakdata");  //defined in memory.h//Its all about safe memory allocation//Its a 2d array (where each element can multiple values) with no of rows=n_break_this_local, no of columns=7, for 7 properties that are defined later in this file//
         	//Please note that the appropriate function for this wasn't available in memory.h in the form of create(int,int,double **)//FixBreakparticleForce::breakdata was added which results 2d array as it can be seen at fix_breakparticle_force.h

   // fill breakage data and remove particles
   

   i = ibreak = 0;
   while (i < nlocal) 	//nlocal is numbewr of particles in current processor//
   {
      if (mask[i] & groupbit && flag[i] == 1.)  //if turns true then data copy, ibreak++, nlocal-- otherwise i++//copying data of original particles//
      {
          //copy data needed for insertion
          vectorCopy3D(x[i],&breakdata[ibreak][0]);	//function defined in vector_liggghts.h//copying from left to right//in a 2d array ibreak is row no and 0 is column no//copying coordinates in first column//
          vectorCopy3D(v[i],&breakdata[ibreak][3]);	//copying from left to right//copying vecoties in 4th column//
          breakdata[ibreak][6] = radius[i];		//putting radius value in 7th column//
          ibreak++;

          // delete existing particle that is going to be broken//
          avec->copy(nlocal-1,i,0);   //command through atomvec.h//avec is a class pointer to AtomVec//after copying coordinates, velocities and radius of parent particle it is going to be deleted//
          nlocal--;					  //after deleteing one particle new no of nlocal is obviously nlocal--//
          //i will remain 0 until if condition turns true//if all cases if is true then nlocal will automatically tends to 0//
          //if in any case if turns false then no data copy rather i++. so, i is the count of number of particles noy eligible for break//
          //So, i increases if particle found with ineligibility of breaking and nlocal reduces if particle with eligibility to be broken. Ultimately, both are tending to become equal at some point of time//
          //That means at that point i=nlocal and therefore false for while loop//
      }
      else i++;
   }

   // update local and global # particles
   
   atom->nlocal = nlocal;		//updating no of particles in current processor//Please note nlocal-- in previous while loop will not change the original value of nlocal, bacause nlocal-- is in while scope//
   double rlocal = static_cast<double>(atom->nlocal);
   MPI_Allreduce(&rlocal,&atom->natoms,1,MPI_DOUBLE,MPI_SUM,world);    	//copies appropriate memory from rlocal to natoms//

   // print stats
   print_stats_breakage_during();     //print stats in terminal//
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

    vectorZeroize3D(omega_ins);		//turning all values to zero//
    vectorZeroize4D(quat_ins);		//turning all values to zero//

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
        vectorCopy3D(&breakdata[iparticle][0],pos_ins);			//copies positions from left to right//
        vectorCopy3D(&breakdata[iparticle][3],v_ins);			//copies velocities from left to right//
        rad_broken = breakdata[iparticle][6];			//copies original particle radius to rad_broken//

        //fprintf(screen,"rad_broken %f\n",rad_broken);

	        // get pti and scale it down with radius of broken particle
        pti = fix_distribution->pti_list[iparticle];     //fixdistribution is a class pointer to FixParticledistributionDiscrete (declared in fix_insert.h, class defined in fix_particledistribution_distribution_discrete.h)//pti_list is double array of class ParticleToInsert (defined in fix_particledistribution_distribution_discrete.h, class defined in particletoinsert.h)//
        pti->scale_pti(rad_broken);  	//scale_pti is declared in particletoinsert.h//function defined in particle to insert.cpp//
										//scale_pti(r_scale) multiplies poistions,velocities and r_bound_ins with r_scale, and volume and mass with r_scale^3//
										//These are relative operations and actual values will be set up in next line//
        nins = pti->set_x_v_omega(pos_ins,v_insert,omega_ins,quat_insert); 	//function defined in particletoinsert.cpp//
						//Copies velocities and omega. Returns values nspheres to nins. Add relative values to position//

        // tally stats
        ninserted_spheres_this_local += nins;   	//Please note it was 0 earlier so now  ninserted_spheres_this_local=nins//
        mass_inserted_this_local += pti->mass_ins;	//Again mass_inserted_this_local=mass_ins
        ninserted_this_local++;		//increment for number of particles inserted////

        iparticle++;			//while loop increment//
    }

    // tally stats, have to do this since operation is locally on each process
    // as opposed to e.g. FixInsertPack::x_v_omega()

  LAMMPS_NS::MPI_Sum_Scalar(ninserted_spheres_this_local,ninserted_spheres_this,world);	//copies memory//
  LAMMPS_NS::MPI_Sum_Scalar(ninserted_this_local,ninserted_this,world);
  LAMMPS_NS::MPI_Sum_Scalar(mass_inserted_this_local,mass_inserted_this,world);
    
   
}

//not sure if we realy need this function, however it needs to be a function with this name here for the parent class fix-instert
double FixBreakparticleForce::insertion_fraction()
{
    
    return volumefraction;
}
