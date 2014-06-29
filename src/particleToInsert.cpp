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

#include "particleToInsert.h"
#include "math.h"
#include "error.h"
#include "update.h"
#include "domain.h"
#include "atom.h"
#include "atom_vec.h"
#include "fix.h"
#include "vector_liggghts.h"
#include "modify.h"

using namespace LAMMPS_NS;

ParticleToInsert::ParticleToInsert(LAMMPS* lmp,int ns) : Pointers(lmp)
{
	if(screen) fprintf(screen ,"\n ===>>> PTI:  ParticleToInsert(LAMMPS* lmp,int ns) \n");
	if(logfile) fprintf(logfile ,"\n ===>>> PTI:  ParticleToInsert(LAMMPS* lmp,int ns) \n");
        groupbit = 0;

        nspheres = ns;

        memory->create(x_ins,nspheres,3,"x_ins");
        radius_ins = new double[nspheres];

    if(screen) fprintf(screen ,"\n PTI: nspheres = %d \n", nspheres);
    if(logfile) fprintf(logfile ,"\n PTI: nspheres = %d \n", nspheres);

    if(screen) fprintf(screen ,"\n <<<=== PTI:  ParticleToInsert(LAMMPS* lmp,int ns) \n");
	if(logfile) fprintf(logfile ,"\n <<<=== PTI:  ParticleToInsert(LAMMPS* lmp,int ns) \n");
}

/* ---------------------------------------------------------------------- */

ParticleToInsert::~ParticleToInsert()
{
	if(screen) fprintf(screen ,"\n ===>>> PTI:  ~ParticleToInsert() \n");
	if(logfile) fprintf(logfile ,"\n ===>>> PTI:  ~ParticleToInsert() \n");
        memory->destroy(x_ins);
        delete []radius_ins;
   if(screen) fprintf(screen ,"\n <<<=== PTI:  ~ParticleToInsert() \n");
	if(logfile) fprintf(logfile ,"\n <<<=== PTI:  ~ParticleToInsert() \n");
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::insert()
{
	if(screen) fprintf(screen ,"\n ===>>> PTI:  insert() \n");
	if(logfile) fprintf(logfile ,"\n ===>>> PTI:  insert() \n");
    // perform the actual insertion
    // add particles, set coordinate and radius
    // set group mask to "all" plus fix groups

    int inserted = 0;
    int nfix = modify->nfix;
    Fix **fix = modify->fix;

    for(int i = 0; i < nspheres; i++)
    {
        
        //if (domain->is_in_extended_subdomain(x_ins[i]))
        //{
                
                inserted++;
                
                atom->avec->create_atom(atom_type,x_ins[i]);
                
                int m = atom->nlocal - 1;		//nlocal is number of particles in current processor after adding up 1 particle by above line. 
				//Since arrays start from 0 we have made m = nlocal - 1
				
                atom->mask[m] = 1 | groupbit;
                vectorCopy3D(v_ins,atom->v[m]);
                vectorCopy3D(omega_ins,atom->omega[m]);
      
				atom->radius[m] = radius_ins[0];
				atom->rmass[m] = mass_ins;

                              
                atom->density[m] = density_ins;
                

                for (int j = 0; j < nfix; j++)
                   if (fix[j]->create_attribute) fix[j]->set_arrays(m);
        //}
    }
    
    
    if(screen) fprintf(screen ,"\n <<<=== PTI:  insert() \n");
	if(logfile) fprintf(logfile ,"\n <<<=== PTI:  insert() \n");
	
	/***************************************/
    return inserted;
    /***************************************/
}



/* ---------------------------------------------------------------------- */

int ParticleToInsert::check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear)
{
//	if(screen) fprintf(screen ,"\n ===>>> PTI:  check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear) \n");
//	if(logfile) fprintf(logfile ,"\n ===>>> PTI:  check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear) \n");
	
    // check sphere against all others in xnear
    // if no overlap add to xnear
    double del[3], rsq, radsum;

    vectorCopy3D(x,x_ins[0]);

    for(int i = 0; i < nnear; i++)
    {
        vectorSubtract3D(x_ins[0],xnear[i],del);
        rsq = vectorMag3DSquared(del);
        
        radsum = radius_ins[0] + xnear[i][3];

        // no success in overlap
        if (rsq <= radsum*radsum) return 0;
    }

    // no overlap with any other - success

    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // add to xnear
    vectorCopy3D(x_ins[0],xnear[nnear]);
    xnear[nnear][3] = radius_ins[0];
    nnear++;

    return 1;
    
//    if(screen) fprintf(screen ,"\n <<<=== PTI:  check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear) \n");
//	if(logfile) fprintf(logfile ,"\n <<<=== PTI:  check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear) \n");
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::set_x_v_omega(double *x, double *v, double *omega, double *quat)
{
	if(screen) fprintf(screen ,"\n ===>>> PTI:  set_x_v_omega(double *x, double *v, double *omega, double *quat) \n");
	if(logfile) fprintf(logfile ,"\n ===>>> PTI:  set_x_v_omega(double *x, double *v, double *omega, double *quat) \n");
    // add insertion position
    // relative position of spheres to each other already stored at this point
    if(screen) fprintf(screen ,"\n PTI: nspheres = %d \n", nspheres);
    if(logfile) fprintf(logfile ,"\n PTI: nspheres = %d \n", nspheres);

    for(int j = 0; j < nspheres; j++)
        vectorAdd3D(x_ins[j],x,x_ins[j]);		//adds relative positions to exact positions//

    // set velocity and omega
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    
    
    if(screen) fprintf(screen ,"\n <<<=== PTI:  set_x_v_omega(double *x, double *v, double *omega, double *quat) \n");
	if(logfile) fprintf(logfile ,"\n <<<=== PTI:  set_x_v_omega(double *x, double *v, double *omega, double *quat) \n");
    
    return nspheres;
}

void ParticleToInsert::set_r_mass_vol_rboundins(double r)
{
//	if(screen) fprintf(screen ,"\n ===>>> PTI:  set_r_mass_vol_rboundins(double *r_sphere) \n");
//	if(logfile) fprintf(logfile ,"\n ===>>> PTI:  set_r_mass_vol_rboundins(double *r_sphere) \n");
    
	radius_ins[0] = r;
	r_bound_ins = radius_ins[0];
    volume_ins = 4.0 * 3.14159265359 * radius_ins[0] * radius_ins[0] * radius_ins[0] / 3.0;
    mass_ins = volume_ins * density_ins;

    if(logfile) fprintf(logfile ,"\n PTI:  set_r_mass_vol_rboundins: r = %f, vol = %f, mass = %f \n", radius_ins[0], volume_ins, mass_ins);
    if(screen) fprintf(screen ,"\n PTI:  set_r_mass_vol_rboundins: r = %f, vol = %f, mass = %f \n", radius_ins[0], volume_ins, mass_ins);
    
//    if(screen) fprintf(screen ,"\n <<<=== PTI:  set_r_mass_vol_rboundins(double *r_sphere) \n");
//	if(logfile) fprintf(logfile ,"\n <<<=== PTI:  set_r_mass_vol_rboundins(double *r_sphere) \n");
    
}

/* ---------------------------------------------------------------------- */

void ParticleToInsert::scale_pti(double r_scale)
{
	if(screen) fprintf(screen ,"\n ===>>> PTI:  scale_pti(double r_scale) \n");
	if(logfile) fprintf(logfile ,"\n ===>>> PTI:  scale_pti(double r_scale) \n");
	
    double r_scale3 = r_scale*r_scale*r_scale;

    for(int i = 0; i < nspheres; i++)
		{
			radius_ins[i] *= r_scale;
			vectorScalarMult3D(x_ins[i],r_scale);
		}

    volume_ins *= r_scale3;
    mass_ins *= r_scale3;

    r_bound_ins *= r_scale;
    
    if(screen) fprintf(screen ,"\n <<<=== PTI:  scale_pti(double r_scale) \n");
	if(logfile) fprintf(logfile ,"\n <<<=== PTI:  scale_pti(double r_scale) \n");
}
