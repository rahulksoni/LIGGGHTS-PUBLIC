/*********************Rahul K Soni************************/

/****************Variables for breakage*******************/

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <cstdlib>		//for accessing rand number generator function
#include <ctime>		//for accessing current time to create a seed value for the random number generator

/********************Functions to create random 3D coordinates of newly generated particles**********************/
/*											   * *			  z | 
											 *---	*			|
									        * |__|	 *			|______	
									        *		 *				  y	
									         *      *
											   *  *
*//**********************Change min/max values as per suitability or discretionary analysis**********************/	
using namespace std;

double x_fill_min = 0.0; 
double x_fill_max = 0.0;
double y_fill_min = 0.0;
double y_fill_max = 0.0;
double z_fill_min = 0.0;
double z_fill_max = 0.0;

double x_fill_min_bin = 0.0; 
double x_fill_max_bin = 0.0;
double y_fill_min_bin = 0.0;
double y_fill_max_bin = 0.0;
double z_fill_min_bin = 0.0;
double z_fill_max_bin = 0.0;

void fill_bounds(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, double bin_x_min, double bin_x_max, double bin_y_min, double bin_y_max, double bin_z_min, double bin_z_max)
	{				
	//	cout<<"mill_axis = "<<mill_axis<<"  x_min ="<<x_min<<"  x_max ="<<x_max<<"  y_min ="<<y_min<<"  y_max ="<<y_max<<"  z_min ="<<z_min<<"  z_max ="<<z_max<<endl<<endl;
	//	double along_axis = 0.9; 
	//	double perpedicular_axis = 0.6;

	double factor = 0.95;	
	/*	if(mill_axis == 1)
			{/*
				x_fill_min = along_axis * x_min; 
				x_fill_max = along_axis * x_max; 
				y_fill_min = 0.0 * y_min; 
				y_fill_max = perpedicular_axis * y_max; 
				z_fill_min = 0.0 * z_min; 
				z_fill_max = perpedicular_axis * z_max; */
				x_fill_min = factor * x_min; 
				x_fill_max = factor * x_max; 
				y_fill_min = factor * y_min; 
				y_fill_max = factor * y_max; 
				z_fill_min = factor * z_min; 
				z_fill_max = factor * z_max;

				x_fill_min_bin = factor * bin_x_min; 
				x_fill_max_bin = factor * bin_x_max; 
				y_fill_min_bin = factor * bin_y_min; 
				y_fill_max_bin = factor * bin_y_max; 
				z_fill_min_bin = factor * bin_z_min; 
				z_fill_max_bin = factor * bin_z_max;
	/*		}else if(mill_axis == 2)
			{
				x_fill_min = perpedicular_axis * x_min; 
				x_fill_max = 0.0 * x_max; 
				y_fill_min = along_axis * y_min; 
				y_fill_max = along_axis * y_max; 
				z_fill_min = 0.0 * z_min; 
				z_fill_max = perpedicular_axis * z_max; 
			}else if(mill_axis == 3)
			{
				x_fill_min = perpedicular_axis * x_min; 
				x_fill_max = 0.0 * x_max; 
				y_fill_min = perpedicular_axis * y_min; 
				y_fill_max = 0.0 * y_max; 
				z_fill_min = along_axis * z_min; 
				z_fill_max = along_axis * z_max; 
			}else if(mill_axis == -1)
			{
				x_fill_min = along_axis * x_min; 
				x_fill_max = along_axis * x_max; 
				y_fill_min = perpedicular_axis * y_min; 
				y_fill_max = 0.0 * y_max; 
				z_fill_min = 0.0 * z_min; 
				z_fill_max = perpedicular_axis * z_max; 
			}else if(mill_axis == -2)
			{
				x_fill_min = 0.0 * x_min; 
				x_fill_max = perpedicular_axis * x_max; 
				y_fill_min = along_axis * y_min; 
				y_fill_max = along_axis * y_max; 
				z_fill_min = 0.0 * z_min; 
				z_fill_max = perpedicular_axis * z_max; 
			}else if(mill_axis == -3)
			{
				x_fill_min = 0.0 * x_min; 
				x_fill_max = perpedicular_axis * x_max; 
				y_fill_min = perpedicular_axis * y_min; 
				y_fill_max = 0.0 * y_max; 
				z_fill_min = along_axis * z_min; 
				z_fill_max = along_axis * z_max; 
			}
			
	//			if(screen) fprintf(screen ,"\n  x_fill_min = %f, x_fill_max = %f, y_fill_min = %f, y_fill_max = %f, z_fill_min = %f, z_fill_max = %f \n",x_fill_min,x_fill_max,y_fill_min,y_fill_max,z_fill_min,z_fill_max);
	//			if(logfile) fprintf(logfile ,"\n  x_fill_min = %f, x_fill_max = %f, y_fill_min = %f, y_fill_max = %f, z_fill_min = %f, z_fill_max = %f \n",x_fill_min,x_fill_max,y_fill_min,y_fill_max,z_fill_min,z_fill_max);
			
	//		cout<<"  x_fill_min ="<<x_fill_min<<"  x_fill_max ="<<x_fill_max<<"  y_fill_min ="<<y_fill_min<<"  y_fill_max ="<<y_fill_max<<"  z_fill_min ="<<z_fill_min<<"  z_fill_max ="<<z_fill_max<<endl<<endl;
	*/		
	}		
	
 
//int random_number;			
			
double random_coordinates_generator(int dimension, int random_number)	//dimension defines whether it is to be done for x or y or z. 0 for x, 1 for y and 2 for z
	{
		random_number *= 233;
		srand(random_number); 			//Current time value as seed for random number generation//	
			
		if(dimension == 0)
			{
			//	srand(count);
				double rrr = ((double)rand()) / ((double)(RAND_MAX));
				rrr = (rrr * (x_fill_max - x_fill_min)) + x_fill_min;
		
				return rrr;
				
			}else if(dimension == 1)
			{
			//	srand(count);
				double rrr = ((double)rand()) / ((double)(RAND_MAX));
				rrr = (rrr * (y_fill_max - y_fill_min)) + y_fill_min;
		
				return rrr;
				
			}else if(dimension == 2)
			{
			//	srand(count);
				double rrr = ((double)rand()) / ((double)(RAND_MAX));
				rrr = (rrr * (z_fill_max - z_fill_min)) + z_fill_min;
		
				return rrr;	
				
			}
			
	//	random_number = random_number + 3;
	}
/***************************************************************************************************************/

double random_coordinates_generator_bin(int dimension, int random_number)	//dimension defines whether it is to be done for x or y or z. 0 for x, 1 for y and 2 for z
	{
		
			random_number *= 233;
		srand(random_number); 			//Current time value as seed for random number generation//	
			
		if(dimension == 0)
			{
			//	srand(count);
				double rrr = ((double)rand()) / ((double)(RAND_MAX));
				rrr = (rrr * (x_fill_max_bin - x_fill_min_bin)) + x_fill_min_bin;
		
				return rrr;
				
			}else if(dimension == 1)
			{
			//	srand(count);
				double rrr = ((double)rand()) / ((double)(RAND_MAX));
				rrr = (rrr * (y_fill_max_bin - y_fill_min_bin)) + y_fill_min_bin;
		
				return rrr;
				
			}else if(dimension == 2)
			{
			//	srand(count);
				double rrr = ((double)rand()) / ((double)(RAND_MAX));
				rrr = (rrr * (z_fill_max_bin - z_fill_min_bin)) + z_fill_min_bin;
		
				return rrr;	
				
			}
			
	//	random_number = random_number + 3;
	}
/***************************************************************************************************************/