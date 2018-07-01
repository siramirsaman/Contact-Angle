/* Fluent Dynamic Contact Angle UDF: In Parallel
 * See https://arxiv.org/abs/1608.02856 for details on correlations
* Amirsaman Farrokhpanah
* The following gives a general idea on how the dynamic contact angle can be implemented in parallel
* in the UDF. It needs to be modified to match your problem.
* Disclaimer: Use at your own risk!
*/
#include "udf.h"

/* Global Variables:
* used to pass information between the two functions below
*/
double dynamic_contact_angle = 60.0; /*dynamic angle: initial value */
double contact_line_velocity = 0; /* capilary velocity */
FILE *file; /* output file, defined here for ease of use */
/* ....
* define other variables as needed
*/

/* This routine calculates the value of contact angle on each time step.
* It is hooked to "Adjust Functions" in Fleunt
*/
DEFINE_ADJUST(Contact_Angle_Update, domain)
{
	#if !RP_HOST    /*  only host process is involved */
	Thread *thread = Lookup_Thread(domain, 2); /* 2 is id number of the wall boundary where the angle is applied to, needs to be updated, you can find it from the boundary tab from Fluent user interface */
	Thread **pt = THREAD_SUB_THREADS(thread);
	cell_t cell;
	face_t f;


	/* calculate capilary velocity here.
	   Read the paper above to get an idea.
	   it needs to be the velocity of the moving contact point
	   not the velocity of fluid:
	   
	   contact_line_velocity = (current_position_of_triple_point - last_known_position_of_triple_point) / (the_time_between_these_two_measurements)
	
	   updating it every timestep will cause accuracy issues, need to jump every couple of time steps (using N_TIME)
	*/
	if (DESIRED_INTERVALS_BETWEEN_TWO_MEASUREMENTS)
	{
		begin_c_loop_all (cell,pt[1])

		{
			C_CENTROID(x,cell,pt[1]);

			if(C_VOF(cell,pt[1])>0.5)
			{
				C_CENTROID(x,cell,pt[1]); /* perform averaging, or other needed calculations on the boundary here */
			}
		}
		end_c_loop_all (cell,pt[1])


		#if RP_NODE     //in parallel mode, sums up the variables here, so values from parallel chuncks are added together
			integer_1 = PRF_GISUM1(integer_1); /* example parallel sum methods */
			array_1[0] = PRF_GRSUM1(array_1[0]);  /* example parallel sum methods */
			float_1 = PRF_GRSUM1(float_1);  /* example parallel sum methods */
		#endif

		/* update here:
			contact_velocity and dynamic_contact_angle
		*/

		if(I_AM_NODE_ZERO_P) /* in parallel */
		{
			if(DO_SOME_SAFETY_CHECKS_HERE)
			dynamic_contact_angle=60;

			file = fopen("output.txt", "a+"); /* write to file for debug purposes */
			fprintf(file, "Time=%f ontact_velocity=%f dynamic_contact_angle=%f \n", CURRENT_TIME, contact_velocity, dynamic_contact_angle);
			fclose(file);
		}

		old_time = CURRENT_TIME; /* Store old values here for next loop, need to be defined in global mode */
		contact_velocity_old = contact_velocity;
	}

}


#endif
}


DEFINE_PROFILE(Contact_Angle_Set_Profile,t,i)
{
	#if !RP_HOST      /*  only the host process is involved */
		face_t f;
		begin_f_loop(f,t)
		{
			F_PROFILE(f,t,i) = dynamic_contact_angle;
		}
		end_f_loop(f,t)
	#endif
}
