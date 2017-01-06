
/* See https://arxiv.org/abs/1608.02856 for details on correlations 
 * Amirsaman Farrokhpanah
 * The following gives a general idea on how the dynamic contact angle can be implemented
 * in the UDF. It needs to be modified to match your problem.
 * Disclaimer: Use at your own risk.
 *
 */

#include "udf.h"

double contact_line_position=0;
double time=0;
double theta_e_radian=(15.0*M_PI/180.0);
double dynamic_contact_angle=15.0;
int first_time=1;
double contact_velocity=0;
FILE *file;
double sum=0.0,R;


DEFINE_ADJUST(Contact_Angle_Update, domain)
{
	Thread *thread = Lookup_Thread(domain, 1);
	Thread **pt = THREAD_SUB_THREADS(thread); 
	cell_t cell;
	face_t f;
	real x[ND_ND];
	double max_x=0, f_Hoff_inverse, x_hoff,temp, Ca, volume;
	int n;

	sum=0.0;
	
	begin_c_loop_all (cell,pt[1])
	{
		if(C_VOF(cell,pt[1])!=0 /*&& ABS(x[1]-0.000125)<1e-6 && ABS(x[2]+0.003)<1e-6 */)
		{
			C_CENTROID(x,cell,pt[1]);
			//printf("## x=%f y=%f z=%f vof=%f\n",x[0],x[1],x[2],C_VOF(cell,pt[1]));
			if(x[0]>max_x)
			max_x=x[0];
		}
		sum+=C_VOF(cell,pt[1]);
	}
	end_c_loop_all (cell,pt[1])

	R=sqrt(sum*1.0e-08/M_PI);
	printf("max_x=%f R=%f contact_velocity=%f time=%f dynamic_contact_angle=%f sum=%f\n",max_x,R,contact_velocity,time,dynamic_contact_angle,sum);
 
	if(first_time==0)
	{
		contact_velocity= (R-contact_line_position)/(RP_Get_Real("flow-time")-time);
		contact_line_position = R;
		time=RP_Get_Real("flow-time"); 
		Ca = contact_velocity*1e-3/0.0728;
		temp= 0.5-0.5 * cos(theta_e_radian);
		temp= 0.5 * log( (1.0+temp)/(1.0-temp) );//atanh(0.5-0.5 * cos(theta_e_radian))
		f_Hoff_inverse = -(9.78546 * pow( temp,1.416430594900850))/(12.819 * pow(temp,1.41643)-100.0);
		x_hoff= Ca +  f_Hoff_inverse ;
		if(contact_velocity>=0)
			dynamic_contact_angle= acos( 1-2*tanh(5.16*pow((x_hoff/(1+1.31*pow(x_hoff,.99))),0.706)   ) ) *180.0/M_PI  ;
		else
			dynamic_contact_angle= 2*theta_e_radian*180.0/M_PI -acos( 1-2*tanh(5.16*pow((x_hoff/(1+1.31*pow(x_hoff,.99))),0.706)   ) ) *180.0/M_PI  ;
		file = fopen("file.txt", "a+");
		fprintf(file,"max_x=%f R=%f contact_velocity=%f time=%f dynamic_contact_angle=%f\n",max_x,R,contact_velocity,time,dynamic_contact_angle);
		fclose(file);
		
	 if(R>1e-18 && first_time==1)
	 {
		contact_line_position = R;
		time=RP_Get_Real("flow-time");
		first_time=0;
		file = fopen("file.txt", "a+");
		fprintf(file,"max_x=%f R=%f contact_velocity=%f time=%f\n",max_x,R,contact_velocity,time);
	 	fclose(file);
	}
}

DEFINE_PROFILE(Contact_Angle_Set_Profile,t,i)
{
	face_t f;
	begin_f_loop(f,t)
	{
		F_PROFILE(f,t,i) = dynamic_contact_angle;
    }
	end_f_loop(f,t)
}