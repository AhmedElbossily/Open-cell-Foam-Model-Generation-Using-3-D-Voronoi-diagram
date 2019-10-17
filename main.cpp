#include <iostream>
#include <math.h>
#include <stdlib.h> 
#include <sstream>
#include <fstream> 
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "CreateDxf.h"
#include "tetgen.h"

using namespace std;
double _pi=3.1415926535897932384626433832795;
double _x0=1;
double _x1=(1.0-_x0)/2.0;
double _min_angle=90;
double _max_angle=90;

double dist_two_points_3d(double x1,double y1,double z1,double x2,double y2,double z2);
double random_number();
void  mps_3d(double *points_x,double *points_y,double *points_z,int num_points,double r,int &num_itr,int initiation);
void  regular_points_seeds(double *points_x,double *points_y,double *points_z,double dist);
double max(double d_1,double d_2);
void trim(float &x1 ,float &y1,float &z1,float &x2,float &y2,float &z2,int &boundary,int &draw_check);
void draw_cells(double *draw_points, int num_draw_points,string &str);
void one_cell_draw(tetgenio &out,int &draw_check,string str);
void V_unite_cube_extraction(tetgenio &out,int &draw_check,double &total_distance, double * draw_points,int * draw_check_array,
							 int &num_draw_points,float &x1,float&y1,float&z1,float&x2,float&y2,float&z2,int &draw_index,
							 int &point_index);
void small_ligm_remove(tetgenio &out,int &draw_check,double beam_rad,int* ch_p_posi,int &p,int &num_small_lig);
void tet_gen(int num_points,double *points_x,double*points_y,double*points_z,double &beam_rad,double denesity,
			 double r_reg,double e,int initiation,double r,int &num_small_lig,int &num_draw_points,double * draw_points);
void tet_gen_honeycomb(int num_points,double *points_x,double*points_y,double*points_z,
					   double &beam_rad,double denesity,int &num_draw_points,double * draw_points);
double min_dist_req_3d(double V,int num_points);
double random_number( double start,double interval );
void  MPS_complete_regular_bycub(double *points_x,double *points_y,double *points_z,double x0, 
								 double y0, double z0,double dist, int &hit, double disk_raduis, int se );
int covered(double *child_cub,double *points_x,double *points_y,double *points_z, double disk_raduis,int hit);
void child_setup(double *child_cub,double x,double y,double z,double dist);
void pearants(int step,double *AAA,int hit,double *points_x,double *points_y,
			  double *points_z,double dist,double *empty_cubs,int &empty_index);
void  cube_corners(double *corner,double x,double y,double z,double dist, int level);
void edit_empty_cubs(double *empty_cubs,int &empty_index);
int seedstracking_lastparent(int &hit,double *points_x,double *points_y, double *points_z, double disk_raduis);
void  mps_3d_for_high_regularity(double *points_x,double *points_y,double *points_z,int num_points,double r,
								 int &num_itr,int initiation);

int main()
{	
	// creat array that holds voroni vertices 
	double * draw_points;
	draw_points= new double [60000000];
	int num_draw_points=0;
	
	//Create a parameter to check whether the user want to generate a honeycomb structure 1 or random voronoi diagram 
	int reg_type = 1;
	
	//the required number of seed inside the unite cube 
	int num_points=559;
	int num_itr=1;
	double* points_x, * points_y, * points_z;
	
	// Create variables that holds Vornoi ligaments characteristics
	double beam_rad=0.0;
	int num_small_lig=0;
	double denesity;

	// calculate the required distance to generate complete regulare voronoi structure for "num_points" seeds in unite cube "_x0*_x0*_x0"
	double r_reg;
	r_reg=min_dist_req_3d(_x0*_x0*_x0,num_points);

	//in case of user want to generate complete regulare structure "honeycomb"
	if(reg_type==1)
	{
		// set the the required density of model
		denesity=0.03;

		// Create arrays holds seeds coordinates 
		points_x = new double [num_points];
		points_y = new double [num_points];
		points_z = new double [num_points];

		// generate regulare distribution seeds 
		regular_points_seeds(points_x,points_y,points_z,r_reg);

		// Creating Vornoi vertices for honeycomb structure using TetGen library 
		tet_gen_honeycomb(num_points,points_x,points_y,points_z,beam_rad, denesity,num_draw_points,draw_points);

		//set the destination for storing DXF file 
		string str="D:\\master\\results\\1reg\\00_Drawing.dxf";

		// Draw Honeycomb structure in dxf file
		draw_cells(draw_points,num_draw_points,str);

		//display voronoi ligaments radiaus in inch and in m
		cout<<"beam_rad	:"<<beam_rad<<"	inch"<<endl;
		cout<<"beam_diameter	:"<<beam_rad*.0254*2<<"	m"<<endl;
	}
	//in case of user want to generete random distributing voronoi 
	else 
	{
		// random number geration function initiation
		int initiation=0;

		//Regularity parameter 0 % complete irrigulare structure and 100% is honeycomb structure
		double  e = 0.9;
		
		//these values are evaluated from running the code before. the code converge more quicly when srand() functions initiated with these values
		int good_initiation[20]={2,28,29,37,45,56,70,95,104,126,162,178,205,229,243,268,309,322,326,357};
		
		//create arrays holps seeds coordinates 
		points_x = new double [num_points];
		points_y = new double [num_points];
		points_z = new double [num_points];

		//location where the Dxf file and statistics .txt will be stored
		string str="D:\\master\\results\\9reg\\00\\00_Drawing.dxf";
		string str_2="D:\\master\\results\\9reg\\00\\00_beam_radius.txt";
		string file_num;
		stringstream convert; // stringstream used for the conversion.
		
		// this variable will hold density itration 
		int e_itra=0;
		double m_rad_1=0.0,m_rad_3=0.0,m_rad_5=0.0,m_rad_7=0.0,m_rad_9=0.0;
		
		//we need to generate 20 model automatically and store them in different distinations 
		while (initiation<20)
		{
			e_itra=0;
			num_itr=0;

			//use the modified MPS to generate vornoi with regularity 90%
			mps_3d_for_high_regularity(points_x,points_y,points_z,num_points,e*r_reg,num_itr,good_initiation[initiation]);
			//mps_3d(points_x,points_y,points_z, num_points, e*r_reg,num_itr,initiation);
			
			while(e_itra<5)
			{
				if(e_itra==0) denesity=.01;
				if(e_itra==1) denesity=.03;
				if(e_itra==2) denesity=.05;
				if(e_itra==3) denesity=.07;
				if(e_itra==4) denesity=.09;

				convert << initiation;		     // add   the value of Number to the characters in the stream
				file_num = convert.str();		//  set   Result to the content of the stream
				
				//change location of the output file based on random points function initiation number
				if (initiation<10)
				{
					str.replace(24,1,file_num);
					str_2.replace(24,1,file_num);
					str.replace(27,1,file_num);
					str_2.replace(27,1,file_num);
				}
				if (initiation>=10)	
				{
					str.replace(23,2,file_num);
					str_2.replace(23,2,file_num);
					str.replace(26,2,file_num);
					str_2.replace(26,2,file_num);
				}
				cout<<str<<endl;
				convert.str("");
				file_num.clear();

				// Get Vornoi vertices using TetGen library
				tet_gen(num_points,points_x,points_y,points_z,beam_rad,denesity,r_reg,e,initiation,e*r_reg,num_small_lig, num_draw_points,draw_points);

				// sum the beam diameter for every density and at the end calculate the mean for 20 model
				if(e_itra==0) m_rad_1=m_rad_1+beam_rad;
				if(e_itra==1) m_rad_3=m_rad_3+beam_rad;
				if(e_itra==2) m_rad_5=m_rad_5+beam_rad;
				if(e_itra==3) m_rad_7=m_rad_7+beam_rad;
				if(e_itra==4)m_rad_9=m_rad_9+beam_rad;

				// Draw Voronoi Seeds for the first itration of densities as the structure will not change by changing densities 
				// only Vornoi legaments diameter will change 
				if(e_itra==0) draw_cells(draw_points,num_draw_points,str);

				if(e_itra==0){ofstream beam_radius (str_2, ios::out);
				beam_radius<<"denesity	"<<denesity<<endl;
				beam_radius<<"beam_rad	"<<beam_rad<<"	inch"<<endl;
				beam_radius<<"num_small_lig	"<<num_small_lig<<endl;
				beam_radius<<"regularity	"<<e<<endl;
				beam_radius<<"***********************"<<endl;}

				else	{ofstream beam_radius (str_2, ios::app);
				beam_radius<<"denesity	"<<denesity<<endl;
				beam_radius<<"beam_rad	"<<beam_rad<<"	inch"<<endl;
				beam_radius<<"num_small_lig	"<<num_small_lig<<endl;
				beam_radius<<"regularity	"<<e<<endl;
				beam_radius<<"***********************"<<endl;}

				e_itra++;
			}
			initiation++;
		}

	
		// calculate the mean of Voronoi ligaments for 20 model and store them in dxf file 
		ofstream beam_radius ("D:\\master\\results\\9reg\\00_beam_radius.txt", ios::out);
		beam_radius<<"regularity	"<<e<<endl;
		beam_radius<<"m_rad_1	"<<m_rad_1/initiation<<endl;
		beam_radius<<"m_rad_3	"<<m_rad_3/initiation<<endl;
		beam_radius<<"m_rad_5	"<<m_rad_5/initiation<<endl;
		beam_radius<<"m_rad_7	"<<m_rad_7/initiation<<endl;
		beam_radius<<"m_rad_9	"<<m_rad_9/initiation<<endl;

	}

	//system("pause");

	return 0;

}

// return the distance between two points in 3d
double dist_two_points_3d(double x1,double y1,double z1,double x2,double y2,double z2)
{
	double ans=sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
	return (ans);
}

//return the random points between _X1 and _X0 interval 
double random_number()
{
	return (_x1+((double) rand() / (RAND_MAX))*_x0);
}

/* Modified MPS for 3d. this function will return when random points reach the desired points "num_points" 
@param
	num_points.....required number of hit pints 
	r.............. minimum distance 
*/
void  mps_3d(double *points_x,double *points_y,double *points_z,int num_points,double r,int &num_itr,int initiation)
{
	// initiate random point function 
	srand(initiation);

	double x,y,z;

	x=random_number();
	y=random_number();
	z=random_number();

	points_x[0]=x;
	points_y[0]=y;
	points_z[0]=z;

	int hit=0;

	while(true)
	{
		// creat another random points 
		x=random_number();
		y=random_number();
		z=random_number();
		double min_distance=1000;
		double distance=0;

		for (int j = 0; j < hit+1; j++)
		{
			distance=dist_two_points_3d(points_x[j],points_y[j],points_z[j],x,y,z);
			if (distance<min_distance) min_distance=distance;

		}

		if (min_distance>r)
		{
			hit++;
			points_x[hit]=x;
			points_y[hit]=y;
			points_z[hit]=z;
		}

		num_itr++;

		if(hit==num_points-1)	break;
	}
}

//generating regulare distriputing points in 3d cubde 
void  regular_points_seeds(double *points_x,double *points_y,double *points_z,double dist)
{
	const int cellperlength=6;
	dist=(1.0/cellperlength);
	double AAA[cellperlength+1],BBB[cellperlength+2];
	double initial_dist_A=(1-cellperlength*dist)/2;
	double initial_dist_B=initial_dist_A+(dist/2);

	for (int i = 0; i < cellperlength+1; i++)		AAA[i]=initial_dist_A+i*dist;
	for (int i = 0; i < cellperlength; i++)			BBB[i]=initial_dist_B+i*dist;

	double x,y,z;
	int index=0;
	for (int i = 0; i < cellperlength+1; i++)
	{
		x= AAA[i];
		for (int j = 0; j < cellperlength+1; j++)
		{
			y= AAA[j];
			for (int k = 0; k < cellperlength+1; k++)
			{
				z= AAA[k];
				points_x[index]=x;
				points_y[index]=y;
				points_z[index]=z;
				index++;
			}

		}
	}

	for (int i = 0; i < cellperlength; i++)
	{
		x= BBB[i];
		for (int j = 0; j < cellperlength; j++)
		{
			y= BBB[j];
			for (int k = 0; k < cellperlength; k++)
			{
				z= BBB[k];
				points_x[index]=x;
				points_y[index]=y;
				points_z[index]=z;
				index++;
			}
		}
	}

}

//return max of two numbers
double max(double d_1,double d_2)
{
	if (d_1>d_2) return d_1;
	else return d_2;

}

//Trim Vornoi boundary ligaments by cube planes 
void trim(float &x1 ,float &y1,float &z1,float &x2,float &y2,float &z2,int &boundary,int &draw_check)
{
	float d=0;
	draw_check=1;
	if (y2>1 && x1>=0 && x1<=1 && y1>=0 && y1<=1 && z1>=0 && z1<=1)
	{
		float px=0.0,py=1.0,pz=1.0;
		d=(py-y1)/(y2-y1);
		x2=x1+(x2-x1)*d;
		y2=y1+(y2-y1)*d;
		z2=z1+(z2-z1)*d;
		boundary=1;
		draw_check=0;
	}
	if (y1>1 && x2>=0 && x2<=1 && y2>=0 && y2<=1 && z2>=0 && z2<=1)
	{	
		float px=0.0,py=1.0,pz=1.0;
		d=(py-y2)/(y1-y2);
		x1=x2+(x1-x2)*d;
		y1=y2+(y1-y2)*d;
		z1=z2+(z1-z2)*d;
		boundary=1;
		draw_check=0;
	}		
	if(y2<0 && x1>=0 && x1<=1 && y1>=0 && y1<=1 && z1>=0 && z1<=1)
	{
		float px=0.0,py=1.0,pz=1.0;
		d=(y1)/(y1-y2);
		x2=x1+(x2-x1)*d;
		y2=y1+(y2-y1)*d;
		z2=z1+(z2-z1)*d;
		boundary=1;
		draw_check=0;
	}
	if (y1<0 && x2>=0 && x2<=1 && y2>=0 && y2<=1 && z2>=0 && z2<=1)
	{	
		float px=0.0,py=1.0,pz=1.0;
		d=(y2)/(y2-y1);
		x1=x2+(x1-x2)*d;
		y1=y2+(y1-y2)*d;
		z1=z2+(z1-z2)*d;
		boundary=1;
		draw_check=0;
	}
	// triming with plan x
	if (x2>1 && x1>=0 && x1<=1 && y1>=0 && y1<=1 && z1>=0 && z1<=1)
	{
		d=(1.0-x1)/(x2-x1);
		x2=x1+(x2-x1)*d;
		y2=y1+(y2-y1)*d;
		z2=z1+(z2-z1)*d;
		boundary=1;
		draw_check=0;
	}
	if (x1>1 && x2>=0 && x2<=1 && y2>=0 && y2<=1 && z2>=0 && z2<=1)
	{	
		d=(1.0-x2)/(x1-x2);
		x1=x2+(x1-x2)*d;
		y1=y2+(y1-y2)*d;
		z1=z2+(z1-z2)*d;
		boundary=1;
		draw_check=0;
	}
	//trimming with plan -x
	if(x2<0 && x1>=0  && x1<=1 && y1>=0 && y1<=1 && z1>=0 && z1<=1)
	{
		d=x1/(x1-x2);
		x2=x1+(x2-x1)*d;
		y2=y1+(y2-y1)*d;
		z2=z1+(z2-z1)*d;
		boundary=1;
		draw_check=0;
	}
	if (x1<0 && x2>=0 && x2<=1 && y2>=0 && y2<=1 && z2>=0 && z2<=1)
	{	
		d=x2/(x2-x1);
		x1=x2+(x1-x2)*d;
		y1=y2+(y1-y2)*d;
		z1=z2+(z1-z2)*d;
		boundary=1;
		draw_check=0;
	}

	// triming with plan z
	if (z2>1 && x1>=0 && x1<=1 && y1>=0 && y1<=1 && z1>=0 && z1<=1)
	{
		d=(1.0-z1)/(z2-z1);
		x2=x1+(x2-x1)*d;
		y2=y1+(y2-y1)*d;
		z2=z1+(z2-z1)*d;
		boundary=1;
		draw_check=1;
	}
	if (z1>1 && x2>=0 && x2<=1 && y2>=0 && y2<=1 && z2>=0 && z2<=1)
	{	
		d=(1.0-z2)/(z1-z2);
		x1=x2+(x1-x2)*d;
		y1=y2+(y1-y2)*d;
		z1=z2+(z1-z2)*d;
		boundary=1;
		draw_check=1;
	}
	//trimming with plan -z
	if(z2<0 && x1>=0 && x1<=1 && y1>=0 && y1<=1 && z1>=0 && z1<=1)
	{
		d=z1/(z1-z2);
		x2=x1+(x2-x1)*d;
		y2=y1+(y2-y1)*d;
		z2=z1+(z2-z1)*d;
		boundary=1;
		draw_check=1;
	}
	if (z1<0 && x2>=0 && x2<=1 && y2>=0 && y2<=1 && z2>=0 && z2<=1)
	{	
		d=z2/(z2-z1);
		x1=x2+(x1-x2)*d;
		y1=y2+(y1-y2)*d;
		z1=z2+(z1-z2)*d;
		boundary=1;
		draw_check=1;
	}
}

// Draw voronoi cells in dxf file 
void draw_cells(double *draw_points, int num_draw_points,string &str)
{
	CreateDxf Draw;
	Draw.DxfBegin( str);
	float x1,y1,z1,x2,y2,z2;
	for (int i = 0; i < num_draw_points; i++)
	{
		x1=draw_points[i*6];
		y1=draw_points[i*6+1];
		z1=draw_points[i*6+2];
		x2=draw_points[i*6+3];
		y2=draw_points[i*6+4];
		z2=draw_points[i*6+5];
		Draw.line(x1,y1,z1, x2,y2,z2, str);
	}
	Draw.DxfEnd(str);
}

// Draw one voronoi cell
void one_cell_draw(tetgenio &out,int &draw_check,string str)
{
	float x1,y1,z1,x2,y2,z2;

	CreateDxf Draw;
	Draw.DxfBegin(str);

	int cell_num=150;
	int num_faces=out.vcelllist[cell_num][0];// number of faces
	int face;
	int num_edg;
	int edg;
	for (int i = 1; i <= num_faces; i++)
	{
		face=out.vcelllist[cell_num][i];

		if(face!=-1)
		{
			num_edg=out.vfacetlist[face].elist[0];// number of edges

			for (int j = 1; j <= num_edg; j++)
			{
				edg=out.vfacetlist[face].elist[j];
				if(edg!=-1)
				{
					x1=out.vpointlist[3*out.vedgelist[edg].v1];
					y1=out.vpointlist[3*out.vedgelist[edg].v1+1];
					z1=out.vpointlist[3*out.vedgelist[edg].v1+2];


					if (out.vedgelist[edg].v2==-1.0)
					{
						x2=out.vedgelist[edg].vnormal[0];
						y2=out.vedgelist[edg].vnormal[1];
						z2=out.vedgelist[edg].vnormal[2];
					}
					else
					{
						x2=out.vpointlist[3*out.vedgelist[edg].v2];
						y2=out.vpointlist[3*out.vedgelist[edg].v2+1];
						z2=out.vpointlist[3*out.vedgelist[edg].v2+2];
					}

					if((x1>=1.0||y1>=1.0||z1>=1.0||x1<=0.0||y1<=0.0||z1<=0.0)&&(x2>=1.0||y2>=1.0||z2>=1.0||x2<=0.0||y2<=0.0||z2<=0.0)){}
					else  if ((out.vedgelist[edg].v2==-1.0)&&(x1>=1.0||y1>=1.0||z1>=1.0||x1<=0.0||y1<=0.0||z1<=0.0))	{}
					else
					{
						int boundary=0;
						trim( x1, y1,z1, x2, y2,z2,boundary,draw_check);

						if(dist_two_points_3d( x1*.25, y1*.25,z1*.25, x2*.25, y2*.25,z2*.25)<_x0/2)
						{
							Draw.line(x1*.25,y1*.25,z1*.25, x2*.25,y2*.25,z2*.25,str);

						}
					}
				}

			}

		}
	}
	x1=out.pointlist[3*cell_num];
	y1=out.pointlist[3*cell_num+1];
	z1=out.pointlist[3*cell_num+2];

	Draw.line(x1*.25,y1*.25,z1*.25, x2*.25,y2*.25,z2*.25,str);
	Draw.DxfEnd(str);
}

void V_unite_cube_extraction(tetgenio &out,int &draw_check,double &total_distance,double * draw_points,int * draw_check_array,int &num_draw_points,
					 float &x1,float&y1,float&z1,float&x2,float&y2,float&z2,int &draw_index,int &point_index)
{
	total_distance=0;
	draw_index=0;
	point_index=0;
	num_draw_points=0;
	double min_strut=100;
	for (int i = 0; i < out.numberofvedges; i++)
	{
		x1=out.vpointlist[3*out.vedgelist[i].v1];
		y1=out.vpointlist[3*out.vedgelist[i].v1+1];
		z1=out.vpointlist[3*out.vedgelist[i].v1+2];

		if (out.vedgelist[i].v2==-1.0)
		{
			x2=out.vedgelist[i].vnormal[0];
			y2=out.vedgelist[i].vnormal[1];
			z2=out.vedgelist[i].vnormal[2];
			//cout<<"edge "<<i<<endl;
		}
		else
		{
			x2=out.vpointlist[3*out.vedgelist[i].v2];
			y2=out.vpointlist[3*out.vedgelist[i].v2+1];
			z2=out.vpointlist[3*out.vedgelist[i].v2+2];
		}

		// check if the line located inside the unit cube
		if((x1>1.0||y1>1.0||z1>1.0||x1<0.0||y1<0.0||z1<0.0)&&(x2>1.0||y2>1.0||z2>1.0||x2<0.0||y2<0.0||z2<0.0)){}
		else  if ((out.vedgelist[i].v2==-1.0)&&(x1>=1.0||y1>=1.0||z1>=1.0||x1<=0.0||y1<=0.0||z1<=0.0))	{}
		else
		{
			int boundary=0;
			trim( x1, y1,z1, x2, y2,z2,boundary,draw_check);

			if(dist_two_points_3d( x1, y1,z1, x2, y2,z2)<_x0/2)
			{
				if (dist_two_points_3d( x1, y1,z1, x2, y2,z2)*0.25<min_strut)
				{
					min_strut=dist_two_points_3d( x1, y1,z1, x2, y2,z2)*0.25;
				}
				total_distance=total_distance + dist_two_points_3d( x1, y1,z1, x2, y2,z2)*0.25;
				draw_points[point_index++]=x1*0.25;
				draw_points[point_index++]=y1*0.25;
				draw_points[point_index++]=z1*0.25;
				draw_points[point_index++]=x2*0.25;
				draw_points[point_index++]=y2*0.25;
				draw_points[point_index++]=z2*0.25;
				num_draw_points++;
				draw_check_array[draw_index++]=draw_check;
			}
		}			
	}

}

void small_ligm_remove(tetgenio &out,int &draw_check,double beam_rad,int* ch_p_posi,int &p,int &num_small_lig)
{
	p=0;
	float x1,y1,z1,x2,y2,z2;
	num_small_lig=0;
	for (int i = 0; i < out.numberofvedges; i++)
	{
		x1=out.vpointlist[3*out.vedgelist[i].v1];
		y1=out.vpointlist[3*out.vedgelist[i].v1+1];
		z1=out.vpointlist[3*out.vedgelist[i].v1+2];

		if (out.vedgelist[i].v2==-1.0)
		{
			x2=out.vedgelist[i].vnormal[0];
			y2=out.vedgelist[i].vnormal[1];
			z2=out.vedgelist[i].vnormal[2];
		}
		else
		{
			x2=out.vpointlist[3*out.vedgelist[i].v2];
			y2=out.vpointlist[3*out.vedgelist[i].v2+1];
			z2=out.vpointlist[3*out.vedgelist[i].v2+2];
		}

		if((x1>=1.0||y1>=1.0||z1>=1.0||x1<=0.0||y1<=0.0||z1<=0.0)&&(x2>=1.0||y2>=1.0||z2>=1.0||x2<=0.0||y2<=0.0||z2<=0.0)){}
		else  if ((out.vedgelist[i].v2==-1.0)&&(x1>=1.0||y1>=1.0||z1>=1.0||x1<=0.0||y1<=0.0||z1<=0.0))	{}
		else
		{
			int boundary=0;
			trim( x1, y1,z1, x2, y2,z2,boundary,draw_check);

			if(dist_two_points_3d( x1, y1,z1, x2, y2,z2)>0 && dist_two_points_3d( x1, y1,z1, x2, y2,z2)<=1.5*beam_rad*2)
			{
				for (int j = 0; j < out.numberofvfacets; j++)
				{
					for (int k = 1; k <= out.vfacetlist[j].elist[0]; k++)
					{
						if (out.vfacetlist[j].elist[k]==i)
						{
							int ch_check_c1=1;
							int ch_check_c2=1;

							for (int m = 0; m < p; m++)
							{
								if(ch_p_posi[m]==out.vfacetlist[j].c1) 
									ch_check_c1=0;
							}
							for (int m = 0; m < p; m++)
							{
								if(ch_p_posi[m]==out.vfacetlist[j].c2) ch_check_c2=0;
							}

							if(ch_check_c1)
							{
								ch_p_posi[p]=out.vfacetlist[j].c1;
								p++;
							}

							if(ch_check_c2)
							{
								ch_p_posi[p]=out.vfacetlist[j].c2;
								p++;
							}

						}
					}

				}
				num_small_lig++;
			}

		}			
	}

	std::cout<<endl<<"num_small_lig "<<num_small_lig<<endl;
	std::cout<<endl<<"p "<<p<<endl;
}

void tet_gen(int num_points,double *points_x,double*points_y,double*points_z,double &beam_rad,double denesity,double r_reg,double e,int initiation,double r,int &num_small_lig,int &num_draw_points,double * draw_points)
{
	num_small_lig=0;
	double total_distance=0.0;
	tetgenio out;
	tetgenbehavior A ;
	tetgenio in;

	A.voroout = 1;

	int draw_check=1;
	float x1,y1,z1,x2,y2,z2;
	int* ch_p_posi ;

	int * draw_check_array;
	int p=0; // number of removed points.
	int min_p=10000;
	int num_itr=0;
	ch_p_posi=new int[num_points] ;

	draw_check_array= new int [100000];
	num_draw_points=0;
	int draw_index=0;
	int point_index=0;
	int best_p_initiation=0;
	int best_lig=0;
	int min_num_small_lig=10000;
	while (true)
	{
		in.numberofpoints=num_points;
		in.pointlist = new double[in.numberofpoints * 3];
		in.pointmarkerlist=new int[in.numberofpoints];
		for (int i = 0; i < num_points; i++)
		{
			in.pointmarkerlist[i]=0;
			in.pointlist[i * 3] =points_x[i];
			in.pointlist[i * 3 + 1] = points_y[i];
			in.pointlist[i * 3 + 2] = points_z[i];
		}

		tetrahedralize(&A,&in,&out);

		V_unite_cube_extraction(out,draw_check,total_distance,draw_points,draw_check_array,num_draw_points,x1,y1,z1,x2,y2,z2, draw_index, point_index);

		beam_rad=sqrt((0.25*0.25*0.25*denesity)/(total_distance*_pi));

		small_ligm_remove( out,draw_check, beam_rad,ch_p_posi,p,num_small_lig);

		// Remove this break if you want to iterate throw srand() initiation to get the best initiation number that not has short ligaments 
		break;

		if(p<min_p)
		{
			min_p=p;
			best_p_initiation=initiation;
		}

		if(num_small_lig<min_num_small_lig)
		{
			min_num_small_lig=num_small_lig;
			best_lig=initiation;
		}

		initiation++;
		mps_3d(points_x,points_y,points_z, num_points, e*r_reg,num_itr,initiation);
		in.deinitialize();
		out.deinitialize();		
	}
}

// Creating Vornoi vertices for honeycomb structure using TetGen library 
void tet_gen_honeycomb(int num_points,double *points_x,double*points_y,double*points_z,double &beam_rad,double denesity,int &num_draw_points,double * draw_points)
{
	double total_distance=0.0;
	tetgenio out;
	tetgenbehavior A ;
	tetgenio in;
	int draw_check=1;
	int * draw_check_array;
	draw_check_array= new int [10000];

	A.voroout = 1;

	float x1,y1,z1,x2,y2,z2;
	int* ch_p_posi ;

	num_draw_points=0;
	int draw_index=0;
	int point_index=0;

	in.numberofpoints=num_points;
	in.pointlist = new double[in.numberofpoints * 3];
	in.pointmarkerlist=new int[in.numberofpoints];
	for (int i = 0; i < num_points; i++)
	{
		in.pointmarkerlist[i]=0;
		in.pointlist[i * 3] =points_x[i];
		in.pointlist[i * 3 + 1] = points_y[i];
		in.pointlist[i * 3 + 2] = points_z[i];
	}

	tetrahedralize(&A,&in,&out);

	V_unite_cube_extraction(out,draw_check,total_distance,draw_points,draw_check_array,num_draw_points,x1,y1,z1,x2,y2,z2, draw_index, point_index);

	beam_rad=sqrt((0.25*0.25*0.25*denesity)/(total_distance*_pi));

	//	in.deinitialize();
	//	out.deinitialize();

}

// min distance required between seeds in 3d
double min_dist_req_3d(double V,int num_points)
{	
	double n=V/(sqrt(2.0)*num_points);
	return((sqrt(6.0)/2.0)*pow(n,1/3.));
}

//generating random numders in interval 
double random_number( double start,double interval )
{
	return (start+((double) rand() / (RAND_MAX))*interval);
}

//MPS inside child cube
void  MPS_complete_regular_bycub(double *points_x,double *points_y,double *points_z,double x0, 
								 double y0, double z0,double dist, int &hit, double disk_raduis, int se )
{
	double x,y,z;
	int itration=0;
	while(true)
	{
		x=random_number(x0,dist);
		y=random_number(y0,dist);
		z=random_number(z0,dist);
		double min_distance=1000;
		double distance=0;

		for (int j = 0; j <= hit; j++)
		{
			distance=dist_two_points_3d(points_x[j],points_y[j],points_z[j],x,y,z);

			if (distance<min_distance) min_distance=distance;

		}

		if (min_distance>=disk_raduis)
		{
			hit++;
			points_x[hit]=x;
			points_y[hit]=y;
			points_z[hit]=z;
			break;
		}

		itration++;
		if (se==0)
		{
			if(itration==100)break;
		}
		if(itration==1000)break;
	}

}

//check child cube is covered or not 
int covered(double *child_cub,double *points_x,double *points_y,double *points_z, double disk_raduis,int hit)
{
	double distance;
	double max_distance;
	for (int m = 0; m <= hit; m++)
	{
		max_distance=0.0;
		for (int i = 0; i < 24; i=i+3)
		{
			distance=dist_two_points_3d(points_x[m],points_y[m],points_z[m],child_cub[i],child_cub[i+1],child_cub[i+2]);
			if (distance>max_distance) max_distance=distance;
		}
		if (max_distance < disk_raduis) return(1);
	}
	return(0);
}

//set up child cube for Simple MPS by Implicit Quad-Trees method
void child_setup(double *child_cub,double x,double y,double z,double dist)
{
	child_cub[0]=x;				child_cub[1]=y;				child_cub[2]=z;

	child_cub[3]=x+dist;		child_cub[4]=y;				child_cub[5]=z;

	child_cub[6]=x+dist;		child_cub[7]=y+dist;		child_cub[8]=z;

	child_cub[9]=x;				child_cub[10]=y+dist;		child_cub[11]=z;


	child_cub[12]=x;			child_cub[13]=y;			child_cub[14]=z+dist;

	child_cub[15]=x+dist;		child_cub[16]=y;			child_cub[17]=z+dist;

	child_cub[18]=x+dist;		child_cub[19]=y+dist;		child_cub[20]=z+dist;

	child_cub[21]=x;			child_cub[22]=y+dist;		child_cub[23]=z+dist;

	for (int i = 0; i < 24; i++)
	{
		if (child_cub[i]>1.0)	child_cub[i]=1.0;
	}

}

//create parent cubes and check if it empty or not for Simple MPS by Implicit Quad-Trees method
void pearants(int step,double *AAA,int hit,double *points_x,double *points_y,
			  double *points_z,double dist,double *empty_cubs,int &empty_index)
{
	double x,y,z;
	int check;
	empty_index=0;
	for (int i = 0; i < step; i++)
	{
		x=AAA[i];
		for (int j = 0; j < step; j++)
		{
			y=AAA[j];
			for (int k = 0; k < step; k++)
			{
				z=AAA[k];
				check=0;
				for (int m = 0; m <= hit; m++)
				{
					if (points_x[m]>= x &&points_x[m]<= x+dist)
					{
						if (points_y[m]>=y &&points_y[m]<= y+dist)
						{
							if (points_z[m]>=z &&points_z[m]<= z+dist)
							{
								check=1;
								break;
							}

						}

					}

				}
				if (check==0)
				{
					empty_cubs[empty_index++]=AAA[i];
					empty_cubs[empty_index++]=AAA[j];
					empty_cubs[empty_index++]=AAA[k];
				}
			}

		}

	}
}

// evaluate cube cornors for Simple MPS by Implicit Quad-Trees method
void  cube_corners(double *corner,double x,double y,double z,double dist, int level)
{
	int index=0;
	for (int i = 0; i < level; i++)
	{
		for (int j = 0; j < level; j++)
		{
			for (int k = 0; k < level; k++)
			{
				corner[index++]=x+dist*i;
				if(corner[index-1]>1)	corner[index-1]=1;

				corner[index++]=y+dist*j;
				if(corner[index-1]>1)	corner[index-1]=1;

				corner[index++]=z+dist*k;
				if(corner[index-1]>1)	corner[index-1]=1;
			}

		}
	}

}

void edit_empty_cubs(double *empty_cubs,int &empty_index)
{
	for (int i = 0; i < empty_index; i=i+3)
	{
		if (empty_cubs[i]==1000.0)
		{
			empty_index=empty_index-3;	
			for (int j = i; j < empty_index; j++)
			{
				empty_cubs[j]=empty_cubs[j+3];

			}
			i=i-3;
		}
	}
}

//tracking empty cube and insert dart inside it 
int seedstracking_lastparent(int &hit,double *points_x,double *points_y, double *points_z, double disk_raduis)
{
	
	int check=0;
	int empty_index=0;
	int level=512;
	int check_hit,check_cover,stop_k,n_cub_points,n_cub,step;
	double *empty_cubs,*AAA;
	double child_cub[24], corner[24];
	double x,y,z, splite_dist,dist;
	
	//determine dimintion of the cube inside sphere dart  
	dist=sqrt((disk_raduis*disk_raduis)/2.0);
	step=int(1/dist)+1;

	//determine number of cubes which will cover the domain
	n_cub_points=(step+1)*(step+1)*(step+1);
	n_cub=(step)*(step)*(step);
	AAA= new double [step+1];

	AAA[step]=1.0;
	for (int i = 0; i < step; i++)		AAA[i]=i*dist;

	empty_cubs= new double [n_cub*3];

	// create parent small cubes and check if it occupied with dart or not and evalute empty cubes  
	pearants(step,AAA, hit,points_x,points_y,points_z, dist,empty_cubs,empty_index);

	// iterate over the empty cubes 
	for (int i = 0; i < empty_index; i=i+3)
	{
		//get the empty cube
		x=empty_cubs[i];y=empty_cubs[i+1];z=empty_cubs[i+2];

		//divid cube length by 2
		splite_dist=dist/2;

		//create child cube from parent cube
		child_setup(corner,x,y,z,splite_dist);

		for (int k = 0; k < 2*2*2*3; k=k+3)
		{
			stop_k=0;
			while (true)
			{
				//get the child coordinates 
				child_setup(child_cub,corner[k],corner[k+1],corner[k+2],splite_dist);

				//check the child is coverd or not
				check_cover=covered(child_cub,points_x,points_y,points_z,disk_raduis,hit);

				//break if the cube is coverd 
				if(check_cover) break;

				//if not covert throw darts inside it 
				if(!check_cover)
				{
					check_hit=hit;

					MPS_complete_regular_bycub(points_x,points_y,points_z,child_cub[0],child_cub[1],child_cub[2]
					,splite_dist, hit,disk_raduis,0);
					if(hit>check_hit)	{stop_k=1;break;}

					// create another child
					x=child_cub[0];y=child_cub[1];z=child_cub[2];
					k=0;
					splite_dist=splite_dist/2;
					child_setup(corner,x,y,z,splite_dist);
				}
			}
			if(stop_k)	break;
		}
		
		if(hit==558) break;
	}

	// iterate over the empty cubes from reverse 
	for (int i = 0; i < empty_index; i=i+3)
	{
		x=empty_cubs[i];y=empty_cubs[i+1];z=empty_cubs[i+2];
		splite_dist=dist/2;
		child_setup(corner,x,y,z,splite_dist);
		for (int k = 0; k < 2*2*2*3; k=k+3)
		{
			stop_k=0;
			while (true)
			{
				//get the child coordinates 
				child_setup(child_cub,corner[23-(k+2)],corner[23-(k+1)],corner[23-k],splite_dist);
				
				//check the child is coverd or not
				check_cover=covered(child_cub,points_x,points_y,points_z,disk_raduis,hit);

				if(check_cover) break;
				if(!check_cover)
				{
					check_hit=hit;

					MPS_complete_regular_bycub(points_x,points_y,points_z,child_cub[0],child_cub[1],child_cub[2]
					,splite_dist, hit,disk_raduis,0);
					if(hit>check_hit)	{stop_k=1;break;}

					// create another child
					x=child_cub[0];y=child_cub[1];z=child_cub[2];
					k=0;
					splite_dist=splite_dist/2;
					child_setup(corner,x,y,z,splite_dist);
				}
			}
			if(stop_k)	break;
		}
		
		if(hit==558) break;
	}

	//create matrix of cubes according to level of iteration 
	pearants(step,AAA, hit,points_x,points_y,points_z, dist,empty_cubs,empty_index);
	double *corner2;
	int check_cover_withot_point=0;
	int uu;
	for (int level = 2; level < 64; level=level*2)
	{	
		
		corner2= new double[level*level*level*3];
		uu=0;
		for (int i = 0; i < empty_index; i=i+3)
		{
			splite_dist=dist/level;
			x=empty_cubs[i];y=empty_cubs[i+1];z=empty_cubs[i+2];
			cube_corners(corner2, x,y,z,splite_dist, level);
			check_cover_withot_point=0;
			for (int j = 0; j < (level*level*level*3); j=j+3)
			{
				child_setup(child_cub,corner2[j],corner2[j+1],corner2[j+2],splite_dist);

				if(!covered(child_cub,points_x,points_y,points_z,disk_raduis,hit))
				{
					check_hit=hit;
					check_cover_withot_point=1;
					MPS_complete_regular_bycub(points_x,points_y,points_z,child_cub[0],child_cub[1],child_cub[2]
					,splite_dist, hit,disk_raduis,0);
					if(hit>check_hit)	break;
				}
				if(hit==558)
				{
					delete(empty_cubs);
					delete(AAA);
					delete(corner2);
					return(1);
				}
			}
			if(check_cover_withot_point==0)
			{
				empty_cubs[i]=1000.0;
				uu++;
				if ((empty_index/3)-uu<(558-hit))
				{
					delete(empty_cubs);
					delete(AAA);
					delete(corner2);
					return(0);

				}
			}
		}
		edit_empty_cubs(empty_cubs,empty_index);

		delete(corner2);
	}
	delete(AAA);
	delete(empty_cubs);
	return(0);

}

//Modified MPS for high regularity of seeds 
void  mps_3d_for_high_regularity(double *points_x,double *points_y,double *points_z,int num_points,double r,int &num_itr,int initiation)
{
	// Initiate srand()
	srand(initiation);

	double x,y,z;

	x=random_number(_x1,_x0);
	y=random_number(_x1,_x0);
	z=random_number(_x1,_x0);

	points_x[0]=x;
	points_y[0]=y;
	points_z[0]=z;

	//initiate number od hit darts
	int hit=0;

	while(true)
	{
		//generate another random seed
		x=random_number(_x1,_x0);
		y=random_number(_x1,_x0);
		z=random_number(_x1,_x0);
		double min_distance=1000;
		double distance=0;

		for (int j = 0; j < hit+1; j++)
		{
			distance=dist_two_points_3d(points_x[j],points_y[j],points_z[j],x,y,z);

			if (distance<min_distance) min_distance=distance;

		}
		
		// chck the dart hit or not 
		if (min_distance>r)
		{
			hit++;
			points_x[hit]=x;
			points_y[hit]=y;
			points_z[hit]=z;
		}

		num_itr++;
		//break if number of dart reaches the desire number 
		if(hit==num_points-1)	break;
		
		//break from this method if number of itrations exceeds 100000 and swithch to seedstracking_lastparent method 
		if(num_itr>100000)	break;
	}

	//seedstracking(hit,points_x,points_y, points_z,r);
	//seedstracking_bycub(hit,points_x,points_y, points_z,r);
	int check_ini;
	check_ini=seedstracking_lastparent(hit,points_x,points_y, points_z,r);

}
