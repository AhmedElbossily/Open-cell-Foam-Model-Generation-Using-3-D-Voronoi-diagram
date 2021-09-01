#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "tetgen.h"
#include "util.h"

using namespace std;

int main()
{
	// creat array that holds voroni vertices
	double *draw_points;
	draw_points = new double[60000000];
	int num_draw_points = 0;

	//Create a parameter to check whether the user want to generate a honeycomb structure 1 or random voronoi diagram
	int reg_type = 1;
	cout << "What strructure you want in the cube? \n";
	cout << "Enter (1) for  Honeycomb Structure \n";
	cout << "Enter (2) for  Irregular Voronoi Structure \n";
	cin >> reg_type;

	//the required number of seed inside the unite cube
	int num_points = 559;
	cout << "Enter the number of seeds you want inside the unite cube \n";
	cin >> num_points;

	int num_itr = 1;
	double *points_x, *points_y, *points_z;

	// Create variables that holds Vornoi ligaments characteristics
	double beam_rad = 0.0;
	int num_small_lig = 0;
	double denesity;

	// calculate the required distance to generate complete regulare voronoi structure for "num_points" seeds in unite cube "_x0*_x0*_x0"
	double r_reg;
	r_reg = min_dist_req_3d(_x0 * _x0 * _x0, num_points);

	//in case of user want to generate complete regulare structure "honeycomb"
	if (reg_type == 1)
	{
		// set the the required density of model
		denesity = 0.03;

		// Create arrays holds seeds coordinates
		points_x = new double[num_points];
		points_y = new double[num_points];
		points_z = new double[num_points];

		// generate regulare distribution seeds
		regular_points_seeds(points_x, points_y, points_z, r_reg);

		// Creating Vornoi vertices for honeycomb structure using TetGen library
		tet_gen_honeycomb(num_points, points_x, points_y, points_z, beam_rad, denesity, num_draw_points, draw_points);

		//set the destination for storing DXF file
		string str = "../output/00_Drawing.dxf";

		// Draw Honeycomb structure in dxf file
		draw_cells(draw_points, num_draw_points, str);

		//display voronoi ligaments radiaus in inch and in m
		cout << "beam_rad	:" << beam_rad << "	inch" << endl;
		cout << "beam_diameter	:" << beam_rad * .0254 * 2 << "	m" << endl;

		delete[] points_x, points_y, points_z;
	}
	//in case of user want to generete random distributing voronoi
	else
	{
		// random number geration function initiation
		int initiation = 0;

		//Regularity parameter 0 % complete irrigulare structure and 100% is honeycomb structure
		double e = 0.9;

		//these values are evaluated from running the code before. the code converge more quicly when srand() functions initiated with these values
		int good_initiation[20] = {2, 28, 29, 37, 45, 56, 70, 95, 104, 126, 162, 178, 205, 229, 243, 268, 309, 322, 326, 357};

		//create arrays holps seeds coordinates
		points_x = new double[num_points];
		points_y = new double[num_points];
		points_z = new double[num_points];

		//location where the Dxf file and statistics .txt will be stored
		string str = "../output/";

		// this variable will hold density itration
		int e_itra = 0;
		double m_rad_1 = 0.0, m_rad_3 = 0.0, m_rad_5 = 0.0, m_rad_7 = 0.0, m_rad_9 = 0.0;

		//we need to generate 20 model automatically and store them in different distinations
		while (initiation < 20)
		{
			e_itra = 0;
			num_itr = 0;

			//use the modified MPS to generate vornoi with regularity 90%
			mps_3d_for_high_regularity(points_x, points_y, points_z, num_points, e * r_reg, num_itr, good_initiation[initiation]);
			//mps_3d(points_x,points_y,points_z, num_points, e*r_reg,num_itr,initiation);

			while (e_itra < 5)
			{
				if (e_itra == 0)
					denesity = .01;
				if (e_itra == 1)
					denesity = .03;
				if (e_itra == 2)
					denesity = .05;
				if (e_itra == 3)
					denesity = .07;
				if (e_itra == 4)
					denesity = .09;

				// Get Vornoi vertices using TetGen library
				tet_gen(num_points, points_x, points_y, points_z, beam_rad, denesity, r_reg, e, initiation, e * r_reg, num_small_lig, num_draw_points, draw_points);

				// sum the beam diameter for every density and at the end calculate the mean for 20 model
				if (e_itra == 0)
					m_rad_1 = m_rad_1 + beam_rad;
				if (e_itra == 1)
					m_rad_3 = m_rad_3 + beam_rad;
				if (e_itra == 2)
					m_rad_5 = m_rad_5 + beam_rad;
				if (e_itra == 3)
					m_rad_7 = m_rad_7 + beam_rad;
				if (e_itra == 4)
					m_rad_9 = m_rad_9 + beam_rad;

				// Draw Voronoi Seeds for the first itration of densities as the structure will not change by changing densities
				// only Vornoi legaments diameter will change
				string path = str + to_string(initiation) + ".dxf";
				if (e_itra == 0)
					draw_cells(draw_points, num_draw_points, path);

				string path_2 = str + to_string(initiation) + ".txt";

				if (e_itra == 0)
				{
					ofstream beam_radius(path_2, ios::out);
					beam_radius << "denesity	" << denesity << endl;
					beam_radius << "beam_rad	" << beam_rad << "	inch" << endl;
					beam_radius << "num_small_lig	" << num_small_lig << endl;
					beam_radius << "regularity	" << e << endl;
					beam_radius << "***********************" << endl;
				}

				else
				{
					ofstream beam_radius(path_2, ios::app);
					beam_radius << "denesity	" << denesity << endl;
					beam_radius << "beam_rad	" << beam_rad << "	inch" << endl;
					beam_radius << "num_small_lig	" << num_small_lig << endl;
					beam_radius << "regularity	" << e << endl;
					beam_radius << "***********************" << endl;
				}

				e_itra++;
			}
			initiation++;
		}

		// calculate the mean of Voronoi ligaments for 20 model and store them in dxf file
		ofstream beam_radius("../output/mean_beam_radius.txt", ios::out);
		beam_radius << "regularity	" << e << endl;
		beam_radius << "m_rad_1	" << m_rad_1 / initiation << endl;
		beam_radius << "m_rad_3	" << m_rad_3 / initiation << endl;
		beam_radius << "m_rad_5	" << m_rad_5 / initiation << endl;
		beam_radius << "m_rad_7	" << m_rad_7 / initiation << endl;
		beam_radius << "m_rad_9	" << m_rad_9 / initiation << endl;

		delete[] points_x, points_y, points_z;
	}
	return 0;
}
