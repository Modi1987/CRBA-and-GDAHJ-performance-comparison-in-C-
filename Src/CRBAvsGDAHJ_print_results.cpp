//============================================================================
// Name        : CRBAvsGDAHJ_print_results.cpp
// Author      : Mohammad SAFEEA
// Version     : 
// Copyright   : Provided under MIT license
// Description : Used to compare execution time required for calculating mass matrix of serially linked robot
// 				 using CRBA and GDAHJ
// Both CRBA and GDAHJ implemented are the most efficient.
//-----------------------------------------------------------------------------
// The Composite Rigid Body Algorithm implemented here pertains to the most efficient algorithm
// As described by Featherstone in the reference:
// [a] Book: "rigid body dynamics algorithms", by Roy Featherstone
// All efficiency tricks described in the appindices of [a] are implemented, resulting
// in CRBA with minimal cost of:
// (10n2+22n-32) mul + (6n2+37n-43) add, where:
// n: is the number of degrees of freedom of serially linked robots
// n2: is square of n
// mul: stands for multiplication floating point operation
// add: stands for addition floating point operation
//-----------------------------------------------------------------------------
// GDAHJ implemented here is as described in the attached article,
// The presneted GDAHJ method acheives a minimal cost in the O(n2), 
// reducing the O(n2) operations required from 16 (for CRBA) to 5.5 (for GDAHJ). 
//============================================================================

#include <chrono>
#include <ctime>
#include <fstream>
#include <string>
#include <cstdlib>
#include <math.h>
#include <random>
#include <iostream>
#include <thread>
//Max number of degrees of freedom of robot
const int MAX_Dimention = 30;
const int MAX_DimentionPlusOne = MAX_Dimention + 1;
const int iterations = 500000; // number of calculations of the mass matrix

/****** IMPORTANT NOTE ******/
// The indexing of all matrices in this program starts from 1. As such, each element with index [0] was left unused.
// The reason, this is a re-implementation of the original code written in MATLAB (MATLAB indexing convention starts from 1)

/* Robot Kinematic PARAMETERS */
int NDF_robot = 2; // number of degrees of freedom of your serially linked robot
double dh_a[MAX_DimentionPlusOne]; // the elements {1 to NDF_robot} of this matrix contani the {a} parameter of modified-DH of your robot, in this code those elements are generated randomly
double dh_alfa[MAX_DimentionPlusOne]; // the elements {1 to NDF_robot} of this matrix contanin the {a} parameters of modified-DH of the robot. The rest of the elements are left un-used. In this code the elements of this matrix are generated randomly.
double dh_d[MAX_DimentionPlusOne]; // the elements {1 to NDF_robot} of this matrix contanin the {d} parameters of modified-DH of the robot. The rest of the elements are left un-used. In this code the elements of this matrix are generated randomly.
double q_JOINTS[MAX_DimentionPlusOne]; // the elements {1 to NDF_robot} of this matrix contanin the JOINT ANGLES of the robot. The rest of the elements are left un-used. In this code the elements of this matrix are generated randomly.

/* VARIABLES FOR  GDAHJ method, which is implemented in function {get_A_GDAHJ}, ARE DECLARED AS GLOBAL VARIABLES HERE */
double A_GDAHJ[MAX_DimentionPlusOne][MAX_DimentionPlusOne]; // Mass matrix calculated in GDAHJ method, this matrix is calculated using the function {get_A_GDAHJ}
double Tgdahj[5][5][MAX_DimentionPlusOne]; // transofrm matrices of the robot at some configuration, this tensor is assigned after calling the function {forward_Kinematics}
/* Inertial constants for GDAHJ method */
double mcii[MAX_DimentionPlusOne]; // constant matrix, each element [i] of this matrix contains the mass of link [i] (here generated randomly, in a real world case elements shall be assigned according to your robot)
double neg_mcii[MAX_DimentionPlusOne]; // constant matrix, each element [i] of this matrix contains the nigative of the mass of link [i]
double sigmamk[MAX_DimentionPlusOne]; // masses summed in a specific way, {if you are interested in more info you can refer to the code}.
double Pcii[4][MAX_DimentionPlusOne]; // constant matrix, elements [1 to 3][i] of this matrix represnet the position of center of mass of link [i] described in its local frame.
double Icii_mini_eig[MAX_DimentionPlusOne]; // constant matrix, each element [i] contains the minimum eigen value of inertial tensor of link [i].
double Icii_Qtilde[4][3][MAX_DimentionPlusOne]; // constant tensor, each matrix [1 to 3][1 to 2][i] contains a special factorization of the ienrtial tensor of link [i], for more info on way of factorization used please refer to submitted article.  
double traces[MAX_DimentionPlusOne]; // constant matrix, each element [i] of this matrix contains the trace of the inertial tensors of link [i]
double sum_of_traces[MAX_DimentionPlusOne]; // constant matrix, each element contains the previous traces summed in a particular way, {if you are interested in more info you can refer to the code}.

/* VARIABLES FOR  CRBA method, which is implemented in function {get_A_CRBA}, ARE DECLARED AS GLOBAL VARIABLES HERE */
double A_CRBA[MAX_DimentionPlusOne][MAX_DimentionPlusOne];
double crba_IC_cp[11][MAX_DimentionPlusOne]; // composite rigid body inertias, calculated when function {CRBA_minimalCost_11} is called
double crba_ct[MAX_DimentionPlusOne]; //  array of cosine of joint angles
double crba_st[MAX_DimentionPlusOne]; // array of sine of joint angles
/* Inertial constants for CRBA method */
double crba_IC[11][MAX_DimentionPlusOne]; // inertia constants of links with CRBA convention [I;Pcii;mcii]
double crba_ca[MAX_DimentionPlusOne]; // array of cosine of dh_alfa
double crba_sa[MAX_DimentionPlusOne]; // array of sine of dh_alfa
double crba_csa[MAX_DimentionPlusOne]; // array of cos(dh_alfa).*sin(dh_alfa) 
double crba_s2a[MAX_DimentionPlusOne]; // array of sin(dh_alfa).*sin(dh_alfa) 
double crba_p2csa[MAX_DimentionPlusOne]; // array of 2.*cos(dh_alfa).*sin(dh_lafa) 
double crba_p1_2s2a[MAX_DimentionPlusOne]; // array of 1-sin(dh_alfa).*sin(dh_lafa) 
double crba_sigmam[MAX_DimentionPlusOne]; // array of masses summed in a particular way, for more info check the code below.
double crba_am[MAX_DimentionPlusOne]; // array of sum of masses multiplied by parameters dh_a, for more info check the code below. 
// for more info about the physical meaning of the previous constants check the reference [a] mentioned above

/* Functions declaration */
double fRand(double fMin, double fMax); // used to generate random number in the range [fMin,fMax]
void print_CRBA_time(double time, int NDF); // used to print timing data for CRBA method
void print_GDAHJ_time(double time, int NDF); // used to print timing data for GDAHJ method
void get_A_GDAHJ(int NDOF); // used for calculating mass matrix using GDAHJ
void forward_Kinematics(int DOFplusOne); // used for calculating transform amtrices {forward kinematics} of the links
void calculateInertialTensors(); // used for calculating generating inertia tensors for links of the robot, and also for storing them in CRBA convention and GDAHJ convention
void generate_A_RandomOrthonormal(); // used to generate a random orthonormal matrix, this orthornomral matrix is stored in matrix {Rxyz declared below}, and is used for calculating the a random inertia tensor
void print_A_GDAHJ(int NDF); // used to print the calculation time for GDAHJ method
void disp(double x); // to display a double on console
void disp(const char* input); // to disply a string on console
void initiate_GDAHJ_Constants(int dof); // used to calculate the offline constants for GDAHJ
void initiate_CRBA_Constants(int dof); // used to calculate the offline constants for CRBA like {cos(alfa),sin(alfa), cos(alfa)sin(alfa), ....etc} for a total list of the constants refer to the body of the function
void CRBA_minimalCost_11(int crba_ndf); // Most efficient CRBA implementation for calculating the joint space inertia matrix
void print_A_CRBA(int NDF); // prints the time required by CRBA to calculate the mass matrix
double get_time_milli(); // used to get the current time in millisconds

/* Auxuliary variables*/
double Rxyz[4][4];// used by function (returnRandomRotation) to store value of randomly generated orthonormal matrix that is used to calculate a random inertia tensor

/* The main function of the C++ program, comparison is done in here*/
int main()
{
	disp("Comparing CRBA vs GDAHJ for calculating JSIM !");
	disp("Inertial data of robot links are generated randomly");
	disp("Results are formatted as:"); 
	disp("[ DOF, time (GDAHJ), time (CRBA), error as: sum(sum(abs(A_CRBA-AGDAHJ))) /n^2]");
	int DOFplusOne = NDF_robot + 1;
	double fMin = 0.0;
	double fMax = 1.0;
	// Generate kinematic parameters (modified DH) randomly, insert your robot's parameters for real-world use.
	for (int i = 1; i < MAX_DimentionPlusOne; i++)
	{
		dh_alfa[i] = fRand(fMin, fMax);
		dh_a[i] = fRand(fMin, fMax);
		dh_d[i] = fRand(fMin, fMax);
	}
	// generate the configuration of the robot, joint angles, randomly. Insert your robot configuration for real-world use.
	for (int i = 1; i < MAX_DimentionPlusOne; i++)
	{
		// joint angle in range [-pi,pi]
		q_JOINTS[i] = fRand(-4 * atan(1), 4 * atan(1));
	}
	// Calculate the forward Kinematics, transform matrices (Tgdahj) of robot.
	forward_Kinematics(MAX_DimentionPlusOne);

	// Generate random Inertia data, in real-wrold scinario, you shall insert your robot inertial data
	for (int i = 1; i < MAX_DimentionPlusOne; i++)
	{
		for (int j = 1; j < 4; j++)
		{
			Pcii[j][i] = fRand(fMin, fMax);
		}
		mcii[i] = fRand(fMin, fMax);
		neg_mcii[i] = -mcii[i];
		// Initiate global variables with zeros before calculating them later
		sigmamk[i] = 0.0;
		sum_of_traces[i] = 0.0;
		traces[i] = 0.0;
	}

	// Calculate, Icii_Qtilde, Icii_mini_eig, traces, IciiMinimal (the data structure of the inertial tensors implemented by GDAHJ)
	calculateInertialTensors();

	// Perform the comparision between CRBA and GDAHJ for serial robots with different DOF
	std::cout << "res=[";
	while (true)
	{
		std::cout << NDF_robot << ","; // displays DOF of the robot
		//********** GDAHJ ***************//
		initiate_GDAHJ_Constants(NDF_robot); // initiate the constants of the GDAHJ method
		double t_GDAHJ_0 = get_time_milli();
		// Here GDAHJ is invoked 500 000  times
		for (int calc_iteration = 0; calc_iteration < iterations; calc_iteration++)
		{
			get_A_GDAHJ(NDF_robot);
		}
		double t_GDAHJ_1 = get_time_milli();
		double dt_GDAHJ = t_GDAHJ_1 - t_GDAHJ_0;
		std::cout << dt_GDAHJ << ",";


		//********** CRBA ***************//
		initiate_CRBA_Constants(NDF_robot);
		double t_CRBA_0 = get_time_milli();
		// Here CRBA is invoked 500 000  times
		for (int calc_iteration = 0; calc_iteration < iterations; calc_iteration++)
		{
			CRBA_minimalCost_11(NDF_robot);
		}
		double t_CRBA_1 = get_time_milli();
		double dt_CRBA = t_CRBA_1 - t_CRBA_0;
		std::cout << dt_CRBA << ",";


		//********** Calculating the error (sum of absolute value of differences) ***************//
		double dasum = 0.0;
		for (int i = 1; i < NDF_robot + 1; i++)
		{
			for (int j = 1; j < NDF_robot + 1; j++)
			{
				dasum = dasum + abs(A_CRBA[i][j] - A_GDAHJ[i][j]);
			}
		}
		double avg = dasum / ((double)NDF_robot * NDF_robot);
		std::cout << avg;
		NDF_robot = NDF_robot +1;
		if (NDF_robot == MAX_DimentionPlusOne)
		{
			std::cout << "]" << "\n";
			break;
		}
		else
		{
			std::cout<< ";\n";
		}		
	}

	std::this_thread::sleep_for(std::chrono::nanoseconds(10));
	print_A_GDAHJ(MAX_Dimention);
	print_A_CRBA(MAX_Dimention);
	return 0;
}
/***************************************************/
void initiate_GDAHJ_Constants(int dof)
{
	// Calculate off-line constants for GDAHJ (here indeces are for robots with NDF degrees of freedom)
	for (int i = dof; i > 0; i--)
	{
		if (i == dof)
		{
			sigmamk[i] = 0.0;
			sum_of_traces[i] = traces[i];
		}
		else
		{
			sigmamk[i] = sigmamk[i + 1] + mcii[i + 1];
			sum_of_traces[i] = sum_of_traces[i + 1] + traces[i];
			//disp(sum_of_traces[i]);
		}
	}
}
void initiate_CRBA_Constants(int dof)
{
	// int dof_plus_one = dof + 1;
	for (int i = 1; i < MAX_DimentionPlusOne; i++)
	{
		crba_ct[i] = cos(q_JOINTS[i]);
		crba_st[i] = sin(q_JOINTS[i]);
		crba_ca[i] = cos(dh_alfa[i]);
		crba_sa[i] = sin(dh_alfa[i]);
		crba_csa[i] = crba_ca[i] * crba_sa[i];
		crba_s2a[i] = crba_sa[i] * crba_sa[i];
		crba_p2csa[i] = 2 * crba_csa[i];
		crba_p1_2s2a[i] = 1 - crba_s2a[i] - crba_s2a[i];
		/* Very important {crba_IC} is  initiated inside the function
		{calculateInertialTensors} */
	}
	crba_sigmam[dof] = mcii[dof];
	for (int i = dof - 1; i > 0; i--)
	{
		crba_sigmam[i] = crba_sigmam[i + 1] + mcii[i];
	}
	crba_am[dof] = 0;
	for (int i = dof; i > 0; i--)
	{
		crba_am[i] = crba_sigmam[i] * dh_a[i];
	}
}

void CRBA_minimalCost_11(int crba_ndf)
{
	double EZET[7];
	for (int i = 0; i < 7; i++) { EZET[i] = 0.0; }
	double Ey[4];
	for (int i = 0; i < 4; i++) { Ey[i] = 0.0; }
	// initiate the array of the composite rigid body inertias:
	for (int i = 1; i < 11; i++)
	{
		for (int j = 1; j < crba_ndf + 1; j++)
		{
			crba_IC_cp[i][j] = crba_IC[i][j];
		}
	}
	//// Composite inertia calculation, O[n] part of the algorithm
	for (int i = crba_ndf; i > 1; i--)
	{
		//// Calculating the articulated rigid bodies inertias, tow 
		// scre transforms are required:
		// The cost of each is 16 mul +15 add
		// Discribtion is available page 249, feathersotnes book,
		// Appindex A.4

		//// First transform:
		// In this case the transofr is: 
		// Xrotz[q[i]]*Xtrans[[0,0,d[i]]]'
		//------------------------------------------
		// vector h
		double h[4];
		double y[4];
		for (int k = 1; k < 4; k++)
		{
			h[k] = crba_IC_cp[7 + k - 1][i];
			y[k] = h[k];
		}
		// calculate vetor y,
		y[3] = y[3] + crba_IC_cp[10][i] * dh_d[i]; // description is at end of page 248, sign flipped 
		// due to Modefied DH convention
		// calculate Matrix Z: 
		// description is at end of page 248, sign flipped due to Modefied DH convention
		double z[7]; for (int k = 1; k < 7; k++) { z[k] = crba_IC_cp[k][i]; };
		double temp = dh_d[i] * (h[3] + y[3]);
		z[1] = z[1] + temp;
		z[2] = z[2] + temp;
		z[5] = z[5] - dh_d[i] * h[1];
		z[6] = z[6] - dh_d[i] * h[2];
		// rotate I_hat matrix around Z axis 
		// equations are described at page 255 of the book, 
		// section title "Rotation of a Symmetric Matrix about a Coordinate Axis"
		double cs = crba_ct[i] * crba_st[i];
		double s2 = crba_st[i] * crba_st[i];
		double p2cs = cs + cs;
		double p1_2s2 = 1 - s2 - s2;
		temp = z[2] - z[1];
		double alfa_z = s2 * temp - p2cs * z[4];
		double beta_z = p1_2s2 * z[4] - cs * temp;
		EZET[1] = z[1] + alfa_z;
		EZET[2] = z[2] - alfa_z;
		EZET[3] = z[3];
		EZET[4] = beta_z;
		EZET[5] = crba_ct[i] * z[5] - crba_st[i] * z[6];
		EZET[6] = crba_ct[i] * z[6] + crba_st[i] * z[5];
		// 2nd rotate y vector around Z axis with angle q[i]
		Ey[1] = crba_ct[i] * y[1] - crba_st[i] * y[2];
		Ey[2] = crba_st[i] * y[1] + crba_ct[i] * y[2];
		Ey[3] = y[3];

		//// Second transform:
		// Xrotx[alfa[i]] * Xtrans[[a[i],0,0]]'
		//------------------------------------------
		// vector h
		for (int k = 1; k < 4; k++) { h[k] = Ey[k]; y[k] = h[k]; }
		// calculate vetor y
		y[1] = y[1] + crba_am[i];
		//     y[1]=y[1]+Ivec[10]*a[i];
			// calculate Matrix Z:
		for (int k = 1; k < 7; k++) { z[k] = EZET[k]; };
		temp = dh_a[i] * (h[1] + y[1]);
		z[2] = z[2] + temp;
		z[3] = z[3] + temp;
		z[4] = z[4] - dh_a[i] * h[2];
		z[5] = z[5] - dh_a[i] * h[3];
		// rotate I_hat amtrix around X axis
		// 1st rotate Z matrix
		temp = z[3] - z[2];
		double alfa_x = crba_s2a[i] * temp - crba_p2csa[i] * z[6];
		double beta_x = crba_p1_2s2a[i] * z[6] - crba_csa[i] * temp;
		EZET[1] = z[1];
		EZET[2] = z[2] + alfa_x;
		EZET[3] = z[3] - alfa_x;
		EZET[4] = crba_ca[i] * z[4] - crba_sa[i] * z[5];
		EZET[5] = crba_ca[i] * z[5] + crba_sa[i] * z[4];
		EZET[6] = beta_x;
		// 2nd rotate y vector around X axis with angle alfa[i]
		Ey[1] = y[1];
		Ey[2] = crba_ca[i] * y[2] - crba_sa[i] * y[3];
		Ey[3] = crba_sa[i] * y[2] + crba_ca[i] * y[3];

		//// Calculate the composite rigid bodies inertia
		crba_IC_cp[1][i - 1] = crba_IC_cp[1][i - 1] + EZET[1];
		crba_IC_cp[2][i - 1] = crba_IC_cp[2][i - 1] + EZET[2];
		crba_IC_cp[3][i - 1] = crba_IC_cp[3][i - 1] + EZET[3];
		crba_IC_cp[4][i - 1] = crba_IC_cp[4][i - 1] + EZET[4];
		crba_IC_cp[5][i - 1] = crba_IC_cp[5][i - 1] + EZET[5];
		crba_IC_cp[6][i - 1] = crba_IC_cp[6][i - 1] + EZET[6];
		crba_IC_cp[7][i - 1] = crba_IC_cp[7][i - 1] + Ey[1];
		crba_IC_cp[8][i - 1] = crba_IC_cp[8][i - 1] + Ey[2];
		crba_IC_cp[9][i - 1] = crba_IC_cp[9][i - 1] + Ey[3];
		crba_IC_cp[10][i - 1] = crba_sigmam[i - 1];
		//// Previous optimized code is equivelent to the following 
		// two compact statements:   
		//     Tf{i}=Xrotx[alfa[i]]*Xtrans[[a[i],0,0]]'*Xrotz[q[i]]*Xtrans[[0,0,d[i]]]';
		//     IC{i-1} = IC{i-1} + Tf{i}*IC{i}*Tf{i}';
	}

	for (int i = 1; i < crba_ndf + 1; i++)
	{
		A_CRBA[i][i] = crba_IC_cp[3][i];
		int j = i;

		double mx = crba_IC_cp[5][i];
		double my = crba_IC_cp[6][i];
		double mz = crba_IC_cp[3][i];
		double fx = -crba_IC_cp[8][i];
		double fy = crba_IC_cp[7][i];
		double fz = 0;
		for (int kkk = (i - 1); kkk > 0; kkk--)
		{
			double temp1 = mx - fy * dh_d[j];
			double temp2 = my + fx * dh_d[j];
			double mx1 = crba_ct[j] * temp1 - crba_st[j] * temp2;
			double my1 = crba_ct[j] * temp2 + crba_st[j] * temp1;
			double mz1 = mz;
			double fx1 = fx * crba_ct[j] - fy * crba_st[j];
			double fy1 = fx * crba_st[j] + fy * crba_ct[j];
			double fz1 = fz;
			mx = mx1;
			my = my1;
			mz = mz1;
			fx = fx1;
			fy = fy1;
			fz = fz1;
			mx1 = mx;
			temp1 = mz + dh_a[j] * fy;
			temp2 = my - dh_a[j] * fz;
			my1 = -crba_sa[j] * temp1 + crba_ca[j] * temp2;
			mz1 = crba_sa[j] * temp2 + crba_ca[j] * temp1;
			fx1 = fx;
			fy1 = crba_ca[j] * fy - crba_sa[j] * fz;
			fz1 = fy * crba_sa[j] + fz * crba_ca[j];

			j = kkk;
			A_CRBA[i][j] = mz1;
			A_CRBA[j][i] = mz1;
			mx = mx1;
			my = my1;
			mz = mz1;
			fx = fx1;
			fy = fy1;
			fz = fz1;
		}
	}


	//// About
	// A C++ implementation of:
	// Featherstone's (CRBA) with all effeciency tricks as described in the book  
	// "Rigid Body Dynamics Algorithms", Springer, 2008
	// The implemented code acheives a minimal cost of:
	// (10n2+22n-32) mul + (6n2+37n-43) add, where:
	// n: is the number of degrees of freedom of serially linked robots
	// n2: is square of n
	// mul: stands for multiplication floating point operation
	// add: stands for addition floating point operation

	// Written by: Mohammad SAFEEA, 17th May 2019


}

/***************************************************/
// Print formatters for the mass matrices
int mask_elements_from_index = 4;
/***************************************************/
void print_A_GDAHJ(int NDF)
{
	int unmask_lements_from_index = NDF - 4;
	int ndf_plus_one = NDF + 1;
	std::cout << "Inertia matrix calculated using GDAHJ for "<< NDF << "DOF serially linked robot is: \n";
	bool print_flag;
	for (int i = 1; i < ndf_plus_one; i++)
	{
		print_flag = true;
		for (int j = 1; j < ndf_plus_one; j++)
		{
			if (j > mask_elements_from_index && j < unmask_lements_from_index)
			{
				if (print_flag == true)
				{
					std::cout << "[ . . . . . . . . . ]";
					print_flag = false;
				}
			}
			else
			{
				std::cout << "[" << A_GDAHJ[i][j] << "]";
			}
		}
		std::cout << "\n";
	}
}
/***************************************************/
void print_A_CRBA(int NDF)
{
	int unmask_lements_from_index = NDF - 4;
	int ndf_plus_one = NDF + 1;
	std::cout << "Inertia matrix calculated using CRBA for " << NDF << "DOF serially linked robot is: \n";
	bool print_flag;
	print_flag = true;
	for (int i = 1; i < ndf_plus_one; i++)
	{
		print_flag = true;
		for (int j = 1; j < ndf_plus_one; j++)
		{
			if (j > mask_elements_from_index && j < unmask_lements_from_index)
			{
				if (print_flag == true)
				{
					std::cout << "[ . . . . . . . . . ]";
					print_flag = false;
				}
			}
			else
			{
				std::cout << "[" << A_CRBA[i][j] << "]";
			}
		}
		std::cout << "\n";
	}
}
/***************************************************/
void calculateInertialTensors()
{
	double a, b, g;
	// A function for calculating, Icii_Qtilde, Icii_mini_eig, traces, IciiMinimal
	for (int count = 1; count < MAX_DimentionPlusOne; count++)
	{
		// Generating a random rotation matrix (and storing it in matrix Rxyz)
		generate_A_RandomOrthonormal();
		double raR[4][4]; // A rotation matrix to be generated randomly, for use in random generation of an Intertial tensor
		for (int i = 1; i < 4; i++)
		{
			for (int j = 1; j < 4; j++)
			{
				raR[i][j] = Rxyz[i][j];
			}
		}
		// generating random eigen values for inertial tensor
		double fMin = 0.0;
		double fMax = 1.0;
		double eig[4];
		eig[1] = fRand(fMin, fMax);
		eig[2] = fRand(fMin, fMax);
		eig[3] = fRand(fMin, fMax);
		// find minimum eigen value
		double min_eig = eig[1];
		int index = 1;
		if (min_eig > eig[2]) { min_eig = eig[2]; index = 2; }
		if (min_eig > eig[3]) { min_eig = eig[3]; index = 3; }
		// calculate constant tensors for GDAHJ:
		traces[count] = eig[1] + eig[2] + eig[3];
		Icii_mini_eig[count] = min_eig;
		int tempIndex = 0;
		for (int kk = 1; kk < 4; kk++)
		{
			if (kk == index)
			{
			}
			else
			{
				double temp = std::sqrt(eig[kk] - min_eig);
				tempIndex = tempIndex + 1;
				for (int jjk = 1; jjk < 4; jjk++)
				{
					Icii_Qtilde[jjk][tempIndex][count] = temp * raR[jjk][kk];
					// disp(Icii_Qtilde[jjk][tempIndex][i]);
				}

			}
		}
		double eigMat_Rt[4][4]; // matrix containing [diag(eig)*R']
		double Icii_link[4][4]; // inertia tensor for each link
		for (int i = 1; i < 4; i++)
		{
			for (int j = 1; j < 4; j++)
			{
				eigMat_Rt[i][j] = eig[i] * raR[j][i];
			}
		}
		for (int i = 1; i < 4; i++)
		{
			for (int j = 1; j < 4; j++)
			{
				Icii_link[i][j] = raR[i][1] * eigMat_Rt[1][j] +
					raR[i][2] * eigMat_Rt[2][j] +
					raR[i][3] * eigMat_Rt[3][j];
			}
		}

		// calculating {Icii(:,:,i)+mcii(i)*skew(Pcii(:,i))*skew(Pcii(:,i))';}
		crba_IC[1][count] = Icii_link[1][1] + mcii[count] * (Pcii[2][count] * Pcii[2][count]
			+ Pcii[3][count] * Pcii[3][count]);
		crba_IC[2][count] = Icii_link[2][2] + mcii[count] * (Pcii[1][count] * Pcii[1][count]
			+ Pcii[3][count] * Pcii[3][count]);
		crba_IC[3][count] = Icii_link[3][3] + mcii[count] * (Pcii[1][count] * Pcii[1][count] +
			Pcii[2][count] * Pcii[2][count]);
		crba_IC[4][count] = Icii_link[1][2] + mcii[count] * (-Pcii[1][count] * Pcii[2][count]);
		crba_IC[5][count] = Icii_link[1][3] + mcii[count] * (-Pcii[1][count] * Pcii[3][count]);
		crba_IC[6][count] = Icii_link[2][3] + mcii[count] * (-Pcii[2][count] * Pcii[3][count]);
		crba_IC[7][count] = mcii[count] * Pcii[1][count];
		crba_IC[8][count] = mcii[count] * Pcii[2][count];
		crba_IC[9][count] = mcii[count] * Pcii[3][count];
		crba_IC[10][count] = mcii[count];
		/*
		// Following code was used for debug, to check that the link inertia (CRBA minimal convenction) are calculated properly
		// using the debug, link inertias are being calculated properly
		std::cout << "The inertia vector of Link " << count << " is: ";
		for (int elementCount = 1; elementCount < 11; elementCount++)
		{
			std::cout << crba_IC[elementCount][count] << "  ";
		}
		std::cout << "\n";
		*/
	}
}
// disp("previous numbsers are Icii_Qtilde")

void generate_A_RandomOrthonormal()
{
	double fMin = 0.0;
	double fMax = 1.0;
	// Generating a random rotation angles
	double a = fRand(fMin, fMax);
	double b = fRand(fMin, fMax);
	double g = fRand(fMin, fMax);

	double Rx[4][4];
	double Ry[4][4];
	double Rxy[4][4];
	double Rz[4][4];
	// Rxyz is defined as a global variable

	for (int i = 1; i < 4; i++)
	{
		for (int j = 1; j < 4; j++)
		{
			Rx[i][j] = 0.0;
			Ry[i][j] = 0.0;
			Rz[i][j] = 0.0;
			Rxy[i][j] = 0.0;
			Rxyz[i][j] = 0.0;
		}
	}
	// Calculate Rx
	Rx[1][1] = cos(a); Rx[1][2] = -sin(a);
	Rx[2][1] = sin(a); Rx[2][2] = cos(a);
	Rx[3][3] = 1;
	// Calculate Ry
	Ry[1][1] = cos(b); Ry[1][3] = sin(b);
	Ry[2][2] = 1;
	Ry[3][1] = -sin(b); Ry[3][3] = cos(b);
	// Calculate Rz
	Rz[1][1] = 1;
	Rz[2][2] = cos(g); Rz[2][3] = -sin(g);
	Rz[3][2] = sin(g); Rz[3][3] = cos(g);

	for (int i = 1; i < 4; i++)
	{
		for (int j = 1; j < 4; j++)
		{
			for (int k = 1; k < 4; k++)
			{
				Rxy[i][k] = Rxy[i][k] + Rx[i][j] * Ry[j][k];
			}

		}
	}

	for (int i = 1; i < 4; i++)
	{
		for (int j = 1; j < 4; j++)
		{
			for (int k = 1; k < 4; k++)
			{
				Rxyz[i][k] = Rxyz[i][k] + Rxy[i][j] * Rz[j][k];
			}

		}
	}
}


double get_time_milli()
{
	double sysTime = std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::system_clock::now().time_since_epoch()).count();
	return sysTime;
}

int separator_line_width = 35;
void print_GDAHJ_time(double time, int NDF)
{
	for (int i = 0; i < separator_line_width; i++)
	{
		std::cout << "-";
	}
	std::cout << "\n";
	std::cout << "GDAHJ-computation time (milliseconds) for " << NDF << " DOF is:" << std::endl;
	std::cout << time << std::endl;
	for (int i = 0; i < separator_line_width; i++)
	{
		std::cout << "-";
	}
	std::cout << "\n";
}

void print_CRBA_time(double time, int NDF)
{
	for (int i = 0; i < separator_line_width; i++)
	{
		std::cout << "-";
	}
	std::cout << "\n";
	std::cout << "CRBA-computation time (milliseconds) for " << NDF << " DOF is:" << std::endl;
	std::cout << time << std::endl;
	for (int i = 0; i < separator_line_width; i++)
	{
		std::cout << "-";
	}
	std::cout << "\n";
}


double fRand(double fMin, double fMax)
{
	srand(time(NULL));
	double f = ((double)rand()) / RAND_MAX;
	return fMin + f * (fMax - fMin);
	//return 2.12358;
}





void forward_Kinematics(int DOFplusOne)
{
	double T1[5][5];

	for (int i = 1; i < 5; i++)
	{
		for (int j = 1; j < 5; j++)
		{
			if (i == j)
			{
				Tgdahj[i][j][0] = 1.0;
			}
			else
			{
				Tgdahj[i][j][0] = 0.0;
			}
		}
	}
	for (int count = 1; count < DOFplusOne; count++)
	{
		double alfaVar = dh_alfa[count];
		double qVar = q_JOINTS[count];
		double aVar = dh_a[count];
		double dVar = dh_d[count];
		T1[1][1] = cos(qVar); T1[1][2] = -sin(qVar); T1[1][3] = 0.0; ; T1[1][4] = aVar;
		T1[2][1] = sin(qVar) * cos(alfaVar); T1[2][2] = cos(qVar) * cos(alfaVar); T1[2][3] = -sin(alfaVar); T1[2][4] = -sin(alfaVar) * dVar;
		T1[3][1] = sin(qVar) * sin(alfaVar); T1[3][2] = cos(qVar) * sin(alfaVar); T1[3][3] = 0.0; T1[3][3] = cos(alfaVar); T1[3][4] = cos(alfaVar) * dVar;
		T1[4][1] = 0.0; T1[4][2] = 0.0; T1[4][3] = 0.0; T1[4][4] = 1;
		// display T1 (DH matrix)
		/**
		// std::cout << "DH matrix for link " << count << " is:\n";
		for (int i = 1; i < 5; i++)
		{
			for (int j = 1; j < 5; j++)
			{
				std::cout << T1[i][j] << "  ";
			}
			std::cout << "\n";
		}
		// previous double loop is to show the DH matrix of each link
		**/

		// disp(T1[1][1]); disp(" previous number is T1[1][1]");
		for (int i = 1; i < 5; i++)
		{
			for (int j = 1; j < 5; j++)
			{
				Tgdahj[i][j][count] = Tgdahj[i][1][count - 1] * T1[1][j] +
					Tgdahj[i][2][count - 1] * T1[2][j] +
					Tgdahj[i][3][count - 1] * T1[3][j] +
					Tgdahj[i][4][count - 1] * T1[4][j];
			}
		}
		// disp(Tgdahj[1][1][count]);
		// disp(" previous numsbers are Tgfahj[1][1]");
	}


}

void get_A_GDAHJ(int NDOF)
{
	/*
	 * GDAHJ code implementation
	 *
	 *  Created on: 17/May/2019
	 *      Author: Mohammad SAFEEA
	 *
	 *
	 * Following code defines a function that calculates the mass matrix of
	 * the manipulator using GDAHJ algorithm.
	 *
	 */

	 /* List of arguments, declared as global variables:
	 double A_GDAHJ[][DOFplusOne],
	 double Tgdahj[][5][DOFplusOne],
	 double Pcii[][DOFplusOne],
	 double Icii_mini_eig[]
	 double Icii_Qtilde[][3][DOFplusOne],
	 double mcii[],
	 double neg_mcii[]
	 double sigmamk[],
	 double sum_of_traces[]
	 */
	int DOFplusOne = NDOF + 1;
	// Declare some variables
	double li[4][MAX_DimentionPlusOne];
	double pcii[4][MAX_DimentionPlusOne];
	double pci[4][MAX_DimentionPlusOne];
	double sigmamkpck[4][MAX_DimentionPlusOne];
	double Di[4][MAX_DimentionPlusOne];
	double kjpj[4][MAX_DimentionPlusOne];
	double Ci[4][4][MAX_DimentionPlusOne];
	double Ei[4][4][MAX_DimentionPlusOne + 1];
	// Calculate (s) vectors Qtilde'*ei
	double s[3][4];
	double Irot[4][4];
	// initiate accumulators with zeros
	for (int count = 0; count < MAX_DimentionPlusOne; count++)
	{
		for (int i = 0; i < 4; i++)
		{
			li[i][count] = 0.0;
			pcii[i][count] = 0.0;
			pci[i][count] = 0.0;
			sigmamkpck[i][count] = 0.0;
			Di[i][count] = 0.0;
			kjpj[i][count] = 0.0;
			for (int j = 0; j < 4; j++)
			{
				Ci[i][j][count] = 0.0;
				Ei[i][j][count] = 0.0;
			}
		}
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			Ei[i][j][MAX_DimentionPlusOne] = 0.0;
		}
	}

	// Rotate inertia tensors
	for (int kkk = NDOF; kkk > 0; kkk--)
	{
		// In following the cost is 12 add + 18 mul
		for (int i = 1; i < 3; i++)
		{
			for (int j = 1; j < 4; j++)
			{
				s[i][j] = Tgdahj[j][1][kkk] * Icii_Qtilde[1][i][kkk]
					+ Tgdahj[j][2][kkk] * Icii_Qtilde[2][i][kkk] +
					Tgdahj[j][3][kkk] * Icii_Qtilde[3][i][kkk];
				// disp(Tgdahj[j][3][kkk]);
			}
		}
		for (int i = 1; i < 3; i++)
		{
			for (int j = i; j < 4; j++)
			{
				Irot[i][j] = s[1][i] * s[1][j] + s[2][i] * s[2][j];
			}
		}
		// In following the cost is 2 add
		Irot[1][1] = Icii_mini_eig[kkk] + Irot[1][1];
		Irot[2][2] = Icii_mini_eig[kkk] + Irot[2][2];
		// In total cost is 19 add + 28 mul
		Ei[1][1][kkk] = Ei[1][1][kkk + 1] + Irot[1][1];
		Ei[1][2][kkk] = Ei[1][2][kkk + 1] + Irot[1][2];
		Ei[1][3][kkk] = Ei[1][3][kkk + 1] + Irot[1][3];
		Ei[2][2][kkk] = Ei[2][2][kkk + 1] + Irot[2][2];
		Ei[2][3][kkk] = Ei[2][3][kkk + 1] + Irot[2][3];
		Ei[3][3][kkk] = sum_of_traces[kkk] - Ei[1][1][kkk] - Ei[2][2][kkk];
		Ei[2][1][kkk] = Ei[1][2][kkk];
		Ei[3][1][kkk] = Ei[1][3][kkk];
		Ei[3][2][kkk] = Ei[2][3][kkk];
	}

	//// variables initiation
	for (int i = NDOF; i > 0; i--)
	{

		pcii[1][i] = Tgdahj[1][1][i] * Pcii[1][i] + Tgdahj[1][2][i] * Pcii[2][i] + Tgdahj[1][3][i] * Pcii[3][i];
		pcii[2][i] = Tgdahj[2][1][i] * Pcii[1][i] + Tgdahj[2][2][i] * Pcii[2][i] + Tgdahj[2][3][i] * Pcii[3][i];
		pcii[3][i] = Tgdahj[3][1][i] * Pcii[1][i] + Tgdahj[3][2][i] * Pcii[2][i] + Tgdahj[3][3][i] * Pcii[3][i];

		pci[1][i] = Tgdahj[1][4][i] + pcii[1][i];
		pci[2][i] = Tgdahj[2][4][i] + pcii[2][i];
		pci[3][i] = Tgdahj[3][4][i] + pcii[3][i];

		kjpj[1][i] = Tgdahj[2][3][i] * Tgdahj[3][4][i] - Tgdahj[3][3][i] * Tgdahj[2][4][i];
		kjpj[2][i] = Tgdahj[3][3][i] * Tgdahj[1][4][i] - Tgdahj[1][3][i] * Tgdahj[3][4][i];
		kjpj[3][i] = Tgdahj[1][3][i] * Tgdahj[2][4][i] - Tgdahj[2][3][i] * Tgdahj[1][4][i];

	}
	//// Other variables initiation
	double mciipcii[4];
	for (int i = NDOF - 1; i > 0; i--)
	{
		li[1][i] = Tgdahj[1][4][i + 1] - Tgdahj[1][4][i];
		li[2][i] = Tgdahj[2][4][i + 1] - Tgdahj[2][4][i];
		li[3][i] = Tgdahj[3][4][i + 1] - Tgdahj[3][4][i];

		sigmamkpck[1][i] = sigmamkpck[1][i + 1] + mcii[i + 1] * pci[1][i + 1];
		sigmamkpck[2][i] = sigmamkpck[2][i + 1] + mcii[i + 1] * pci[2][i + 1];
		sigmamkpck[3][i] = sigmamkpck[3][i + 1] + mcii[i + 1] * pci[3][i + 1];
	}

	mciipcii[1] = neg_mcii[NDOF] * pcii[1][NDOF];
	mciipcii[2] = neg_mcii[NDOF] * pcii[2][NDOF];
	mciipcii[3] = neg_mcii[NDOF] * pcii[3][NDOF];
	double d11, d22, d33, b11, b22, b33;
	d11 = mciipcii[1] * pci[1][NDOF];
	d22 = mciipcii[2] * pci[2][NDOF];
	d33 = mciipcii[3] * pci[3][NDOF];

	b11 = li[1][NDOF] * sigmamkpck[1][NDOF];
	b22 = li[2][NDOF] * sigmamkpck[2][NDOF];
	b33 = li[3][NDOF] * sigmamkpck[3][NDOF];

	Ci[1][1][NDOF] = b33 + b22 - d33 - d22;
	Ci[2][1][NDOF] = mciipcii[1] * pci[2][NDOF] - li[1][NDOF] * sigmamkpck[2][NDOF];
	Ci[3][1][NDOF] = mciipcii[1] * pci[3][NDOF] - li[1][NDOF] * sigmamkpck[3][NDOF];
	Ci[1][2][NDOF] = mciipcii[2] * pci[1][NDOF] - li[2][NDOF] * sigmamkpck[1][NDOF];
	Ci[2][2][NDOF] = b33 + b11 - d33 - d11;
	Ci[3][2][NDOF] = mciipcii[2] * pci[3][NDOF] - li[2][NDOF] * sigmamkpck[3][NDOF];
	Ci[1][3][NDOF] = mciipcii[3] * pci[1][NDOF] - li[3][NDOF] * sigmamkpck[1][NDOF];
	Ci[2][3][NDOF] = mciipcii[3] * pci[2][NDOF] - li[3][NDOF] * sigmamkpck[2][NDOF];
	Ci[3][3][NDOF] = b22 + b11 - d22 - d11;

	Di[1][NDOF] = mciipcii[1] - sigmamk[NDOF] * li[1][NDOF];
	Di[2][NDOF] = mciipcii[2] - sigmamk[NDOF] * li[2][NDOF];
	Di[3][NDOF] = mciipcii[3] - sigmamk[NDOF] * li[3][NDOF];

	for (int i = NDOF - 1; i > 0; i--)
	{
		mciipcii[1] = neg_mcii[i] * pcii[1][i];
		mciipcii[2] = neg_mcii[i] * pcii[2][i];
		mciipcii[3] = neg_mcii[i] * pcii[3][i];

		d11 = mciipcii[1] * pci[1][i];
		d22 = mciipcii[2] * pci[2][i];
		d33 = mciipcii[3] * pci[3][i];

		b11 = li[1][i] * sigmamkpck[1][i];
		b22 = li[2][i] * sigmamkpck[2][i];
		b33 = li[3][i] * sigmamkpck[3][i];

		Ci[1][1][i] = Ci[1][1][i + 1] - d33 - d22 + b33 + b22;
		Ci[2][1][i] = Ci[2][1][i + 1] + mciipcii[1] * pci[2][i] - li[1][i] * sigmamkpck[2][i];
		Ci[3][1][i] = Ci[3][1][i + 1] + mciipcii[1] * pci[3][i] - li[1][i] * sigmamkpck[3][i];
		Ci[1][2][i] = Ci[1][2][i + 1] + mciipcii[2] * pci[1][i] - li[2][i] * sigmamkpck[1][i];
		Ci[2][2][i] = Ci[2][2][i + 1] - d33 - d11 + b33 + b11;
		Ci[3][2][i] = Ci[3][2][i + 1] + mciipcii[2] * pci[3][i] - li[2][i] * sigmamkpck[3][i];
		Ci[1][3][i] = Ci[1][3][i + 1] + mciipcii[3] * pci[1][i] - li[3][i] * sigmamkpck[1][i];
		Ci[2][3][i] = Ci[2][3][i + 1] + mciipcii[3] * pci[2][i] - li[3][i] * sigmamkpck[2][i];
		Ci[3][3][i] = Ci[3][3][i + 1] - d22 - d11 + b22 + b11;

		Di[1][i] = Di[1][i + 1] + mciipcii[1] - sigmamk[i] * li[1][i];
		Di[2][i] = Di[2][i + 1] + mciipcii[2] - sigmamk[i] * li[2][i];
		Di[3][i] = Di[3][i + 1] + mciipcii[3] - sigmamk[i] * li[3][i];
	}

	// calculating the coefficient vector of the invarient moment
	double fi[4][MAX_DimentionPlusOne];
	double gi[4][MAX_DimentionPlusOne];
	for (int i = NDOF; i > 0; i--)
	{
		fi[1][i] = (Ci[1][1][i] + Ei[1][1][i]) * Tgdahj[1][3][i] + (Ci[2][1][i] + Ei[2][1][i]) * Tgdahj[2][3][i]
			+ (Ci[3][1][i] + Ei[3][1][i]) * Tgdahj[3][3][i];
		fi[2][i] = (Ci[1][2][i] + Ei[1][2][i]) * Tgdahj[1][3][i] + (Ci[2][2][i] + Ei[2][2][i])
			* Tgdahj[2][3][i] + (Ci[3][2][i] + Ei[3][2][i]) * Tgdahj[3][3][i];
		fi[3][i] = (Ci[1][3][i] + Ei[1][3][i]) * Tgdahj[1][3][i] + (Ci[2][3][i] + Ei[2][3][i])
			* Tgdahj[2][3][i] + (Ci[3][3][i] + Ei[3][3][i]) * Tgdahj[3][3][i];

		gi[1][i] = Tgdahj[2][3][i] * Di[3][i] - Tgdahj[3][3][i] * Di[2][i];
		gi[2][i] = Tgdahj[3][3][i] * Di[1][i] - Tgdahj[1][3][i] * Di[3][i];
		gi[3][i] = Tgdahj[1][3][i] * Di[2][i] - Tgdahj[2][3][i] * Di[1][i];

	}
	// calculating the inertia tensor

	for (int i = 1; i < (NDOF + 1); i++)
	{
		for (int j = 1; j < i + 1; j++)
		{
			A_GDAHJ[i][j] = Tgdahj[1][3][j] * fi[1][i]
				+ Tgdahj[2][3][j] * fi[2][i] + Tgdahj[3][3][j] * fi[3][i] +
				kjpj[1][j] * gi[1][i]
				+ kjpj[2][j] * gi[2][i] + kjpj[3][j] * gi[3][i];
			A_GDAHJ[j][i] = A_GDAHJ[i][j];

		}
	}


}

void disp(double x)
{
	std::cout << x;
}

void disp(const char* input)
{
	std::cout << input << "\n";
}