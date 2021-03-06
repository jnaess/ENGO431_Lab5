#pragma once

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace Eigen;
using namespace std;

struct Model {

	double id = 0,x = 0, y = 0, X = 0, Y = 0, Z= 0;
};
struct EOP_model {

	double Xc, Yc, Zc, omega, phi, kappa;
};
struct param {
	double X, Y, Z;
	MatrixXd v_hat;
};

class Intersection {

public:
	vector <Model> imcoord;
	vector <EOP_model> eop;
	vector <param> results;
	//vector <Model> objcoord;
	MatrixXd x_hat;
	MatrixXd A, w;
	int num_pt = 0;
	int count = 0;
	//double z_ave = 0;
	//double Xc, Yc, Zc, omega, phi, kappa;
	double c;
	double stdv;
	//double stdv = 0.015; // sample data
	//double stdv = 0.006; // our data
	//double c = 152.15; // sample data
	//double c = 153.358; // our data
	bool criteria = false;
	Intersection();
	void getImage(string filename);
	void getEOP(string filename);
	//void getObject(string filename);
	void approximate();
	void update(int point_index);
	MatrixXd rotate(double angle, int axis);
	void designA(int point_index);
	void getxhat(int point_index);
	void OutputResults(const string& filename);
};
