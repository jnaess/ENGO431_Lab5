#pragma once

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

using namespace Eigen;
using namespace std;

struct Model {

	double id = 0,x = 0, y = 0, X = 0, Y = 0, Z= 0;
};

class Resection {

public:
	vector <Model> imcoord;
	vector <Model> objcoord;
	MatrixXd xhat, A;
	int num_pt = 0;
	double Xc, Yc, Zc, omega=0, phi=0, kappa;
	double c = 152.15; // sample data
	//double c = 153.358; // our data
	Resection();
	void getImage(string filename);
	void getObject(string filename);
	void approximate();
	void update();
	MatrixXd rotate(double angle, int axis);
	void designA();
};
