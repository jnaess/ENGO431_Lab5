#include "Intersection.h"

Intersection::Intersection(){}

void Intersection::getImage(string filename) {
	ifstream file(filename);
	int i = 0;

	while (file.good()) {
		string line;

		//these are used to compare the type of input to read

		while (getline(file, line)) {

			imcoord.emplace_back();

			//splits the line from spaces into an array of three
			stringstream ssin(line);
			while (ssin.good()) {
				ssin >> imcoord[i].id;
				ssin >> imcoord[i].x;
				ssin >> imcoord[i].y;
			
			}
			i++;
		}


	}
	num_pt = imcoord.size() / 2;
	file.close();
}

void Intersection::getEOP(string filename) {
	ifstream file(filename);
	int i = 0;

	while (file.good()) {
		string line;

		//these are used to compare the type of input to read

		while (getline(file, line)) {

			eop.emplace_back();

			//splits the line from spaces into an array of three
			stringstream ssin(line);
			while (ssin.good()) {
				ssin >> eop[i].Xc;
				ssin >> eop[i].Yc;
				ssin >> eop[i].Zc;
				ssin >> eop[i].omega;
				ssin >> eop[i].phi;
				ssin >> eop[i].kappa;

				// CONVERTING EOP VALUES FROM DECIMAL DEGREES TO RADIANS
				eop[i].omega *= ((atan(1) * 4) / 180);
				eop[i].phi *= ((atan(1) * 4) / 180);
				eop[i].kappa *= ((atan(1) * 4) / 180);


				cout << eop[i].Xc << " " << eop[i].Yc << " " << eop[i].Zc << " " << eop[i].omega << " " << eop[i].phi << " " << eop[i].kappa << endl;
			}
			i++;
		}


	}

	file.close();

}


/*
void Intersection::getObject(string filename) {
	ifstream file(filename);
	int i = 0;

	while (file.good()) {
		string line;

		//these are used to compare the type of input to read

		while (getline(file, line)) {

			objcoord.emplace_back();

			//splits the line from spaces into an array of three
			stringstream ssin(line);
			while (ssin.good()) {
				ssin >> objcoord[i].id;
				ssin >> objcoord[i].X;
				ssin >> objcoord[i].Y;
				ssin >> objcoord[i].Z;
			}
			i++;
		}


	}

	file.close();
}
*/

void Intersection::approximate() {
	MatrixXd A(4, 3);
	int row = A.rows();
	A.setZero();
	int j = 0;

	MatrixXd b(4, 1);
	b.setZero();

	MatrixXd x(3, 1);
	x.setZero();

	
	results.resize(num_pt);
	// NEED TO CONVERT ANGLES TO RADIANS
	MatrixXd M_L(3, 3);
	MatrixXd M_R(3, 3);
	M_L = rotate(eop[0].kappa, 3) * rotate(eop[0].phi, 2) * rotate(eop[0].omega, 1);
	M_R = rotate(eop[1].kappa, 3) * rotate(eop[1].phi, 2) * rotate(eop[1].omega, 1);
	

	// Approximate X, Y, and Z for every point
	for (unsigned int i = 0; i < num_pt; i++)
	{
		// First two rows of A and b (x and y from left image)
		A(0, 0) = imcoord[i].x * M_L(2, 0) + c * M_L(0, 0);
		A(0, 1) = imcoord[i].x * M_L(2, 1) + c * M_L(0, 1);
		A(0, 2) = imcoord[i].x * M_L(2, 2) + c * M_L(0, 2);
		b(0, 0) = A(0, 0) * eop[0].Xc + A(0, 1) * eop[0].Yc + A(0, 2) * eop[0].Zc;

		A(1, 0) = imcoord[i].y * M_L(2, 0) + c * M_L(1, 0);
		A(1, 1) = imcoord[i].y * M_L(2, 1) + c * M_L(1, 1);
		A(1, 2) = imcoord[i].y * M_L(2, 2) + c * M_L(1, 2);
		b(1, 0) = A(1, 0) * eop[0].Xc + A(1, 1) * eop[0].Yc + A(1, 2) * eop[0].Zc;

		// Last two rows of A and b (x and y from right image)

		A(2, 0) = imcoord[num_pt + i].x * M_R(2, 0) + c * M_R(0, 0);
		A(2, 1) = imcoord[num_pt + i].x * M_R(2, 1) + c * M_R(0, 1);
		A(2, 2) = imcoord[num_pt + i].x * M_R(2, 2) + c * M_R(0, 2);
		b(2, 0) = A(2, 0) * eop[1].Xc + A(2, 1) * eop[1].Yc + A(2, 2) * eop[1].Zc;

		A(3, 0) = imcoord[num_pt + i].y * M_R(2, 0) + c * M_R(1, 0);
		A(3, 1) = imcoord[num_pt + i].y * M_R(2, 1) + c * M_R(1, 1);
		A(3, 2) = imcoord[num_pt + i].y * M_R(2, 2) + c * M_R(1, 2);
		b(3, 0) = A(3, 0) * eop[1].Xc + A(3, 1) * eop[1].Yc + A(3, 2) * eop[1].Zc;

		x = (A.transpose() * A).inverse() * A.transpose() * b;

		results[i].X = x(0, 0);
		results[i].Y = x(1, 0);
		results[i].Z = x(2, 0);
		cout << results[i].X << " " << results[i].Y << " " << results[i].Z << endl;
	}
	//x_hat = x;
}

// TAKES ANGLE AS RADIANS
MatrixXd Intersection::rotate(double angle, int axis) {
	MatrixXd R(3, 3);
	if (axis == 1) {
		R << 1, 0, 0,
			0, cos(angle), sin(angle),
			0, -sin(angle), cos(angle);
	}
	else if (axis == 2) {
		R << cos(angle), 0, -sin(angle),
			0, 1, 0,
			sin(angle), 0, cos(angle);
	}
	else if (axis == 3) {
		R << cos(angle), sin(angle), 0,
			-sin(angle), cos(angle), 0,
			0, 0, 1;
	}
	else cout << "invalid axis" << endl;

	return R;

}


void Intersection::update(int point_index) {
	results[point_index].X = x_hat(0, 0);
	results[point_index].Y = x_hat(1, 0);
	results[point_index].Z = x_hat(2, 0);
}


void Intersection::designA(int point_index) {
	cout << "\n\n----- Now adjusting point: " << imcoord[point_index].id << " -----\n";
	
	A.resize(4, 3);
	w.resize(4, 1);
	
	MatrixXd M_L(3, 3);
	MatrixXd M_R(3, 3);
	M_L = rotate(eop[0].kappa, 3) * rotate(eop[0].phi, 2) * rotate(eop[0].omega, 1);
	M_R = rotate(eop[1].kappa, 3) * rotate(eop[1].phi, 2) * rotate(eop[1].omega, 1);


	double dxL = results[point_index].X - eop[0].Xc;
	double dyL = results[point_index].Y - eop[0].Yc;
	double dzL = results[point_index].Z - eop[0].Zc;

	double dxR = results[point_index].X - eop[1].Xc;
	double dyR = results[point_index].Y - eop[1].Yc;
	double dzR = results[point_index].Z - eop[1].Zc;
	//cout << Xc << " " << Yc << " " << Zc << endl;
	//cout << "dx = " << dxL << endl;
	//cout << "dy = " << dyL << endl;
	//cout << "dz = " << dzL << endl;
	
	// First two rows of design matrix (left image point)
	double UL = M_L(0, 0) * dxL + M_L(0, 1) * dyL + M_L(0, 2) * dzL;
	double VL = M_L(1, 0) * dxL + M_L(1, 1) * dyL + M_L(1, 2) * dzL;
	double WL = M_L(2, 0) * dxL + M_L(2, 1) * dyL + M_L(2, 2) * dzL;
	double cWL = -c / (WL * WL);

	A(0, 0) = -cWL * (M_L(2, 0) * UL - M_L(0, 0) * WL);
	A(0, 1) = -cWL * (M_L(2, 1) * UL - M_L(0, 1) * WL);
	A(0, 2) = -cWL * (M_L(2, 2) * UL - M_L(0, 2) * WL);
	A(1, 0) = cWL * (M_L(2, 0) * VL - M_L(1, 0) * WL);
	A(1, 1) = cWL * (M_L(2, 1) * VL - M_L(1, 1) * WL);
	A(1, 2) = cWL * (M_L(2, 2) * VL - M_L(1, 2) * WL);

	// Last two rows of design matrix (right image point)
	double UR = M_R(0, 0) * dxR + M_R(0, 1) * dyR + M_R(0, 2) * dzR;
	double VR = M_R(1, 0) * dxR + M_R(1, 1) * dyR + M_R(1, 2) * dzR;
	double WR = M_R(2, 0) * dxR + M_R(2, 1) * dyR + M_R(2, 2) * dzR;
	double cWR = -c / (WR * WR);

	A(2, 0) = -cWR * (M_R(2, 0) * UR - M_R(0, 0) * WR);
	A(2, 1) = -cWR * (M_R(2, 1) * UR - M_R(0, 1) * WR);
	A(2, 2) = -cWR * (M_R(2, 2) * UR - M_R(0, 2) * WR);
	A(3, 0) = cWR * (M_R(2, 0) * VR - M_R(1, 0) * WR);
	A(3, 1) = cWR * (M_R(2, 1) * VR - M_R(1, 1) * WR);
	A(3, 2) = cWR * (M_R(2, 2) * VR - M_R(1, 2) * WR);

	w(0, 0) = -(c * UL / WL) - imcoord[point_index].x;
	w(1, 0) = -(c * VL / WL) - imcoord[point_index].y;
	w(2, 0) = -(c * UR / WR) - imcoord[point_index+num_pt].x;
	w(3, 0) = -(c * VR / WR) - imcoord[point_index+num_pt].y;
	cout << A << endl;
	cout << "misclosure" << endl << w << endl;
	
}

void Intersection::getxhat(int point_index) {
	
	double S = ceil((eop[0].Zc - results[point_index].Z) /c*1000);
	double tol_c= S * stdv / 10000.0;
	MatrixXd tol(1, 2);
	tol << stdv / (10.0*c),stdv / (10.0 * 162.0);
	double tol_min = tol.minCoeff();
	x_hat.resize(3, 1);
	x_hat(0, 0) = results[point_index].X;
	x_hat(1, 0) = results[point_index].Y;
	x_hat(2, 0) = results[point_index].Z;
	cout << "tolerance is  " << tol_c << endl;
	//cout << "tolerance is  " << tol_c << "    "<< tol <<endl;;
	MatrixXd P(4, 4);
	P.setIdentity();
	P = P / (stdv*stdv);

	MatrixXd delta(3, 1);
	delta = -(A.transpose()*P * A).inverse() * A.transpose() * P*w;
	cout <<"delta  " << endl << delta << endl;
	x_hat = x_hat + delta;
	//cout << "final xhat" << endl << x_hat << endl;

	vector <double> d;
	for (int i= 0; i < delta.size(); i++) {
		d.push_back(abs(delta(i, 0)));
	}
	double min1 = *max_element(d.begin(),d.begin() + 3);
	//double min2 = *max_element(d.begin()+3, d.end());
	//cout << " min   " << min1 << "min2   " << min2<<endl;
	cout << " min   " << min1 <<  endl;
	if (tol_c < min1){
		count++;
		cout << "counter  " << count << endl;
	}
	else {
		criteria = true;
		cout << "final xhat " << endl << x_hat << endl;
		results[point_index].v_hat = A * delta + w;
		cout << "residuals " << endl << results[point_index].v_hat << endl;
		MatrixXd C_xhat = (A.transpose() * P * A).inverse();
		cout << "C_xhat " << endl << C_xhat << endl;
		MatrixXd C_vhat = A * C_xhat * A.transpose();
		MatrixXd I(4, 4);
		I.setIdentity();
		cout << "C_vhat " << endl << C_vhat << endl;
		cout << "R " << endl << I - (C_vhat * P) << endl;

	}
	/*
	if (tol_c < min1 || tol_min < min2) {
		count++;
		cout << "counter  " << count << endl;
	}
	else {
		criteria = true;
		cout << "final xhat" << endl << x_hat << endl;
	}
	*/
	//exit(1);
}

void Intersection::OutputResults(const string& filename)
{
	// Check if the output file can be opened
	std::ofstream outFile(filename);
	if (outFile.fail())
	{
		std::cout << "Could not open output file" << std::endl;
		exit(1);
	}

	MatrixXd V_RMSE(1, 4);
	V_RMSE.setZero();

	double Z_ave = 0;
	VectorXd v_hattemp;
	cout << "\n\n--- Final Results, Coordinates and Residuals (X [m], Y[m], Z[m], V_xL [mm], V_yL [mm], V_xR [mm], V_yR [mm], RMSE [mm]) ---\n";
	for (int i = 0; i < num_pt; i++)
	{
		cout << setprecision(0) << fixed << imcoord[i].id << setprecision(8) << "  " << results[i].X << "  " << results[i].Y << "  " << results[i].Z << "  ";
		cout << results[i].v_hat(0, 0) << "  " << results[i].v_hat(1, 0) << "  " << results[i].v_hat(2, 0) << "  " << results[i].v_hat(3, 0) << "  ";
		v_hattemp = results[i].v_hat;
		cout << sqrt(v_hattemp.dot(v_hattemp) / 4) << endl << setprecision(6);
		cout.unsetf(ios_base::floatfield);


		for (int j = 0; j < 4; j++) 
		{
			V_RMSE(0, j) += pow(results[i].v_hat(j, 0), 2);
		}
		




		outFile << setprecision(0) << fixed << imcoord[i].id << setprecision(8) << "  " << results[i].X << "  " << results[i].Y << "  " << results[i].Z << "  ";
		outFile << results[i].v_hat(0, 0) << "  " << results[i].v_hat(1, 0) << "  " << results[i].v_hat(2, 0) << "  " << results[i].v_hat(3, 0) << "  ";
		outFile << sqrt(v_hattemp.dot(v_hattemp) / 4) << endl;
		Z_ave += results[i].Z;
		
	}
	Z_ave /= num_pt;
	
	for (int j = 0; j < 4; j++)
	{
		V_RMSE(0, j) = sqrt(V_RMSE(0, j) / num_pt);
	}
	cout << "Average height = " << Z_ave << endl;
	cout << "Residual RMSE (x_L, y_L, x_R, y_R) = " << V_RMSE << endl;
}