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

	MatrixXd M_L(3, 3);
	MatrixXd M_R(3, 3);
	x_hat.resize(num_pt);
	// NEED TO CONVERT ANGLES TO RADIANS
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

		x_hat[i].X = x(0, 0);
		x_hat[i].Y = x(1, 0);
		x_hat[i].Z = x(2, 0);
		cout << x_hat[i].X << " " << x_hat[i].Y << " " << x_hat[i].Z << endl;
	}
	
}

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

/*
void Intersection::approximate() {
	MatrixXd A(num_pt * 2, 4);
	int row = A.rows();
	A.setZero();
	int j = 0;
	for (unsigned int i = 0; i < A.rows() - 1; i += 2) {
		A(i, 0) = imcoord[j].x;
		A(i, 1) = -imcoord[j].y;
		A(i, 2) = 1;
		A(i + 1, 0) = imcoord[j].y;
		A(i + 1, 1) = imcoord[j].x;
		A(i + 1, 3) = 1;
		j++;
	}
	//cout << A;
	MatrixXd param(4, 1); // a,b,tx,ty
	MatrixXd l(row, 1);
	int k = 0;

	for (unsigned int i = 0; i < row - 1; i += 2) {
		l(i, 0) = objcoord[k].X;
		l(i + 1, 0) = objcoord[k].Y;
		z_ave += objcoord[k].Z;
		k++;
	}
	z_ave = z_ave / num_pt;
	param = (A.transpose() * A).inverse() * A.transpose() * l;
	//cout << "obs"<<l<<endl<<"param"  <<param << endl;
	double lamda = 0;
	lamda = sqrt(pow(param(0, 0), 2) + pow(param(1, 0), 2));
	Xc = param(2, 0);
	Yc = param(3, 0);
	Zc = c * lamda + z_ave;
	kappa = atan2(param(1, 0), param(0, 0));

	//cout << "param    " << Xc << "    "<< Yc << "    " << Zc << "    " << kappa << endl;
	xhat.resize(6, 1);
	xhat(0, 0) = Xc;
	xhat(1, 0) = Yc;
	xhat(2, 0) = Zc;
	xhat(3, 0) = omega;
	xhat(4, 0) = phi;
	xhat(5, 0) = kappa;
	//cout << "'approximate" << xhat << endl;
}
*/

/*
void Intersection::update() {
	Xc = xhat(0, 0);
	Yc = xhat(1, 0);
	Zc = xhat(2, 0);
	omega = xhat(3, 0);
	phi = xhat(4, 0);
	kappa = xhat(5, 0);
}


void Intersection::designA() {
	A.resize(num_pt * 2, 3);
	w.resize(num_pt * 2, 1);
	MatrixXd M(3,3);
	cout << omega << " " << phi << " " << kappa << -18.854 * ((atan(1) * 4) / 180) <<  endl;
	cout << kappa - (-18.854 * ((atan(1) * 4) / 180)) << endl;
	kappa = -18.854 * ((atan(1) * 4) / 180);

	M = rotate(kappa, 3) * rotate(phi, 2) * rotate(omega, 1);
	
	
	int j = 0;
	for (unsigned int i = 0; i < num_pt*2 - 1; i += 2) {
		Xc = 6338.6;
		Yc = 3984.6;
		Zc = 1453.1;
		//kappa = -18.854 * ((atan(1) * 4) / 180);
		double dx = objcoord[j].X - Xc;
		double dy = objcoord[j].Y - Yc;
		double dz = objcoord[j].Z - Zc;
		//cout << Xc << " " << Yc << " " << Zc << endl;

		double U = M(0, 0) * dx + M(0, 1) * dy + M(0, 2) * dz;
		double W = M(2, 0) * dx + M(2, 1) * dy + M(2, 2) * dz;
		double V = M(1, 0) * dx + M(1, 1) * dy + M(1, 2) * dz;
		double cW = -c / (W * W);
		
		A(i, 0) = cW * (M(2, 0) * U - M(0, 0) * W);
		A(i, 1) = cW * (M(2, 1) * U - M(0, 1) * W);
		A(i, 2) = cW * (M(2, 2) * U - M(0, 2) * W);
		//A(i, 3) = cW * (dy * (U * M(2, 2) - W * M(0, 2)) - dz*(U * M(2, 1) - W * M(0, 1)));
		//A(i, 3) = dy * A(i, 2) - dz * A(i, 1);
		//A(i, 4) = cW * (dx * (-W * sin(phi) * cos(kappa) - U * cos(phi)) + dy * (W * sin(omega) * cos(phi) * cos(kappa) - U * sin(omega) * sin(phi))
			//+ dz * (-W * cos(omega) * cos(phi) * cos(kappa) + U * cos(omega) * sin(phi)));
		//A(i, 5) = -c * (V / W);

		A(i+1, 0) = -cW * (M(2, 0) * V- M(1, 0) * W);
		A(i+1, 1) = -cW * (M(2, 1) * V - M(1, 1) * W);
		A(i+1, 2) = -cW * (M(2, 2) * V - M(1, 2) * W);
		//A(i+1, 3) = cW * (dy * (V * M(2, 2) - W * M(1, 2)) - dz * (V * M(2, 1) - W * M(1, 1)));
		//A(i+1, 4) = cW * (dx * (W * sin(phi) * sin(kappa) - V * cos(phi)) + dy * (-W * sin(omega) * cos(phi)*sin(kappa) - V * sin(omega)*sin(phi))
			//+ dz * (W * cos(omega) * cos(phi) * sin(kappa) + V * cos(omega) * sin(phi)));
		//A(i+1, 5) = c * U / W;
	
		// VV WHAT IS UP WITH THIS VV
		//w(i, 0) = -0.006 - (c * U / W) - imcoord[j].x;
		//w(i + 1, 0) = 0.006 - (c * V / W) - imcoord[j].y;
		w(i, 0) = -(c * U / W) - imcoord[j].x;
		w(i + 1, 0) = -(c * V / W) - imcoord[j].y;
		j++;
	}

	cout << A<< endl;
	cout <<"misclosure" << endl <<w << endl;

}

void Intersection::getxhat() {
	double stdv = 0.015;
	double S = ceil((Zc - z_ave) /c*1000);
	double tol_c= S * stdv / 10000;
	MatrixXd tol(1, 2);
	tol << stdv / (10*c),stdv / (10 * 162);
	double tol_min = tol.minCoeff();
	

	cout << "tolerance is  " << tol_c << "    "<< tol <<endl;;
	MatrixXd P(num_pt*2, num_pt*2);
	P.setIdentity();
	P = P / (stdv/stdv);

	MatrixXd delta(6, 1);
	delta = -(A.transpose()*P * A).inverse() * A.transpose() * P*w;
	cout <<"delta  " << endl << delta << endl;
	xhat = xhat + delta;

	vector <double> d;
	for (int i= 0; i < delta.size(); i++) {
		d.push_back(abs(delta(i, 0)));
	}
	double min1 = *max_element(d.begin(),d.begin() + 3);
	double min2 = *max_element(d.begin()+3, d.end());
	cout << " min   " << min1 << "min2   " << min2<<endl;
	if (tol_c < min1 || tol_min < min2){
		count++;
		cout << "counter  " << count << endl;
	}
	else {
		criteria = true;
		cout << "final xhat" << endl << xhat << endl;
	}
	
	//exit(1);
}


*/