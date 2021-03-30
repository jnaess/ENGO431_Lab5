#include "Resection.h"

Resection::Resection(){}

void Resection::getImage(string filename) {
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
	num_pt = imcoord.size();
	file.close();
}

void Resection::getObject(string filename) {
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

void Resection::approximate() {
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
	double z_ave = 0;
	for (unsigned int i = 0; i < row - 1; i += 2) {
		l(i, 0) = objcoord[k].X;
		l(i + 1, 0) = objcoord[k].Y;
		z_ave += objcoord[k].Z;
		k++;
	}
	z_ave = z_ave /num_pt;
	param = (A.transpose() * A).inverse() * A.transpose() * l;
	//cout << "obs"<<l<<endl<<"param"  <<param << endl;
	double lamda = 0;
	lamda = sqrt(pow(param(0, 0), 2) + pow(param(1, 0), 2));
	Xc = param(2, 0);
	Yc = param(3, 0);
	Zc = c * lamda + z_ave;
	kappa = atan(param(1, 0) / param(0, 0));
	
}
void Resection::update() {
	Xc = xhat(0, 0);
	Yc = xhat(1, 0);
	Zc = xhat(2, 0);
	omega = xhat(3, 0);
	phi = xhat(4, 0);
	kappa = xhat(5, 0);
}
MatrixXd Resection::rotate(double angle, int axis) {
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

void Resection::designA() {
	A.resize(num_pt * 2, 6);
	MatrixXd M(3,3);
	M = rotate(kappa, 3) * rotate(phi, 2) * rotate(omega, 1);
	int j = 0;
	for (unsigned int i = 0; i < num_pt*2 - 1; i += 2) {
		double dx = objcoord[j].X - Xc;
		double dy = objcoord[j].Y - Yc;
		double dz = objcoord[j].Z - Zc;
	

		double U = M(0, 0) * dx + M(0, 1) * dy + M(0, 2) * dz;
		double W = M(2, 0) * dx + M(2, 1) * dy + M(2, 2) * dz;
		double V = M(1, 0) * dx + M(1, 1) * dy + M(1, 2) * dz;
		double cW = -c / (W * W);

		A(i, 0) = cW * (M(2, 0) * U - M(0, 0) * W);
		A(i, 1) = cW * (M(2, 1) * U - M(0, 1) * W);
		A(i, 2) = cW * (M(2, 2) * U - M(0, 2) * W);
		A(i, 3) = cW * (dy * (U * M(2, 2) - W * M(0, 2)) - dz*(U * M(2, 1) - W * M(0, 1)));
		A(i, 4) = cW* (dx * (-W * sin(phi) * cos(kappa) - U * cos(phi)) + dy * (W * sin(omega) * cos(phi)*cos(kappa)- U * sin(omega) * sin(phi))
			+ dz * (-W * cos(omega) * cos(phi) * cos(kappa) + U * cos(omega) * sin(phi)));
		A(i, 5) = -c * V / W;

		A(i+1, 0) = cW * (M(2, 0) * V- M(1, 0) * W);
		A(i+1, 1) = cW * (M(2, 1) * V - M(1, 1) * W);
		A(i+1, 2) = cW * (M(2, 2) * V - M(1, 2) * W);
		A(i+1, 3) = cW * (dy * (V * M(2, 2) - W * M(1, 2)) - dz * (V * M(2, 1) - W * M(1, 1)));
		A(i+1, 4) = cW * (dx * (-W * sin(phi) * sin(kappa) - V * cos(phi)) + dy * (-W * sin(omega) * cos(phi)*sin(kappa) - V * sin(omega)*sin(phi))
			+ dz * (W * cos(omega) * cos(phi) * sin(kappa) + V * cos(omega) * sin(phi)));
		A(i+1, 5) = c * U / W;
		j++;
	}
	cout << A;

}