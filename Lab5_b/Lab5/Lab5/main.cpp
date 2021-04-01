
#include "Intersection.h"

int main() {
	Intersection res;
	res.getImage("sample_image_b.txt");
	res.getEOP("sample_eop.txt");
	res.approximate();

	/*
	while (!res.criteria) {
		res.designA();
		res.getxhat();
		res.update();
	}

	*/
	return 0;
}