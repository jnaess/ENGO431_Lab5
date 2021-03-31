
#include "Resection.h"

int main() {
	Resection res;
	res.getObject("object.txt");
	res.getImage("image.txt");
	res.approximate();
	while (!res.criteria) {
		res.designA();
		res.getxhat();
		res.update();
	}
	return 0;
}