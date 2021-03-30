
#include "Resection.h"

int main() {
	Resection res;
	res.getObject("object.txt");
	res.getImage("image.txt");
	res.approximate();
	res.designA();
	return 0;
}