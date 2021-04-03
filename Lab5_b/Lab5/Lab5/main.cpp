
#include "Intersection.h"

int main() {
	Intersection sample_int;
	sample_int.getImage("sample_image_b.txt");
	sample_int.getEOP("sample_eop.txt");
	sample_int.approximate();

	for (int i = 0; i < sample_int.num_pt; i++)
	{

		//sample_int.designA(i);
		//sample_int.getxhat(i);
		//sample_int.update(i);
		
		while (!sample_int.criteria) {
			sample_int.designA(i);
			sample_int.getxhat(i);
			sample_int.update(i);
		}
		sample_int.criteria = false;
		
	}
	
	return 0;
}