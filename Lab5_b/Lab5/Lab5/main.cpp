
#include "Intersection.h"

int main() {
	cout << "---------- Performing Intersection on Sample Data ----------\n\n";
	Intersection sample_int;
	sample_int.stdv = 0.015;
	sample_int.c = 152.15;
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
	sample_int.OutputResults("sample_results.txt");
	
	cout << "\n\n\n---------- Performing Intersection on Our Data ----------\n\n";

	Intersection our_int;
	our_int.stdv = 0.006;
	our_int.c = 152.358;
	our_int.getImage("image_b.txt");
	our_int.getEOP("eop.txt");
	our_int.approximate();

	for (int i = 0; i < our_int.num_pt; i++)
	{

		//our_int.designA(i);
		//our_int.getxhat(i);
		//our_int.update(i);

		while (!our_int.criteria) {
			our_int.designA(i);
			our_int.getxhat(i);
			our_int.update(i);
		}
		our_int.criteria = false;

	}
	our_int.OutputResults("final_results.txt");
	return 0;
}