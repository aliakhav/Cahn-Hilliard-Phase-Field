#include "PhaseField.h"
#include <cstdlib>
#include <time.h>

void Init_RND(int n, double x, double rnd_init[]){

	int rnd;

	x = x - 0.01;

	time_t seconds;
	time(&seconds);
	srand((unsigned int) seconds);

	for (int i = 0; i < n; i++) {

		rnd = rand() % 1000;
		rnd_init[i] = x + rnd * 2.0e-5;

	}

	return;
}
