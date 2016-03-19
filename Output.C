#include "PhaseField.h"

void Output_(int N, int NN, int step, double phi[]) {

	int np = 0;

	char filename[50];
	sprintf(filename, "solution_%06d.m", step);

	FILE * fp;
	fp = fopen(filename,"w+");
	fprintf(fp, "clear\nclc\nphi1 = [");

	for (int i = 0; i < N; i++) {

		for (int j = 0; j < N; j++) {

			np = i + N * j;

			if (np == NN - 1)
				fprintf(fp, "%.15f", phi[np]);
			else if (j == N - 1)
				fprintf(fp, "%.15f\n", phi[np]);
			else
				fprintf(fp, "%.15f, ", phi[np]);

		}
	}

	fprintf(fp, "];\nfigure\ncontour(phi1,'fill','on')");
	fclose(fp);

	return;
}

