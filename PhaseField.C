#include <iostream>
#include <stdio.h>
#include <math.h>
#include "PhaseField.h"

int main() {

	int N, NN;

	N = 128;	// # of grid points in each direction
	NN = N * N;

	double *phi_new, *phi, *mu;

	phi_new = new double[NN];
	phi = new double[NN];
	mu = new double[NN];

	Init_RND(NN, 0.3, phi);

	double *phi_xx, *phi_yy, *mu_xx, *mu_yy;

	phi_xx = new double[NN];
	mu_xx = new double[NN];

	phi_yy = new double[NN];
	mu_yy = new double[NN];

	double *bound_b_o, *bound_r_o, *bound_u_o, *bound_l_o;
	double *bound_b_i, *bound_r_i, *bound_u_i, *bound_l_i;

	bound_b_o = new double[N];
	bound_r_o = new double[N];
	bound_u_o = new double[N];
	bound_l_o = new double[N];
	bound_b_i = new double[N];
	bound_r_i = new double[N];
	bound_u_i = new double[N];
	bound_l_i = new double[N];

	double gp_x[2];
	double gp_y[2];

	double dt = 2.0e-3;
	double dg, laplacian_phi, laplacian_mu;
	int max_nt = 250000;

//	double counter = dt * max_nt;

	for (int nt = 1; nt < max_nt + 1; nt++) {

		for (int n = 0; n < N; n++) {

			bound_b_o[n] = phi[n];
			bound_r_o[n] = phi[(n + 1) * N - 1];
			bound_u_o[n] = phi[n + N * (N - 1)];
			bound_l_o[n] = phi[n * N];

			bound_b_i[n] = phi[n + N];
			bound_r_i[n] = phi[(n + 1) * N - 2];
			bound_u_i[n] = phi[n + N * (N - 2)];
			bound_l_i[n] = phi[n * N + 1];

		}

//			Computing second order derivatives

//########################################################
//				  1. x direction
//########################################################

		for (int j = 0; j < N; j++) {

			int k = j * N;

			gp_x[0] = bound_r_o[j];
			gp_x[1] = phi[k + 1];

			d_dxx(phi[k], gp_x, &(phi_xx[k]));

		}

		for (int i = 1; i < N - 1; i++) {

			for (int j = 0; j < N; j++) {

				int k = i + j * N;

				gp_x[0] = phi[k - 1];
				gp_x[1] = phi[k + 1];

				d_dxx(phi[k], gp_x, &(phi_xx[k]));
			}
		}

		for (int j = 0; j < N; j++) {

			int k = (N - 1) + j * N;

			gp_x[0] = phi[k - 1];
			gp_x[1] = bound_l_o[j];

			d_dxx(phi[k], gp_x, &(phi_xx[k]));

		}
//########################################################
//					2. y direction
//########################################################

		for (int i = 0; i < N; i++) {

			int k = i;

			gp_y[0] = bound_u_o[i];
			gp_y[1] = phi[k + N];

			d_dxx(phi[k], gp_y, &(phi_yy[k]));

		}

		for (int i = 0; i < N; i++) {

			for (int j = 1; j < N - 1; j++) {

				int k = i + j * N;

				gp_y[0] = phi[k - N];
				gp_y[1] = phi[k + N];

				d_dxx(phi[k], gp_y, &(phi_yy[k]));
			}
		}

		for (int i = 0; i < N; i++) {

			int k = i + (N - 1) * N;

			gp_y[0] = phi[k - N];
			gp_y[1] = bound_b_o[i];

			d_dxx(phi[k], gp_y, &(phi_yy[k]));

		}

//########################################################


		for (int np = 0; np < NN; np++) {

			dg = 64.0 * phi[np] * phi[np] * phi[np] - 96.0 * phi[np] * phi[np] + 32.0 * phi[np];

			laplacian_phi = phi_xx[np] + phi_yy[np];

			mu[np] = dg - 10.0 * laplacian_phi;

		}


//########################################################
//########################################################
//########################################################


		for (int n = 0; n < N; n++) {

			bound_b_o[n] = mu[n];
			bound_r_o[n] = mu[(n + 1) * N - 1];
			bound_u_o[n] = mu[n + N * (N - 1)];
			bound_l_o[n] = mu[n * N];

			bound_b_i[n] = mu[n + N];
			bound_r_i[n] = mu[(n + 1) * N - 2];
			bound_u_i[n] = mu[n + N * (N - 2)];
			bound_l_i[n] = mu[n * N + 1];

		}

//			Computing second order derivatives

//########################################################
//				  1. x direction
//########################################################

		for (int j = 0; j < N; j++) {

			int k = j * N;

			gp_x[0] = bound_r_o[j];
			gp_x[1] = mu[k + 1];

			d_dxx(mu[k], gp_x, &(mu_xx[k]));

		}

		for (int i = 1; i < N - 1; i++) {

			for (int j = 0; j < N; j++) {

				int k = i + j * N;

				gp_x[0] = mu[k - 1];
				gp_x[1] = mu[k + 1];

				d_dxx(mu[k], gp_x, &(mu_xx[k]));
			}
		}

		for (int j = 0; j < N; j++) {

			int k = (N - 1) + j * N;

			gp_x[0] = mu[k - 1];
			gp_x[1] = bound_l_o[j];

			d_dxx(mu[k], gp_x, &(mu_xx[k]));

		}
//########################################################
//					2. y direction
//########################################################

		for (int i = 0; i < N; i++) {

			int k = i;

			gp_y[0] = bound_u_o[i];
			gp_y[1] = mu[k + N];

			d_dxx(mu[k], gp_y, &(mu_yy[k]));

		}

		for (int i = 0; i < N; i++) {

			for (int j = 1; j < N - 1; j++) {

				int k = i + j * N;

				gp_y[0] = mu[k - N];
				gp_y[1] = mu[k + N];

				d_dxx(mu[k], gp_y, &(mu_yy[k]));
			}
		}

		for (int i = 0; i < N; i++) {

			int k = i + (N - 1) * N;

			gp_y[0] = mu[k - N];
			gp_y[1] = bound_b_o[i];

			d_dxx(mu[k], gp_y, &(mu_yy[k]));

		}

//########################################################
//########################################################
//########################################################


		for (int np = 0; np < NN; np++) {

			laplacian_mu = mu_xx[np] + mu_yy[np];

			phi_new[np] = phi[np] + dt * laplacian_mu;

		}

		for (int np = 0; np < NN; np++) {

			phi[np] = phi_new[np];

		}

	}

	int np = 0;

	FILE * fp;
	fp = fopen("phi.m","w+");
	fprintf(fp, "clear\nclc\nphi1 = [");

	for (int i = 0; i < N; i++) {

		for (int j = 0; j < N; j++) {

			np = i + N * j;

			if (np == NN - 1)
				fprintf(fp, "%.15f", phi_new[np]);
			else if (j == N - 1)
				fprintf(fp, "%.15f\n", phi_new[np]);
			else
				fprintf(fp, "%.15f, ", phi_new[np]);

		}
	}

	fprintf(fp, "];\nfigure\ncontour(phi1,'fill','on')");
	fclose(fp);

	/* free the memory */
	delete[] phi_new;
	delete[] phi;
	delete[] mu;

	delete[] bound_b_o;
	delete[] bound_r_o;
	delete[] bound_u_o;
	delete[] bound_l_o;
	delete[] bound_b_i;
	delete[] bound_r_i;
	delete[] bound_u_i;
	delete[] bound_l_i;

	delete[] phi_xx;
	delete[] mu_xx;

	delete[] phi_yy;
	delete[] mu_yy;

	return 0;

}

