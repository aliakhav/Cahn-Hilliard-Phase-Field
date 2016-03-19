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

	double concent = 0.7;

	int cc;

	if (concent == 0.3)
		cc = 1;
	else if (concent == 0.5)
		cc = 2;
	else
		cc = 3;

	Init_RND(NN, concent, phi);

	double *phi_xx, *phi_yy, *mu_xx, *mu_yy;
	double *phi_x, *phi_y;

	phi_xx = new double[NN];
	mu_xx = new double[NN];

	phi_yy = new double[NN];
	mu_yy = new double[NN];

	phi_x = new double[NN];
	phi_y = new double[NN];

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
	int max_nt = 250000, nt;

	double t = dt * max_nt;

	double *mix_E, *Interface_E, *r_c;

	mix_E = new double[max_nt];
	Interface_E = new double[max_nt];
	r_c = new double[max_nt];


	Output_(N, NN, 1, phi);

	char filename[50];
	sprintf(filename, "Records_%03d.m", cc);

	FILE * fp1;
	fp1 = fopen(filename,"w+");
	fprintf(fp1, "clear\nclc\nrecs = [");


//################################################################
//################################################################
//						Time Stepping
//################################################################
//################################################################


//					Applying B.C. for phi

	for (nt = 1; nt < max_nt + 1; nt++) {

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

//	  Computing 1st and 2nd order derivatives of phi

//########################################################
//				  1. x direction
//########################################################

		for (int j = 0; j < N; j++) {

			int k = j * N;

			gp_x[0] = bound_r_o[j];
			gp_x[1] = phi[k + 1];

			d_dx(gp_x, &(phi_x[k]));
			d_dxx(phi[k], gp_x, &(phi_xx[k]));

		}

		for (int i = 1; i < N - 1; i++) {

			for (int j = 0; j < N; j++) {

				int k = i + j * N;

				gp_x[0] = phi[k - 1];
				gp_x[1] = phi[k + 1];

				d_dx(gp_x, &(phi_x[k]));
				d_dxx(phi[k], gp_x, &(phi_xx[k]));
			}
		}

		for (int j = 0; j < N; j++) {

			int k = (N - 1) + j * N;

			gp_x[0] = phi[k - 1];
			gp_x[1] = bound_l_o[j];

			d_dx(gp_x, &(phi_x[k]));
			d_dxx(phi[k], gp_x, &(phi_xx[k]));

		}
//########################################################
//					2. y direction
//########################################################

		for (int i = 0; i < N; i++) {

			int k = i;

			gp_y[0] = bound_u_o[i];
			gp_y[1] = phi[k + N];

			d_dx(gp_y, &(phi_y[k]));
			d_dxx(phi[k], gp_y, &(phi_yy[k]));

		}

		for (int i = 0; i < N; i++) {

			for (int j = 1; j < N - 1; j++) {

				int k = i + j * N;

				gp_y[0] = phi[k - N];
				gp_y[1] = phi[k + N];

				d_dx(gp_y, &(phi_y[k]));
				d_dxx(phi[k], gp_y, &(phi_yy[k]));
			}
		}

		for (int i = 0; i < N; i++) {

			int k = i + (N - 1) * N;

			gp_y[0] = phi[k - N];
			gp_y[1] = bound_b_o[i];

			d_dx(gp_y, &(phi_y[k]));
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
//				Applying B.C. for mu
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

//		Computing 1st and 2nd order derivatives of mu

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

		mix_E[nt] = 0.0;
		Interface_E[nt] = 0.0;

		double phi_xy;
		double den = 0.0, area = 0.0;

		for (int np = 0; np < NN; np++) {

			phi_xy = phi_x[np] * phi_x[np] + phi_y[np] * phi_y[np];

			mix_E[nt] += (phi[np] - 1.0) * (phi[np] - 1.0) * phi[np] * phi[np];

			Interface_E[nt] += phi_xy;

			den += sqrt(phi_xy);
			area += phi[np];

			laplacian_mu = mu_xx[np] + mu_yy[np];

			phi_new[np] = phi[np] + dt * laplacian_mu;

		}

		mix_E[nt] *= 16.0;
		Interface_E[nt] *= 5.0;

		r_c[nt] = area / den;

		for (int np = 0; np < NN; np++) {

			phi[np] = phi_new[np];

		}


		if (nt != max_nt)
			fprintf(fp1, "%.4f, %.4f, %.4f\n", mix_E[nt], Interface_E[nt], r_c[nt]);
		else
			fprintf(fp1, "%.4f, %.4f, %.4f", mix_E[nt], Interface_E[nt], r_c[nt]);

		if (nt == 2000)
			Output_(N, NN, nt, phi_new);
		else if (nt == 5000)
			Output_(N, NN, nt, phi_new);
		else if (nt == 25000)
			Output_(N, NN, nt, phi_new);
		else if (nt == 100000)
			Output_(N, NN, nt, phi_new);
		else if (nt == 250000)
			Output_(N, NN, nt, phi_new);

	}

	fprintf(fp1, "];\n\n");
	fprintf(fp1, "t = %f:%f:%f;\n\n",dt,dt,t);
	fprintf(fp1, "figure \nloglog(t,recs(:,1),t,recs(:,2),t,recs(:,1)+recs(:,2)) \n");
	fprintf(fp1, "grid on \n\nfigure \nloglog(t,recs(:,3)); \ngrid on");

	fclose(fp1);

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
	delete[] phi_x;
	delete[] mu_xx;

	delete[] phi_yy;
	delete[] phi_y;
	delete[] mu_yy;

	delete[] mix_E;
	delete[] Interface_E;
	delete[] r_c;

	return 0;

}

