#include "PhaseField.h"

void d_dxx(double P_gp, double gp2[], double* deriv) {

	*deriv = gp2[1] - 2.0 * P_gp + gp2[0];

	return;
}

