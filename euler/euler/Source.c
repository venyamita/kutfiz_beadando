#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*Nagy P�ter M07ILF
2017.05.10
*/

struct Variable {


	size_t size;
	double *Coords;
};

struct StrMtx {

	size_t width;
	size_t height;
	double* Mtx;
};
void PrintVariable(struct Variable var, FILE * stream) {

	size_t i;
	for (i = 0; i<var.size; i++) {
		fprintf(stream, "%lf ", var.Coords[i]);
	}
	fprintf(stream, "\n");
}

void PrintMatrix(struct StrMtx inputmatrix, char* filename) {
	FILE *outstream;
	outstream = fopen(filename, "w");
	int i, j;
	for (i = 0; i<inputmatrix.height; i++) {
		for (j = 0; j<inputmatrix.width; j++) {
			fprintf(outstream, "%f ", inputmatrix.Mtx[inputmatrix.width*i + j]);
		}
		fprintf(outstream, "\n");
	}
	fclose(outstream);
}

void Gravity(struct Variable in, double t) {
	double GAMMA = 5.9736E24*6.67384E-11;
	double x = in.Coords[0];
	double y = in.Coords[1];
	double vx = in.Coords[2];
	double vy = in.Coords[3];
	in.Coords[0] = vx;
	in.Coords[1] = vy;
	double distrecip = 1. / sqrt(x*x + y*y);
	in.Coords[2] = -GAMMA*x*distrecip*distrecip*distrecip;
	in.Coords[3] = -GAMMA*y*distrecip*distrecip*distrecip;
}

double Energy(struct Variable var) {
	double GAMMA = 5.9736E24*6.67384E-11;
	double m = 7.349E22;
	size_t i;
	double energy = 0.;
	energy = (var.Coords[2] * var.Coords[2] + var.Coords[3] * var.Coords[3])*m*0.5;
	energy = energy + GAMMA*1. / sqrt(var.Coords[0] * var.Coords[0] + var.Coords[1] * var.Coords[1]);
	return energy;
}

void Add(struct Variable A1, struct Variable A2) {
	if (A1.size != A2.size) {
		printf("n ");
		exit(1);
	}
	size_t i;
	for (i = 0; i<A1.size; i++) {
		A1.Coords[i] = A1.Coords[i] + A2.Coords[i];
	}
}

void Multiply(double coeff, struct Variable A) {
	size_t i;
	for (i = 0; i<A.size; i++) {
		A.Coords[i] = coeff*A.Coords[i];
	}
}

void EulerStep(void(*RightHandSide)(struct Variable, double), struct Variable Initial, double t, double dt, struct Variable dummy) {
	size_t i;
	for (i = 0; i < dummy.size; i++)
	{
		dummy.Coords[i] = Initial.Coords[i];
	}
	RightHandSide(dummy, t);
	Multiply(dt, dummy);
	Add(Initial, dummy);
}

void EulerWrapper(void(*RightHandSide)(struct Variable, double), struct Variable Initial, double t0, double tmax, double dt, char* filename) {
	struct Variable dummy;

	dummy.size = Initial.size;
	dummy.Coords = (double*)malloc(dummy.size * sizeof(double));
	double t = t0;
	FILE *out;
	out = fopen(filename, "w");
	FILE *energy;
	energy = fopen("energy_t_euler.txt", "w");
	while (t<tmax) {
		EulerStep(RightHandSide, Initial, t, dt, dummy);
		t = t + dt;
		fprintf(energy, "%f \t %lf\n", t, Energy(Initial));
		PrintVariable(Initial, out);
	}
	fclose(energy);
	fclose(out);
	free(dummy.Coords);

}

int main(int argc, char const *argv[]) {

	struct Variable Initgrav = { 4 };
	Initgrav.Coords = (double *)malloc(4 * sizeof(double));
	Initgrav.Coords[0] = 405500000.;
	Initgrav.Coords[1] = 0.;
	Initgrav.Coords[2] = 0.;
	Initgrav.Coords[3] = -964.;
	double t0 = 0.;
	double period = 3600.*28.*24.;
	double tmax = 50 * period;

	EulerWrapper(Gravity, Initgrav, t0, tmax, 1000., "gravity_moon_euler.txt");
	return 0;
}
