#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void sijk(double a[], double b[], double c[], int n) {
	clock_t begin, end, diff;
	double time;
	int i, j, k;

	begin = clock();
	for (i = 0; i < n; ++i)
		for (j = 0; j < n; ++j) {
			register double r = c[i*n + j];
			for (k = 0; k < n; ++k) {
				r += a[i*n + k] * b[k*n + j];
			}
			c[i*n + j] = r;
		}
	end = clock();
	diff = end - begin;
	time = (double)diff / CLOCKS_PER_SEC;
	printf("Simple ijk, n=%d, Time:%.3f Seconds.\n", n, diff, time);
}

void dgemm2(double a[], double b[], double c[], int n) {
	clock_t begin, end, diff;
	double time;
	int i, j, k;

	begin = clock();
	for (i = 0; i < n; i += 2)
		for (j = 0; j < n; j += 2) {
			register double cij = c[i*n + j], ci1j = c[(i + 1)*n + j],
				cij1 = c[i*n + (j + 1)], ci1j1 = c[(i + 1)*n + (j + 1)];
			for (k = 0; k < n; k += 2) {
				register double aik = a[i*n + k], aik1 = a[i*n + k + 1],
					ai1k = a[(i + 1)*n + k], ai1k1 = a[(i + 1)*n + k + 1];
				register double bkj = b[k*n + j], bk1j = b[(k + 1)*n + j],
					bkj1 = b[k*n + (j + 1)], bk1j1 = b[(k + 1)*n + (j + 1)];
				cij = aik * bkj + aik1 * bk1j + cij;
				ci1j = ai1k * bkj + ai1k1 * bk1j + ci1j;
				cij1 = aik * bkj1 + aik1 * bk1j1 + cij1;
				ci1j1 = ai1k * bkj1 + ai1k1 * bk1j1 + ci1j1;
			}
			c[i*n + j] = cij;
			c[(i + 1)*n + j] = ci1j;
			c[i*n + (j + 1)] = cij1;
			c[(i + 1)*n + (j + 1)] = ci1j1;
		}
	end = clock();
	diff = end - begin;
	time = (double)diff / CLOCKS_PER_SEC;
	printf("Degmm2, n=%d, Time:%.3f Seconds.\n", n, diff, time);
}

void dgemm3(double a[], double b[], double c[], int n) {
	clock_t begin, end, diff;
	double time;
	int i, j, k;

	begin = clock();
	for (i = 0; i < n; i += 2)
		for (j = 0; j < n; j += 2) {
			register int in = i*n, i1n = in + n;
			register double cij = c[in + j], ci1j = c[i1n + j],
				cij1 = c[in + (j + 1)], ci1j1 = c[i1n + (j + 1)];
			for (k = 0; k < n; k += 2) {
				register double aik = a[in + k], aik1 = a[in + k + 1],
					ai1k = a[i1n + k], ai1k1 = a[i1n + k + 1];
				register int kn = k*n, k1n = kn + n;
				register double bkj = b[kn + j], bk1j = b[k1n + j],
					bkj1 = b[kn + (j + 1)], bk1j1 = b[k1n + (j + 1)];
				cij = aik * bkj + aik1 * bk1j + cij;
				ci1j = ai1k * bkj + ai1k1 * bk1j + ci1j;
				cij1 = aik * bkj1 + aik1 * bk1j1 + cij1;
				ci1j1 = ai1k * bkj1 + ai1k1 * bk1j1 + ci1j1;
			}
			c[in + j] = cij;
			c[i1n + j] = ci1j;
			c[in + (j + 1)] = cij1;
			c[i1n + (j + 1)] = ci1j1;
		}
	end = clock();
	diff = end - begin;
	time = (double)diff / CLOCKS_PER_SEC;
	printf("Degmm3, n=%d, Time:%.3f Seconds.\n", n, diff, time);
}

double CheckMaxDiff(double a[], double b[], int n) {
	double ret = 0;
	int i;

	for (i = 0; i < n*n; ++i) {
		if (a[i] - b[i] > ret) {
			ret = a[i] - b[i];
		}
		if (b[i] - a[i] > ret) {
			ret = b[i] - a[i];
		}
	}

	return ret;
}


int main(int argc, char* argv[]) {

	int blcokSize[6] = { 2,4,8,16,32,64 };
	int i, j, n, B;
	double *a, *b, *c, *c1;
	double maxDiff;
	srand(419);
	n = 2048;
	a = malloc(n * n * sizeof(double));
	b = malloc(n * n * sizeof(double));
	c = malloc(n * n * sizeof(double));
	c1 = malloc(n * n * sizeof(double));
	for (j = 0; j < n*n; ++j) {
		a[j] = (double)((rand() << 15) | rand()) / (double)rand();
		b[j] = (double)((rand() << 15) | rand()) / (double)rand();
	}
	memset(c, 0, sizeof(double)*n*n);
	dgemm1(a, b, c, n);


	for (i = 0; i < 6; ++i) {
		B = blcokSize[i];
	}
	free(a);
	free(b);
	return 0;
}