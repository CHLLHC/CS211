#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

//HDdiff is from Stackoverflow
struct timespec HDdiff(struct timespec start, struct timespec end)
{
	struct timespec temp;
	if ((end.tv_nsec - start.tv_nsec) < 0) {
		temp.tv_sec = end.tv_sec - start.tv_sec - 1;
		temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
	}
	else {
		temp.tv_sec = end.tv_sec - start.tv_sec;
		temp.tv_nsec = end.tv_nsec - start.tv_nsec;
	}
	return temp;
}

void sijk(double a[], double b[], double c[], int n) {
	struct timespec begin, end, diff;
	double time;
	int i, j, k;

	clock_gettime(CLOCK_MONOTONIC, &begin);
	for (i = 0; i < n; ++i)
		for (j = 0; j < n; ++j) {
			register double r = c[i*n + j];
			for (k = 0; k < n; ++k) {
				r += a[i*n + k] * b[k*n + j];
			}
			c[i*n + j] = r;
		}
	clock_gettime(CLOCK_MONOTONIC, &end);
	diff = HDdiff(begin, end);
	printf("Simple ijk, n=%d, Time:%ld seconds and %ld nanoseconds.\n", n, diff.tv_sec, diff.tv_nsec);
}

void bijk(double a[], double b[], double c[], int n, int B) {
	struct timespec begin, end, diff;
	double time;
	int i, j, k, i1, j1, k1;

	clock_gettime(CLOCK_MONOTONIC, &begin);

	for (i = 0; i < n; i += B)
		for (j = 0; j < n; j += B)
			for (k = 0; k < n; k += B) {
				/* B x B mini matrix multiplications */
				for (i1 = i; i1 < i + B; i1 += 2)
					for (j1 = j; j1 < j + B; j1 += 2) {
						register int in = i1*n, i1n = in + n;
						register double cij = c[in + j1], ci1j = c[i1n + j1],
							cij1 = c[in + (j1 + 1)], ci1j1 = c[i1n + (j1 + 1)];
						for (k1 = k; k1 < k + B; k1 += 2) {
							register double aik = a[in + k1], aik1 = a[in + k1 + 1],
								ai1k = a[i1n + k1], ai1k1 = a[i1n + k1 + 1];
							register int kn = k1*n, k1n = kn + n;
							register double bkj = b[kn + j1], bk1j = b[k1n + j1],
								bkj1 = b[kn + (j1 + 1)], bk1j1 = b[k1n + (j1 + 1)];
							cij = aik * bkj + aik1 * bk1j + cij;
							ci1j = ai1k * bkj + ai1k1 * bk1j + ci1j;
							cij1 = aik * bkj1 + aik1 * bk1j1 + cij1;
							ci1j1 = ai1k * bkj1 + ai1k1 * bk1j1 + ci1j1;
						}
						c[in + j1] = cij;
						c[i1n + j1] = ci1j;
						c[in + (j1 + 1)] = cij1;
						c[i1n + (j1 + 1)] = ci1j1;
					}
			}

	clock_gettime(CLOCK_MONOTONIC, &end);
	diff = HDdiff(begin, end);
	printf("Blocked cache and register ijk, n=%d, B=%d, Time:%ld seconds and %ld nanoseconds.\n", n, B, diff.tv_sec, diff.tv_nsec);
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

	int i, j, n, B = 16;
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
	sijk(a, b, c, n);
	printf("This is the reference result.\n");

	//BLOCKED ALGORITHM
	memset(c1, 0, sizeof(double)*n*n);
	bijk(a, b, c1, n, B);
	maxDiff = CheckMaxDiff(c, c1, n);
	printf("Maximum difference to the reference result is: %e\n", maxDiff);


	free(a);
	free(b);
	free(c);
	free(c1);
	return 0;
}