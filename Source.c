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

void sjik(double a[], double b[], double c[], int n) {
	struct timespec begin, end, diff;
	double time;
	int i, j, k;

	clock_gettime(CLOCK_MONOTONIC, &begin);
	for (j = 0; j < n; ++j)
		for (i = 0; i < n; ++i) {
			register double r = c[i*n + j];
			for (k = 0; k < n; ++k) {
				r += a[i*n + k] * b[k*n + j];
			}
			c[i*n + j] = r;
		}
	clock_gettime(CLOCK_MONOTONIC, &end);
	diff = HDdiff(begin, end);
	printf("Simple jik, n=%d, Time:%ld seconds and %ld nanoseconds.\n", n, diff.tv_sec, diff.tv_nsec);
}

void sikj(double a[], double b[], double c[], int n) {
	struct timespec begin, end, diff;
	double time;
	int i, j, k;

	clock_gettime(CLOCK_MONOTONIC, &begin);
	for (i = 0; i < n; ++i)
		for (k = 0; k < n; ++k) {
			register double r = a[i*n + k];
			for (j = 0; j < n; ++j) {
				c[i*n + j] += r * b[k*n + j];
			}
		}
	clock_gettime(CLOCK_MONOTONIC, &end);
	diff = HDdiff(begin, end);
	printf("Simple ikj, n=%d, Time:%ld seconds and %ld nanoseconds.\n", n, diff.tv_sec, diff.tv_nsec);
}

void skij(double a[], double b[], double c[], int n) {
	struct timespec begin, end, diff;
	double time;
	int i, j, k;

	clock_gettime(CLOCK_MONOTONIC, &begin);
	for (k = 0; k < n; ++k)
		for (i = 0; i < n; ++i) {
			register double r = a[i*n + k];
			for (j = 0; j < n; ++j) {
				c[i*n + j] += r * b[k*n + j];
			}
		}
	clock_gettime(CLOCK_MONOTONIC, &end);
	diff = HDdiff(begin, end);
	printf("Simple kij, n=%d, Time:%ld seconds and %ld nanoseconds.\n", n, diff.tv_sec, diff.tv_nsec);
}

void sjki(double a[], double b[], double c[], int n) {
	struct timespec begin, end, diff;
	double time;
	int i, j, k;

	clock_gettime(CLOCK_MONOTONIC, &begin);
	for (j = 0; j < n; ++j)
		for (k = 0; k < n; ++k) {
			register double r = b[k*n + j];
			for (i = 0; i < n; ++i) {
				c[i*n + j] += a[i*n + k] * r;
			}
		}
	clock_gettime(CLOCK_MONOTONIC, &end);
	diff = HDdiff(begin, end);
	printf("Simple jki, n=%d, Time:%ld seconds and %ld nanoseconds.\n", n, diff.tv_sec, diff.tv_nsec);
}

void skji(double a[], double b[], double c[], int n) {
	struct timespec begin, end, diff;
	double time;
	int i, j, k;

	clock_gettime(CLOCK_MONOTONIC, &begin);
	for (k = 0; k < n; ++k)
		for (j = 0; j < n; ++j) {
			register double r = b[k*n + j];
			for (i = 0; i < n; ++i) {
				c[i*n + j] += a[i*n + k] * r;
			}
		}
	clock_gettime(CLOCK_MONOTONIC, &end);
	diff = HDdiff(begin, end);
	printf("Simple kji, n=%d, Time:%ld seconds and %ld nanoseconds.\n", n, diff.tv_sec, diff.tv_nsec);
}

void bijk(double a[], double b[], double c[], int n, int B) {
	struct timespec begin, end, diff;
	double time;
	int i, j, k, i1, j1, k1;

	clock_gettime(CLOCK_MONOTONIC, &begin);

	for (i = 0; i < n; i += B)
		for (j = 0; j < n; j += B)
			for (k = 0; k < n; k += B)
				/* B x B mini matrix multiplications */
				for (i1 = i; i1 < i + B; i1++)
					for (j1 = j; j1 < j + B; j1++) {
						register double r = c[i1*n + j1];
						for (k1 = k; k1 < k + B; k1++)
							r += a[i1*n + k1] * b[k1*n + j1];
						c[i1*n + j1] = r;
					}

	clock_gettime(CLOCK_MONOTONIC, &end);
	diff = HDdiff(begin, end);
	printf("Blocked ijk, n=%d, B=%d, Time:%ld seconds and %ld nanoseconds.\n", n, B, diff.tv_sec, diff.tv_nsec);
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
	n = 512;
	a = malloc(n * n * sizeof(double));
	b = malloc(n * n * sizeof(double));
	c = malloc(n * n * sizeof(double));
	c1 = malloc(n * n * sizeof(double));
	for (j = 0; j < n*n; ++j) {
		a[j] = (double)((rand() << 15) | rand()) / (double)rand();
		b[j] = (double)((rand() << 15) | rand()) / (double)rand();
	}
	/*
	memset(c, 0, sizeof(double)*n*n);
	sijk(a, b, c, n);
	printf("This is the reference result.\n");
	memset(c1, 0, sizeof(double)*n*n);
	sikj(a, b, c1, n);
	CheckMaxDiff(c, c1, n);
	memset(c1, 0, sizeof(double)*n*n);
	sjik(a, b, c1, n);
	CheckMaxDiff(c, c1, n);
	memset(c1, 0, sizeof(double)*n*n);
	sjki(a, b, c1, n);
	CheckMaxDiff(c, c1, n);
	memset(c1, 0, sizeof(double)*n*n);
	skij(a, b, c1, n);
	CheckMaxDiff(c, c1, n);
	memset(c1, 0, sizeof(double)*n*n);
	skji(a, b, c1, n);
	CheckMaxDiff(c, c1, n);
	*/
	for (i = 0; i < 6; ++i) {
		B = blcokSize[i];
		memset(c1, 0, sizeof(double)*n*n);
		bijk(a, b, c1, n, B);
	}

	free(a);
	free(b);
	free(c);
	free(c1);
	return 0;
}