/*
CHL for CS211 UCR 2017F
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void dgemm0(double a[], double b[], double c[], int n) {
	clock_t begin, end;
	long diff;
	double time;
	int i, j, k;

	begin = clock();
	for (i = 0; i < n; ++i)
		for (j = 0; j < n; ++j)
			for (k = 0; k < n; ++k)
				c[i*n + j] += a[i*n + k] * b[k*n + j];
	end = clock();
	diff = end - begin;
	time = (double)diff / CLOCKS_PER_SEC;
	printf("Degmm0, n=%d, Time:%.3f Seconds.\n", n, diff, time);
}

void dgemm1(double a[], double b[], double c[], int n) {
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
	printf("Degmm1, n=%d, Time:%.3f Seconds.\n", n, diff, time);
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

/*
CHL for CS211 UCR 2017F
*/
int main(int argc, char* argv[]) {
	
	int task[6] = { 64,128,256,512,1024,2048 };
	int i, j, n;
	double *a, *b, *c0, *c1;
	double maxDiff;
	srand(419);

	for (i = 0; i < 6; ++i) {
		n = task[i];
		a = malloc(n * n * sizeof(double));
		b = malloc(n * n * sizeof(double));
		c0 = malloc(n * n * sizeof(double));
		c1 = malloc(n * n * sizeof(double));
		for (j = 0; j < n*n; ++j) {
			a[j] = (double)((rand() << 15) | rand()) / (double)rand();
			b[j] = (double)((rand() << 15) | rand()) / (double)rand();
			c0[j] = 0;
			c1[j] = 0;
		}

		dgemm0(a, b, c0, n);
		dgemm1(a, b, c1, n);

		maxDiff = CheckMaxDiff(c0, c1, n);

		printf("Maximum difference is: %e\n", maxDiff);

		free(a);
		free(b);
		free(c0);
		free(c1);
	}

	return 0;
}