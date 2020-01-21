#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const int STEP = 5;
int n = 1000, nTests = 4;
unsigned seed = 10;

void makeTest(int n, int nTest);

int main(int argc, char **argv) {
	if (argc > 1) {
		if (!strcmp(argv[1], "--nTest")) sscanf(argv[2], "%d", &nTests);
		if (!strcmp(argv[3], "--dim")) sscanf(argv[4], "%d", &n);
		if (!strcmp(argv[5], "--seed")) sscanf(argv[6], "%u", &seed);
	}

	srand(seed);

	int bcd = 120;
	int N[] = {bcd*6, bcd*10, bcd*40};
	nTests = sizeof(N) / sizeof(int);
	for (int i = 0; i < nTests; i++) {
		makeTest(N[i], i);
	}
}

void makeTest(int n, int testID) {
	FILE *fp;
	char filename[256];
	sprintf(filename, "./data/test_%d", n);
	fp = fopen(filename, "w");

	fprintf(fp, "SHAPE = %d\n", n);
	int random_number; 
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			random_number = (rand() % n) % 1000 + 1;
			if (j == n-1) fprintf(fp, "%d\n", random_number);
			else fprintf(fp, "%d ", random_number);
		}
	}

	fclose(fp);
}
