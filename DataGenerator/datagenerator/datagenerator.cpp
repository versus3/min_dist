// generates dataset (points) to be used for testing in Min-dist Location Selection problem
//
// Yuan (Andy) Xue 
// Last Modified: 02/10/2011

#include "stdafx.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define N_C 100000		// number of clients
#define N_P 5000		// nubmer of potential locations
#define N_F 5000		// number of existing facilities

// boundary of the square data space
#define MIN 0
#define MAX 1000

// Output filename
#define FILE_N "../tN.txt"  // store n_c, n_p, and n_f
#define FILE_C "../tC.txt"	// store C
#define FILE_P "../tP.txt"	// store P
#define FILE_F "../tF.txt"	// store F

float randRange(int min, int max);


int main(int argc, char* argv[])
{
	int n_c = N_C, n_p = N_P, n_f = N_F, min = MIN, max = MAX;

	// start the random number generator
	srand(0);

	// read input. use default values if they're not entered
	if (argc == 6) {
		sscanf(argv[1], "%d", &n_c);
		sscanf(argv[2], "%d", &n_p);
		sscanf(argv[3], "%d", &n_f);
		sscanf(argv[4], "%d", &min);
		sscanf(argv[5], "%d", &max);
	} else if (argc == 4) {
		sscanf(argv[1], "%d", &n_c);
		sscanf(argv[2], "%d", &n_p);
		sscanf(argv[3], "%d", &n_f);
	}

	fprintf(stderr,
		"n_c=%d, n_p=%d, n_f=%d\nboundary: min=%d, max=%d\n", n_c,
		n_p, n_f, min, max);

    /* OUTPUT: output three sets of data into three files for c, p, and f
     * The format is "%d %.4f %.4f %.4f %.4f".
     * (index, range in x-axis [x, x], range in y-axis [y y])
     * The format is compatible with the RTree implementation
     */
	int i;
    double x, y;

    // output n_c, n_p, n_f
    FILE *file = fopen(FILE_N, "w");
	fprintf(file, "%d %d %d\n", n_c, n_p, n_f);
    
    // output clients
    file = fopen(FILE_C, "w");
	for (i = 0; i < n_c; i++) {
        x = randRange(min, max); y = randRange(min, max);
		fprintf(file, "%d %.4f %.4f %.4f %.4f\n", i+1, x, x, y, y);
	}

    // output candidate locations
    file = fopen(FILE_P, "w");
	for (i = 0; i < n_p; i++) {
        x = randRange(min, max); y = randRange(min, max);
		fprintf(file, "%d %.4f %.4f %.4f %.4f\n", i+1, x, x, y, y);
	}

    // output facilities
    file = fopen(FILE_F, "w");
	for (i = 0; i < n_f; i++) {
        x = randRange(min, max); y = randRange(min, max);
		fprintf(file, "%d %.4f %.4f %.4f %.4f\n", i+1, x, x, y, y);
	}

	return EXIT_SUCCESS;
}

// generate a random double number in a specified range
float randRange(int min, int max) {
	if (min > max) {
		fprintf(stderr, "randRange() failed\n");
		exit(EXIT_FAILURE);
	} else {
		return rand() / ((float) (RAND_MAX) + 1) * (max - min) + min;
	}
}

