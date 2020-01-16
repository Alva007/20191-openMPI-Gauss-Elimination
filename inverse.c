#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

const int TAG = 100;
const double _zeros[1000000] = {0};
int rank, size;

void _quit (const char *msg) {
	printf("%s\n", msg);
	exit(EXIT_FAILURE);
	MPI_Finalize();
}

void read_matrix_data(const char *path, double **_buff, int *N) {
	FILE *fp;
	fp = fopen(path, "r");
	fscanf(fp, "SHAPE = %d\n", N);
	const int n = *N;
	*_buff = (double*) malloc(n* 2*n *sizeof(double)); // matrix n row and 2*n col (origin matrix append identify matrix based-on Gausse Elimination)
	double *buff = *_buff;
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 2*n; j++) {
			if (j < n) {
				if (j == n-1) fscanf(fp, "%lf\n", &buff[i * 2*n + j]);
				else fscanf(fp, "%lf ", &buff[i * 2*n + j]);
			} else {
				buff[i * 2*n + j] = 0;
				if (j % n == i) buff[i * 2*n + j] = 1;
			}
		}
	}

	fclose(fp);
}

void _toString(double *mtx, int R, int C, int fromR, int fromC) {
	for (int i = fromR; i < R; i++) {
		for (int j = fromC; j < C; j++) {
			printf("%.1f ", mtx[i * C + j]);
		}
		printf("\n");
	}
}

void toString(double *mtx, int R, int C) {
	_toString(mtx, R, C, 0, 0);
}

int main(int argc, char **argv) {
	double t_start, t_end, t_total;

	int code = MPI_Init(&argc, &argv);
	if (code != MPI_SUCCESS) _quit("Error");

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double *matrix;
	int N;

	if (rank == 0) {
		int test;
		if (argc < 2) _quit("Missing param: --test, EX: --test 5");
		else sscanf(argv[2], "%d", &test);

		char path2source[256]; 
		sprintf(path2source, "./data/test_%d", test);
		read_matrix_data(path2source, &matrix, &N);
	} 

	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// number of rows for each rank
	int num_rows = (int) N / size;

	double  *sub_matrix = (double*) malloc(num_rows * 2*N * sizeof(double));
	MPI_Scatter(matrix, num_rows * 2*N, MPI_DOUBLE, sub_matrix, num_rows * 2*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// row var be used for elimination in each rank
	double *row = (double*) malloc(2*N * sizeof(double));
	if (rank == 0) {
		t_start = MPI_Wtime();
	}

	// Gauss Elimination phase
	int pivot;
	int scale;
	int column;
	int start_row;

	// current rank receivers 
	start_row = rank * num_rows;
	for (int i = 0; i < start_row; i++) {
		// Wait for the preceeding ranks (0 to rank-1) to forward this rank a row
		// Each rank,  call MPI_Bcast num_rows times for 
		MPI_Bcast(row, 2*N, MPI_DOUBLE, i / num_rows, MPI_COMM_WORLD);

		// Eliminate from this element from all rows mapped to this rank
		for (int j = 0; j < num_rows; j++) {
			scale = sub_matrix[j * 2*N + i];
			for (int k = i + 1; k < 2*N; k++) {
				sub_matrix[j * 2*N + k] -= scale * row[k];
			}

			sub_matrix[j * 2*N + i] = 0;
		}
	}

	// current rank senders
	for (int i = 0; i < num_rows; i++) {
		// column index of pivot
		column = rank * num_rows + i;
		pivot = sub_matrix[i * 2*N + column];
		// check pivot
		// if (!pivot) {
			// Normalize every other element in this row to the pivot
			for (int j = column + 1; j < 2*N; j++) {
				sub_matrix[i * 2*N + j] /= pivot;
			}

			sub_matrix[i * 2*N + column] = 1;
			
			// Fill row to be sent
			memcpy(row, &sub_matrix[i * 2*N], 2*N * sizeof(double));

			// Bcast the rest of the rows for this rank
			MPI_Bcast(row, 2*N, MPI_DOUBLE, rank, MPI_COMM_WORLD);

			for (int j = i+1; j < num_rows; j++) {
				scale = sub_matrix[j * 2*N + column];

				for (int k = column + 1; k < 2*N; k++) {
					sub_matrix[j * 2*N + k] -= scale * row[k];
				}

				sub_matrix[j * 2*N + column] = 0;
			}

		// } else { // pivot == 0
		//	memcpy(row, &_zeros[0], N * sizeof(double));
		// }

	}


	/*
	for (int i = rank * num_rows + 1; i < 2*N; i++) {
		MPI_Bcast(row, 2*N, MPI_DOUBLE, i / num_rows, MPI_COMM_WORLD);
	}
	*/

	MPI_Barrier(MPI_COMM_WORLD);

	// Back elimination phase

	if (rank == 0) {
		t_end = MPI_Wtime();
		t_total = t_end - t_start;
	}

	MPI_Gather(sub_matrix, num_rows * 2*N, MPI_DOUBLE, matrix, num_rows * 2*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		printf("\nDONE: \n");
		toString(matrix, N, 2*N);

		printf("Totol time: %f\n", t_total);
		free(matrix);
	}

	free(sub_matrix);
	free(row);


	MPI_Finalize();

	return 0;
}
