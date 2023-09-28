#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

// 1.) Create Random Numbers based on command line prompt * check *
// 2.) Create non zero elements * check *
// 3.) Do the formulations

// For passing for threads
typedef struct PAR{
    int vectorStart;
    int vectorEnd;
}parameter;


// Global Variables for easy access
double** matrix;
double* vector;
double* correct_vector;

int threadSize;
int n;

int* proper_distribute;

// Functions to be used
// n is equal to the size given in prompt
double** createMatrix(int n){
	double** matrix = (double**)malloc(n* sizeof(double*));
	for(int i = 0; i < n; i++) matrix[i] = malloc(n*sizeof(double));
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			// For random variables
			matrix[i][j] = (double) ((rand() % 100) + 1);
		}
	}
	return(matrix);
}

// Same as createMatrix but single array instead of 2d
double* createVector(int n){
	double* vector = (double*)malloc(n*sizeof(double));
	for(int i = 0; i < n; i++){
		vector[i] = (double) ((rand() % 100) + 1); ; 
	}
	
	return(vector);
}

double** testMatrix(int n){
	double** matrix = (double**)malloc(n* sizeof(double*));
	for(int i = 0; i < n; i++) matrix[i] = malloc(n*sizeof(double));
	
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			// For random variables
			matrix[i][j] = (i + 1);	
		}
	}
	return(matrix);
}

double* testVector(int n){
	double* vector = (double*)malloc(n*sizeof(double));
	for(int i = 0; i < n; i++){
		vector[i] = (i + 1); 
	}
	
	return(vector);
}

void printMatrix(double** matrix, int n){
	printf("Test Matrix\n");
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			printf("%f ", matrix[i][j]);
		}
		printf("\n");
	}
}

void printVector(double* vector, int n){
	printf("Test Vector\n");
	for(int i = 0; i < n; i++){
		printf("%f ", vector[i]);
	}
	printf("\n");
}

// Pearson Coefficient Formulas
double xj(double** matrix, double* vector, int j, int m){
	double sum = 0;
	for(int i = 0; i < m; i++){
		sum += matrix[i][j];
	}
	return sum;
}

double x2j(double** matrix, double* vector, int j, int m){
	double sum = 0;
	for(int i = 0; i < m; i++){
		sum += (matrix[i][j] * matrix[i][j]);
	}
	return sum;
}


double y(double** matrix, double* vector, int j, int m){
	double sum = 0;
	for(int i = 0; i < m; i++){
		sum += vector[i];
	}
	return sum;
}

double y2(double** matrix, double* vector, int j, int m){
	double sum = 0;
	for(int i = 0; i < m; i++){
		sum += (vector[i] * vector[i]);
	}
	return sum;
}

double x2y(double** matrix, double* vector, int j, int m){
	double sum = 0;
	for(int i = 0; i < m; i++){
		sum += (matrix[i][j] * vector[i]);
	}
	return sum;
}	


void pearson_cor(double** matrix, double* vector, int m){
	// for i:=1 to n
	for (int i = 0; i < m; i++){
		// Convert answers to float to be able to have final answer as float
		double x = xj(matrix, vector, i,m);
		double x2 = x2j(matrix, vector, i, m);
		double yj = y(matrix, vector, i, m);
		double y2j = y2(matrix, vector, i, m);
		double xy = x2y(matrix, vector, i, m);
		vector[i] = ((m * xy) - (x * yj)) / sqrt( ((m * x2) - (x*x))  * ((m * y2j) - (yj*yj)) );
	}
}

void pearson_corr(double** matrix, double* vector, int start, int m, int end){
	// for i:=1 to n
	for (int i = start; i < m; i++){
		// Convert answers to float to be able to have final answer as float
		double x = xj(matrix, vector, i,m);
		double x2 = x2j(matrix, vector, i, m);
		double yj = y(matrix, vector, i, m);
		double y2j = y2(matrix, vector, i, m);
		double xy = x2y(matrix, vector, i, m);
		correct_vector[i] = ((m * xy) - (x * yj)) / sqrt( ((m * x2) - (x*x))  * ((m * y2j) - (yj*yj)) );
		// printf("answer = %f\n", vector[i]);
	}
}

int* distribute_arrays (int n, int t){ 
	int quotient = n / t;
	int remainder = n % t;

	int* array = (int*)malloc(t*sizeof(int));

	for (int i = 0; i < t; i++){
		array[i] = quotient;
		if(i < remainder) array[i] += 1;
	}

	return array;

}

double** transpose_matrix(double** matrix, int n) {
  // Create a new matrix to store the transposed matrix.
  double** transposed_matrix = malloc(n * sizeof(double*));
  for (int i = 0; i < n; i++) {
    transposed_matrix[i] = malloc(n * sizeof(double));
  }

  // Iterate over the rows and columns of the original matrix.
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      // Copy the element at the current row and column to the corresponding row and column in the transposed matrix.
      transposed_matrix[j][i] = matrix[i][j];
    }
  }

  return transposed_matrix;
}

void print_array(int* array, int quotient){
	for(int i = 0; i < quotient; i++){
		printf("%d ", array[i]);
	}
	printf("\n");
}

void* thread_pearson(void* arg) {
    parameter* temp = (parameter*) arg;
	int start = temp->vectorStart;
	int end = temp->vectorEnd;
	pearson_corr(matrix, vector, start, n, end);
    pthread_exit(NULL);
}

// main program
int main(int argc, char* argv[]){
	// struct timeval start_time, end_time;
	// time_t t;
	uint64_t diff;
	struct timespec start,end;

    // Remember to dumbproof this later |||
	// number which nxn is based off
	n = strtol(argv[1], NULL, 10);
    // Equals to the number of threads
    threadSize = strtol(argv[2], NULL, 10);

	if(threadSize >= n){
		printf("Thread Size should never be greater than or equal to N");
		return 0;
	}

    // creation of number of threads
    pthread_t threads[threadSize];
	parameter params[threadSize];

	// for random numbers per run
	srand(time(NULL));

    // divide submatrices by n x n/t
    int submatrices =  n * n/threadSize;
	int prevStart = 0;
	int prevEnd = 0;

	correct_vector = (double*)malloc(n*sizeof(double));

	// For even distribution to threads
	proper_distribute = distribute_arrays(n, threadSize);

	for(int i = 0; i < threadSize; i++){
		prevStart = prevEnd;
		prevEnd = prevEnd + (proper_distribute[i]);
		if(prevStart < 0) prevStart = 0;
		params[i].vectorEnd = prevEnd;
		params[i].vectorStart = prevStart;
		// printf("start = %d\n", params[i].vectorStart);
		// printf("end = %d\n", params[i].vectorEnd);
		// printf("\n");
	}

	// print_array(proper_distribute, threadSize);

	
	// Generate the matrix and vector
	matrix = createMatrix(n);
	vector = createVector(n);
	
	// For debugging
	// matrix = testMatrix(5);
	// vector = testVector(5);
	
	// test the matrix and vector if generating
	// printMatrix(matrix, n);
	// printVector(vector, n);

	// clock_t begin = clock();
	clock_gettime(CLOCK_MONOTONIC, &start);	
	
	// pearson_cor(matrix, vector, n);

    for (int i = 0; i < threadSize; i++) {
        pthread_create(&threads[i], NULL, thread_pearson, (void*) &params[i]);
    }

    // For synchronization
    for (int i = 0; i < threadSize; i++) {
        pthread_join(threads[i], NULL);
    }
    
	// clock_t end = clock();
	clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    // printf("Main Code finished\n");
	// printFloat(correct_vector, n);

	double execution_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
	printf("Execution time: %f seconds\n", execution_time);

	// double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	
	// printf("\nTime Spent == %f\n", time_spent);
	// printf("end\n");

	// printVector(vector, n);
	// printVector(correct_vector, n);
		
	return 0;
}