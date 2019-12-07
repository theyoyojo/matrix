#include <stdio.h>
#include <string.h>
#include <assert.h>


#define ROWMAX 64UL
#define COLMAX 64UL

#define str(__arg) #__arg

#define ROWNUMBER_MAXWIDTH str(2)
#define COLNUMBER_MAXWIDTH str(2)

typedef struct matrix {
	double numbers[ROWMAX][COLMAX] ;
	size_t rows ;
	size_t cols ;
} matrix_t ;

void matrix_init(matrix_t * matrix, size_t cols) {
	matrix->cols = cols ;
	matrix->rows = 0 ;
}

int matrix_add_row(matrix_t * matrix, double * row, size_t row_size) {
	assert(matrix) ;
	if (row_size != matrix->cols) {
		fprintf(stderr, "Error: new row is too long (%zu != %zu)\n", matrix->rows, row_size) ;
		return  -1 ;
	}
	else if (matrix->rows + 1 >= ROWMAX) {
		fprintf(stderr, "Error: row max exceeded (%zu)\n", ROWMAX) ;
		return -2 ;
	}

	memcpy(matrix->numbers[matrix->rows++], row, matrix->cols * sizeof(double)) ;
	return 0 ;
}

#define ROW_DECOR(__number) do { printf("[%" ROWNUMBER_MAXWIDTH "zu]", __number) ; } while (0)
#define COL_DECOR(__number) do { printf("   [%" COLNUMBER_MAXWIDTH "zu]   ", __number) ; } while (0)

void matrix_print(matrix_t * matrix) {
	assert(matrix) ;

	printf("[%" ROWNUMBER_MAXWIDTH "zux%-" COLNUMBER_MAXWIDTH "zu]", matrix->rows, matrix->cols) ;
	for (size_t i = 0; i < matrix->cols; ++i) {
		COL_DECOR(i) ;
	}
	printf("\n") ;

	for (size_t i = 0; i < matrix->rows; ++i) {
		ROW_DECOR(i) ;
		for (size_t j = 0; j < matrix->cols; ++j) {
			printf("%8.2lg |", matrix->numbers[i][j]) ;
		}
		printf("\n") ;
	}
}

void matrix_scan(matrix_t * matrix) {
	assert(matrix) ;
	int ex_rows, ex_cols ;
	double buf[ROWMAX] ;
	scanf("%d", &ex_rows) ;
	scanf("%d", &ex_cols) ;

	matrix_init(matrix, ex_cols) ;

	for (int i = 0; i < ex_rows; ++i) {
		for (int j = 0; j < ex_cols; ++j) {
			scanf("%lf", &buf[j]) ;
		}
		matrix_add_row(matrix, buf, ex_cols) ;
	}
}

int matrix_multiply(matrix_t * A, matrix_t * B, matrix_t * C) {
	assert(A), assert(B) ;

	if (A->cols != B->rows) {
		fprintf(stderr, "Error: number of cols in first matrix must match number of rows in second\n") ;
		return -1 ;
	}

	C->cols = A->rows ;
	C->rows = B->cols ;
	memset(C->numbers, 0, C->rows * C->cols * sizeof(double)) ;

	for (size_t i = 0; i < C->rows; ++i) {
		for (size_t j = 0; j < C->cols; ++j) {
			for (size_t k = 0; k < C->rows; k++) {
				C->numbers[i][j] += A->numbers[i][k] * B->numbers[k][j] ;
			}
		}
	}

	return 0 ;

}

int matrix_augment(matrix_t * A, matrix_t * B, matrix_t * C) {
	return 0 ;
}

int matrix_row_add(matrix_t * A, size_t row_start, size_t row_target, double weight, int trace) ;
inline int matrix_row_add(matrix_t * A, size_t row_start, size_t row_target, double weight, int trace) {
	if (trace) {
		printf(">> R%zu += %lg * R%zu\n", row_target, weight, row_start) ;
		printf("(\?\?) +") ;
	}
	for (size_t i = 0; i < A->cols; ++i) {
		A->numbers[row_target][i] += A->numbers[row_start][i] * weight ;
		if (trace) {
			printf("%8.2lg  ", A->numbers[row_start][i] * weight) ;
		}

	}
	if (trace) {
		printf("to row %zu\n", row_target) ;
	}

	return 0 ;
}

int matrix_row_interchange(matrix_t * A, size_t row1, size_t row2, int trace) {

	if (row1 >= A->rows || row2 >= A->rows) {
		fprintf(stderr, "Error: invalid row(s): cannot interchange\n") ;
		return -1 ;
	}
	double tmp ;
	if (trace) {
		printf(">> R%zu <---> R%zu\n", row1, row2) ;
	}
	if (row1 != row2) {
		for (size_t i = 0; i < A->cols; ++i) {
			if (A->numbers[row1][i] != A->numbers[row2][i]) {
				tmp = A->numbers[row2][i] ;
				A->numbers[row2][i] = A->numbers[row1][i] ;
				A->numbers[row1][i] = tmp ;
			}
		}
	}
	 
	if (trace) {
		matrix_print(A) ;
	}
	
	return 0 ;
}

int matrix_row_scale(matrix_t * A, size_t row, double factor, int trace) {
	if (row >= A->rows) {
		fprintf(stderr, "Error: invalid rows: cannot interchange\n") ;
		return -1 ;
	}

	if (trace) {
		printf(">> R%zu *= %8.2lg\n", row, factor) ;
	}

	for (size_t i = 0; i < A->cols; ++i) {
		A->numbers[row][i] *= factor ;
	}

	if (trace) {
		matrix_print(A) ;
	}

	return 0 ;
}

void matrix_to_echelon(matrix_t * A, int trace) {
	assert(A) ;
	int offset ;

	for  (size_t i = 0; i < A->cols; ++i) {
		offset = 0 ;
		while (A->numbers[i + offset][i] == 0) {
			
			++offset ;

			if (i + offset >= A->rows) {
				goto breakcontinue ;
			}
		}
		if (i != i + offset) {
			matrix_row_interchange(A, i, i + offset, trace) ;
		}
		for (size_t j = i + offset + 1; j < A->rows; ++j) {
			matrix_row_add(A, i, j, -1 * A->numbers[j][i]/A->numbers[i][i], trace) ;
			if (trace) {
				matrix_print(A) ;
			}
		}
breakcontinue:
		(void)0 ;
	}
	
}

long matrix_row_find_pivot(matrix_t * A, size_t row) {
	assert(A) ;
	if (row >= A->rows) {
		return -2 ;
	}

	for(size_t i = 0; i < A->cols; ++i) {
		if (A->numbers[row][i] != 0) {
			return i ;
		}
	}

	return -1 ;
}

void matrix_to_echelon_reduced(matrix_t * A, int trace) {
	assert(A) ;

	long pivot_col ;
	matrix_to_echelon(A, trace) ;

	for (int i = A->rows - 1; i > 0; --i) {
		if ((pivot_col = matrix_row_find_pivot(A, i)) < 0) {
			continue ;
		}
		else {
			for (int j = i - 1; j >= 0; --j) {
				if (A->numbers[j][pivot_col]) {
					matrix_row_add(A, i, j, -1 * A->numbers[j][pivot_col]/A->numbers[i][pivot_col], trace) ;
					if (trace) {
						matrix_print(A) ;
					}
				}
			}
		}
	}

	for (size_t i = 0; i < A->rows; ++i) {
		if ((pivot_col = matrix_row_find_pivot(A, i)) >= 0) {
			matrix_row_scale(A, i, 1 / A->numbers[i][pivot_col], trace) ;
			if (trace) {
				matrix_print(A) ;
			}
		}
	}

}

void matrix_copy(matrix_t * src, matrix_t * dest) {
	assert(src), assert(dest) ;
	memcpy(dest, src, sizeof(matrix_t)) ;
}

int matrix_row_zeros(matrix_t * A, size_t row) {
	if (row >= A->rows) {
		fprintf(stderr, "Error: invalid row: cannot count zeros\n") ;
		return -1 ;
	}

	int count ;

	count = 0 ;
	for (size_t i = 0; i < A->cols; ++i) {
		if (!A->numbers[row][i]) {
			++count ;
		}
	}

	return count ;

}

int matrix_col_zeros(matrix_t * A, size_t col) {
	if (col >= A->cols) {
		fprintf(stderr, "Error: invalid col: cannot count zeros\n") ;
		return -1 ;
	}

	int count ;

	count = 0 ;
	for (size_t i = 0; i < A->cols; ++i) {
		if (!A->numbers[i][col]) {
			++count ;
		}
	}

	return count ;

}

int main(int argc, char * argv[]) {
	matrix_t A/*, B, C */ ;

	matrix_scan(&A) ;
	/* matrix_scan(&B) ; */

	matrix_print(&A) ;
	/* matrix_print(&B) ; */

	/* matrix_multiply(&A, &B, &C) ; */

	/* matrix_print(&C) ; */

	matrix_to_echelon_reduced(&A, 1) ;
	matrix_print(&A) ;

	return 0 ;
}
