#define SAY_HI printf("Hi from process %d of total %d processes\n", current_pr, total_pr);
//#define WO_PIVOT_SEARCH_ATALL
#define W_REVERSE_FLOW
#define W_FULL_PIVOT_SEARCH

double f(int, int);
double identity(int, int);
int readMatrixByRows(const char*,double*, int, int, int, int);
int readMatrixByColumns(const char*,double*, int, int, int, int);
int initMatrixByRows(double*, int, int, int, int);
int initMatrixByColumns(double*, int, int, int, int);
int printUpperLeftBlock(double*, int, int);
int blocksMultiply_baad(double *a, int j, int l, 
							double *b, int h, int k, 
							double *out, 
							int block_string_size, int block_size, 
							int short_block_string_size, int small_block_side, 
							int block_side, int total_block_rows);
double MPI_getResidual(double *reinitialized, double *reversed, int matrix_side, int block_side, int total_pr, int current_pr, const int* blocks_order_reversed);
int initParameters(int matrix_side, int block_side, int total_pr, int current_pr, 
	int *total_block_rows, int *total_full_block_rows, 
	int *block_size, int *block_string_size, 
	int *max_block_rows_pp, int *max_rows_pp, int *short_block_string_size,
	int *last_block_row_proc_id, int *last_block_row_in_current_pr,
	int *small_block_row_width, int *small_block_size,
	int *current_pr_full_rows, int *last_block_row_width,
	int *matrix_size_current_pr);
void MPI_printUpperLeftBlock(double *a, int matrix_side, int block_side, int total_pr, int current_pr, const int *blocks_order_reversed);
