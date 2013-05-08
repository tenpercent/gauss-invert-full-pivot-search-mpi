#include <cmath>
#include <cstdio>
#include <mpi.h>
#include "tools.h"
#include "io.h"
double f (int i,int j){
	return fabs(i-j);
	//return (i>j)?i:j;
}
double identity (int i, int j){
    if (i==j) return 1.;
    else return 0.;
}
int readMatrixByRows(const char* name, double *a, int matrix_side, int block_side, int total_pr, int current_pr){
   	FILE *fp = 0;
   	int ret=0, loc=0;
   	int i, j, count=0;
	double trash;
	int block_number = 0;
	int block_number_proc_id = 0;
	int block_number_per_process = 0;
	int total_block_rows = (matrix_side + block_side - 1)/block_side;
	int current_pr_block_rows = (total_block_rows + total_pr - 1)/total_pr;
	int block_string_size = matrix_side * block_side;
	int matrix_size_current_pr = current_pr_block_rows *block_string_size;

   	fp = fopen(name, "r");
   	if (!fp){
   		loc = 1;
   	}
   	MPI_Allreduce(&loc, &ret, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   	if (ret){
       	return -1;
   	}
   	loc=0;
	for(i=0; i<matrix_side; i++){
		block_number = i/block_side;
		block_number_proc_id = block_number % total_pr;
		block_number_per_process = block_number / total_pr;
		if(block_number_proc_id==current_pr){
			for(j=0; j<matrix_side; j++){
				loc+=fscanf(fp, "%lf", a + block_number_per_process*block_string_size + (i%block_side)*matrix_side + j);
				count++;
			}
		}
		else{
			for(j=0; j<matrix_side; j++){
				loc+=fscanf(fp, "%lf", &trash);
			}
		}
	}
	if(loc<matrix_side*matrix_side){
		loc=1;
	}
	else{
		loc=0;
	}
	MPI_Allreduce(&loc, &ret, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   	if (ret){
   		if(current_pr==0){
   			printf("Failed to read whole matrix\n");
   		}
       	return -1;
   	}
	fclose(fp);
	for(j=count;j<matrix_size_current_pr;j++){
		a[j]=0.;
	}
    return 0;
}
int readMatrixByColumns(const char* name, double *a, int matrix_side, int block_side, int total_pr, int current_pr){
    FILE *fp = 0;
   	int ret=0, loc=0;
	int i, j;
	int count = 0;
	double trash;
	int block_number = 0;
	int block_number_proc_id = 0;
	int block_number_per_process = 0;

	int total_block_rows, total_full_block_rows, block_size, block_string_size;
	int max_block_rows_pp, max_rows_pp, short_block_string_size, last_block_row_proc_id, last_block_row_in_current_pr;
	int small_block_row_width, small_block_size, current_pr_full_rows, last_block_row_width, matrix_size_current_pr;
	initParameters(matrix_side, block_side, total_pr, current_pr, 

	&total_block_rows, &total_full_block_rows, 
	&block_size, &block_string_size, 
	&max_block_rows_pp, &max_rows_pp, &short_block_string_size,
	&last_block_row_proc_id, &last_block_row_in_current_pr,
	&small_block_row_width, &small_block_size,
	&current_pr_full_rows, &last_block_row_width,
	&matrix_size_current_pr);

   	fp = fopen(name, "r");
   	if (!fp){
   		loc = 1;
   	}
   	MPI_Allreduce(&loc, &ret, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   	if (ret){
   		if(current_pr==0){
   			printf("Failed to open file\n");
   		}
       	return -1;
   	}
   	loc=0;
	for(i=0; i<matrix_side; i++){
		for(j=0;j<matrix_side;j++){
			block_number = j/block_side;
			block_number_proc_id = block_number % total_pr;
			block_number_per_process = block_number / total_pr;
			if(block_number_proc_id==current_pr){
				loc+=fscanf(fp, "%lf", a + i*max_rows_pp + block_number_per_process*block_side + j%block_side);
				count++;
			}
			else{
				loc+=fscanf(fp, "%lf", &trash);
			}
		}
		for(;count<(i+1)*max_rows_pp;count++){
			a[count]=0.;
		}
	}
	if(loc<matrix_side*matrix_side){
		loc=1;
	}
	else{
		loc=0;
	}
	MPI_Allreduce(&loc, &ret, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   	if (ret){
   		if(current_pr==0){
   			printf("Failed to read whole matrix\n");
   		}
       	return -1;
   	}
	fclose(fp);
   	return 0;
}

int initMatrixByRows(double *a, int n, int m, int p, int k){
	int i,j;
    	int pos_i=0;
    	int total_block_rows = (n+m-1)/m;
    	int max_block_rows = ((total_block_rows + p-1)/p);
    	int max_rows = max_block_rows*m;
    	for (i=0; i<max_rows; i++){
        	for (j=0; j<n; j++){
	           	pos_i = (i/m);
				pos_i = pos_i*p*m + i%m;
				pos_i += k*m;
	           	if (pos_i<n){
        	       	a[i*n+j]=f(pos_i,j);
            	}
	           	else a[i*n+j]=0.;
        	}
    	}
    	return 0;
}
int initMatrixByColumns(double *a, int n, int m, int p, int k){
	int i,j;
   	int pos_j=0;
   	int total_block_rows = (n+m-1)/m;
   	int max_block_rows = ((total_block_rows + p-1)/p);
   	int max_rows = max_block_rows*m;
   	for (i=0; i<n; i++){
       	for (j=0; j<max_rows; j++){
            pos_j = (j/m);
   			pos_j = pos_j*p*m + j%m;
			pos_j += k*m;
	        if (pos_j<n){
           		a[i*max_rows+j]=f(i,pos_j);
            }
	        else a[i*max_rows+j]=0.;
        }
    }
    return 0;
}

double MPI_getResidual(double *reinitialized, double *reversed, int matrix_side, int block_side, int total_pr, int current_pr, const int *blocks_order_reversed){
	double residual = -1.;
	double local_residual = -1.;

	double *temp_vect = 0;
	double *temp_block = 0;
	double *temp_block_for_multiplication = 0;
	double *buf = 0;
	int *blocks_order = 0;	

	int i, j, k, l, m, pos_j, pos_k;

	int total_block_rows, total_full_block_rows, block_size, block_string_size;
	int max_block_rows_pp, max_rows_pp, short_block_string_size, last_block_row_proc_id, last_block_row_in_current_pr;
	int small_block_row_width, small_block_size, current_pr_full_rows, last_block_row_width, matrix_size_current_pr;
	initParameters(matrix_side, block_side, total_pr, current_pr, 

	&total_block_rows, &total_full_block_rows, 
	&block_size, &block_string_size, 
	&max_block_rows_pp, &max_rows_pp, &short_block_string_size,
	&last_block_row_proc_id, &last_block_row_in_current_pr,
	&small_block_row_width, &small_block_size,
	&current_pr_full_rows, &last_block_row_width,
	&matrix_size_current_pr);
	
	int current_pr_rows = current_pr_full_rows * block_side;
	if ((small_block_row_width)&&(current_pr==last_block_row_proc_id)){
		current_pr_rows += small_block_row_width;
	}
	temp_vect = new double[max_rows_pp];
	temp_block = new double[block_size];
	temp_block_for_multiplication = new double[block_size];
	buf = new double[matrix_size_current_pr];
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	blocks_order = new int[total_block_rows];
	for(i=0; i<total_block_rows;i++){
		blocks_order[i]=blocks_order_reversed[blocks_order_reversed[i]];
	}

	for(i=0;i<max_rows_pp;i++){
		temp_vect[i]=0.;
	}
	
	for(i=0;i<total_pr;i++){
		if(current_pr==i){
			for(j=0;j<matrix_size_current_pr;j++){
				buf[j]=reinitialized[j];
			}
		}
		MPI_Bcast(buf, matrix_size_current_pr, MPI_DOUBLE, i, MPI_COMM_WORLD);

		for(j=0; j<max_block_rows_pp; j++){
			pos_j = blocks_order_reversed[j*total_pr + current_pr];
			for(k=0; k<max_block_rows_pp;k++){
				
				for(l=0;l<block_size;l++){
					temp_block[l]=0.;
				}
				
				for(l=0;l<total_full_block_rows;l++){
					simpleMatrixMultiply(reversed + j*block_string_size + l*block_size, 
					buf + l*short_block_string_size + k*block_size,
					temp_block_for_multiplication, block_side, block_side, block_side);
					addToMatrix(temp_block, temp_block_for_multiplication, block_size);
				}
				
				if(small_block_row_width){
					simpleMatrixMultiply(reversed + j*block_string_size + total_full_block_rows*block_size,
					buf + total_full_block_rows*short_block_string_size + k*small_block_size,
					temp_block_for_multiplication, block_side, small_block_row_width, block_side);
					addToMatrix(temp_block, temp_block_for_multiplication, block_size);
				}
				
				//try not to forget to subtract Id matrix!				
				pos_k = k*total_pr + i;
				if(pos_j==pos_k){
					if(pos_j<total_full_block_rows){
						for(l=0;l<block_side;l++){
							for(m=0;m<block_side;m++){
								if(l==m){
									temp_block[l*block_side+m]-=1;
								}
							}
						}
					}
					else{
						for(l=0;l<small_block_row_width;l++){
							for(m=0;m<small_block_row_width;m++){
								if(l==m){
									temp_block[l*block_side+m]-=1.;
								}
							}
						}	
					}
				}

				if((current_pr!=last_block_row_proc_id)||(pos_j<total_full_block_rows)){
					for(l=0; l<block_side; l++){
						for(m=0; m<block_side; m++){//5X combo!!!
							temp_vect[j*block_side+l]+=fabs(temp_block[l*block_side+m]);
						}
					}
				}
				else{
					for(l=0; l<small_block_row_width; l++){
						for(m=0; m<block_side; m++){
							temp_vect[j*block_side+l]+=fabs(temp_block[l*block_side+m]);
						}
					}
				}
			}
		}
	}
	
	local_residual = -1.;
	for(i=0; i<current_pr_rows;i++){
		if(local_residual<temp_vect[i]){
			local_residual = temp_vect[i];
		}
	}
#ifdef RESIDAL_PRINT	
	for(i=0;i<current_pr_rows;i++){
		pos_j = current_pr*block_side + (i/block_side)*total_pr*block_side + i%block_side;
		printf("%d row %.3le residual\n", pos_j, temp_vect[i]);
	}
#endif
    MPI_Allreduce(&local_residual, &residual, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	
	//delete[] temp_vect;
	//delete[] temp_block;
	//delete[] temp_block_for_multiplication;
	//delete[] buf;

	//delete[] blocks_order;

	return residual;
}

int initParameters(int matrix_side, int block_side, int total_pr, int current_pr, 
	int *total_block_rows, int *total_full_block_rows, 
	int *block_size, int *block_string_size, 
	int *max_block_rows_pp, int *max_rows_pp, int *short_block_string_size,
	int *last_block_row_proc_id, int *last_block_row_in_current_pr,
	int *small_block_row_width, int *small_block_size,
	int *current_pr_full_rows, int *last_block_row_width,
	int *matrix_size_current_pr){

	*total_block_rows = (matrix_side+block_side-1)/block_side;
   	*total_full_block_rows = matrix_side/block_side;

   	*block_size = block_side*block_side;
   	*block_string_size = matrix_side*block_side;

   	*max_block_rows_pp = (*total_block_rows + total_pr - 1)/total_pr;
   	*max_rows_pp = (*max_block_rows_pp) * block_side;
   	*short_block_string_size = (*max_block_rows_pp) * (*block_size);

   	*last_block_row_proc_id = (*total_block_rows-1)%total_pr;
  	*last_block_row_in_current_pr = (current_pr<=(*last_block_row_proc_id)) ? (*max_block_rows_pp) : (*max_block_rows_pp-1);
   	
   	*small_block_row_width = matrix_side - (*total_full_block_rows)*block_side;
   	*small_block_size = block_side*(*small_block_row_width);

   	*current_pr_full_rows = (((*last_block_row_proc_id)!=current_pr)||((*small_block_row_width)==0)) ? (*last_block_row_in_current_pr) : (*last_block_row_in_current_pr-1);//¿¿¿¿¿ ¿¿¿¿¿¿¿ ¿¿¿¿¿¿¿ ¿¿ ¿¿¿¿¿¿¿ ¿¿¿¿¿¿¿
   	*last_block_row_width = ((*current_pr_full_rows)!=(*last_block_row_in_current_pr)) ? (*small_block_row_width) : block_side;
   	
   	*matrix_size_current_pr = (*max_rows_pp)*matrix_side;	
   	return 0;
}


void MPI_printUpperLeftBlock(double *a, int matrix_side, int block_side, int total_pr, int current_pr, const int *blocks_order_reversed, 
	double *buf_string, double *buf_string_2){

	double *buf = 0;
	double *recvbuf = 0;
	double *sendbuf = 0;

	int corner_side = 7;
	MPI_Status status;
	int i, j, real_i;
	int pos_i, i_proc_id;
	int total_block_rows, total_full_block_rows, block_size, block_string_size;
	int max_block_rows_pp, max_rows_pp, short_block_string_size, last_block_row_proc_id, last_block_row_in_current_pr;
	int small_block_row_width, small_block_size, current_pr_full_rows, last_block_row_width, matrix_size_current_pr;
	initParameters(matrix_side, block_side, total_pr, current_pr, 

	&total_block_rows, &total_full_block_rows, 
	&block_size, &block_string_size, 
	&max_block_rows_pp, &max_rows_pp, &short_block_string_size,
	&last_block_row_proc_id, &last_block_row_in_current_pr,
	&small_block_row_width, &small_block_size,
	&current_pr_full_rows, &last_block_row_width,
	&matrix_size_current_pr);

	buf = new double[corner_side*corner_side]; 
	recvbuf = buf_string;
	sendbuf = buf_string_2;

	if (matrix_side<corner_side){
		corner_side = matrix_side;
	}

	int corner_total_block_rows = (corner_side + block_side - 1)/block_side;
	//int current_pr_block_rows = (corner_total_block_rows + total_pr - 1)/total_pr;
	//int sendbuf_size = corner_side*block_side;

	for(i = 0; i < corner_total_block_rows; i++){
		real_i = blocks_order_reversed[i];
		//pos_i = i/total_pr;
		//i_proc_id = i%total_pr;
		pos_i = real_i/total_pr;
		i_proc_id = real_i%total_pr;
		if(i_proc_id != 0){
			if(current_pr == i_proc_id){
				for(j=0;j<block_string_size;j++){
					sendbuf[j]=a[pos_i*block_string_size+j];
				}
				MPI_Send(sendbuf, block_string_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			}
			if(current_pr == 0){
				MPI_Recv(recvbuf, block_string_size, MPI_DOUBLE, i_proc_id, 0, MPI_COMM_WORLD, &status);
				for(j=0; j<block_string_size; j++){
					buf[i*block_string_size+j]=recvbuf[j];
				}
			}
		}
		else{
			for(j=0; j<block_string_size; j++){
					buf[i*block_string_size+j]=a[pos_i*block_string_size+j];
			}
		}
	}
	if(current_pr==0){
		printUpperLeftBlock(buf, matrix_side, block_side);
	}
	delete[] buf;
}