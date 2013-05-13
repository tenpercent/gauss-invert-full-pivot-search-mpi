#include "io.h"
#include "tools.h"
#include <math.h>

double MPI_getResidual_rewrite(double *reinitialized, double *reversed, int matrix_side, int block_side, int total_pr, int current_pr, const int *blocks_order_reversed,
	double *buf_string, double *buf_1, double *buf_2){
	double residual = -1.;
	double local_residual = -1.;

	double *temp_vect = 0;
	double *temp_block = 0;
	double *temp_block_for_multiplication = 0;
	double *buf = 0;
	int *blocks_order = 0;	

	int i, j, k, l, m, pos_j, pos_k;

	int current_row_proc_id, local_i;

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

	temp_vect = buf_string;
	//careful here
	temp_block = buf_1;
	temp_block_for_multiplication = buf_2;

	buf = new double[block_string_size]; 

	blocks_order = new int[total_block_rows];
	for(i=0; i<total_block_rows;i++){
		blocks_order[i]=blocks_order_reversed[blocks_order_reversed[i]];
	}

	for(i=0;i<max_rows_pp;i++){
		temp_vect[i]=0.;
	}
	
	for(i=0;i<total_block_rows;i++){
		current_row_proc_id = i % total_pr;
		local_i = (i+total_pr-1-current_pr)/total_pr;
		
		if(current_pr==current_row_proc_id){
			for(j=0;j<total_full_block_rows;j++){
				for(k=0; k<block_size; k++){
					buf[j*block_size+k]=reinitialized[j*short_block_string_size+local_i*block_size+k];
				}
			}
			if(small_block_row_width){
				for(k=0; k<small_block_size; k++){
					buf[total_full_block_rows*block_size+k]=reinitialized[total_full_block_rows*short_block_string_size+local_i*small_block_size+k];
				}
			}
		}
		
		MPI_Bcast(buf, block_string_size, MPI_DOUBLE, current_row_proc_id, MPI_COMM_WORLD);
		
		for(j=0; j<max_block_rows_pp; j++){
			pos_j = blocks_order_reversed[j*total_pr + current_pr];
			
			for(k=0;k<block_size;k++){
				temp_block[k]=0.;
			}
				
			for(k=0;k<total_full_block_rows;k++){
				simpleMatrixMultiply(reversed + j*block_string_size + k*block_size, 
				buf + k*block_size,
				temp_block_for_multiplication, block_side, block_side, block_side);
				addToMatrix(temp_block, temp_block_for_multiplication, block_size);
			}
				
			if(small_block_row_width){
				simpleMatrixMultiply(reversed + j*block_string_size + total_full_block_rows*block_size,
				buf + total_full_block_rows*block_size,
				temp_block_for_multiplication, block_side, small_block_row_width, block_side);
				addToMatrix(temp_block, temp_block_for_multiplication, block_size);
			}
				
			if(pos_j==i){
				if(pos_j<total_full_block_rows){
					for(l=0;l<block_side;l++){
						for(m=0;m<block_side;m++){
							if(l==m){
								temp_block[l*block_side+m]-=1.;	
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
			
			//if((current_pr!=last_block_row_proc_id)||(pos_j<total_full_block_rows)){
			if((pos_j<total_full_block_rows)){
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

	local_residual = -1.;
	for(i=0; i<current_pr_rows;i++){
		if(local_residual<temp_vect[i]){
			local_residual = temp_vect[i];
		}
	}
#ifdef RESIDAL_PRINT	
	for(i=0;i<current_pr_rows;i++){
		pos_j = current_pr*block_side + (i/block_side)*total_pr*block_side + i%block_side;
		printf("%d row %.3le residual %d process\n", pos_j, temp_vect[i], current_pr);
	}
#endif
    MPI_Allreduce(&local_residual, &residual, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	
	delete[] buf;
	delete[] blocks_order;
	return residual;
}