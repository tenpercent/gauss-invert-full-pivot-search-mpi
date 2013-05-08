#include <cstdlib>
#include <unistd.h>
#include "mpi.h"
#include "io.h"
#include "tools.h"
#include "gauss-invert.h"
using namespace std;
int main(int argc, char *argv[]){
        double *a=0;
        double *b=0;
        int *blocks_order_reversed = 0;
        int *blocks_order = 0;
        double *buf_1 = 0;
        double *buf_2 = 0;
        double *buf_string = 0;
        double *buf_string_2 = 0;

        double time_w1=0.;
        double time_w2=0.;
        int current_pr;
        int total_pr;
        int init_method = 0;
        double residual = -1.;
        int i, j, matrix_side, block_side;
        char *name = 0;
        int total_block_rows, total_full_block_rows, block_size, block_string_size;
        int max_block_rows_pp, max_rows_pp, short_block_string_size, last_block_row_proc_id, last_block_row_in_current_pr;
        int small_block_row_width, small_block_size, current_pr_full_rows, last_block_row_width, matrix_size_current_pr;
        int pos_i=0;

    	MPI_Init(&argc, &argv);

    	MPI_Comm_rank(MPI_COMM_WORLD, &current_pr);//номер узла
    	MPI_Comm_size(MPI_COMM_WORLD, &total_pr);//число узлов

    	/*
    	Usage:  ./a.out matrix_size block_size
            ./a.out filename matrix_size block_size
    	*/
    	if (argc==3){
        	matrix_side = atoi(argv[1]);
        	block_side = atoi(argv[2]);
            init_method = 1;
    	}
    	else if (argc==4){
        	name = argv[1];
        	matrix_side = atoi(argv[2]);
        	block_side = atoi(argv[3]);
            init_method = 2;
    	}
    	else{
        	printf("Usage:\n %s <matrix_size> <block_size> <total_processes>\n %s file_name <matrix_size> <block_size>\n",
                argv[0],
                argv[0]);
    	return -1;
    	}
        
        initParameters(matrix_side, block_side, total_pr, current_pr, 
        &total_block_rows, &total_full_block_rows, 
        &block_size, &block_string_size, 
        &max_block_rows_pp, &max_rows_pp, &short_block_string_size,
        &last_block_row_proc_id, &last_block_row_in_current_pr,
        &small_block_row_width, &small_block_size,
        &current_pr_full_rows, &last_block_row_width,
        &matrix_size_current_pr);

    	a = new double[matrix_size_current_pr];
    	b = new double[matrix_size_current_pr];

        blocks_order_reversed = new int[total_block_rows];
        blocks_order = new int[total_block_rows];

        buf_1 = new double[block_size];
        buf_2 = new double[block_size];

        buf_string = new double[2*block_string_size];
        buf_string_2 = new double[2*block_string_size];

        for (i = 0; i < total_block_rows; i++)
        {
            blocks_order_reversed[i]=i;
            blocks_order[i]=i;
        }
    	
        for (i=0; i<matrix_size_current_pr; i++){
           	    a[i]=0.;
            	b[i]=0.;
    	}
    	if (init_method==1){
        	initMatrixByRows(a, matrix_side, block_side, total_pr, current_pr);
    	}
    	else if (init_method==2){
            readMatrixByRows(name, a, matrix_side, block_side, total_pr, current_pr);
    	}
    	for (i=0; i<max_rows_pp; i++){
        	for (j=0; j<matrix_side; j++){
            		pos_i = (i/block_side);
                    pos_i = pos_i*total_pr*block_side + i%block_side;
                    pos_i += current_pr*block_side;
	            	if (pos_i<matrix_side){
        	        	b[i*matrix_side+j]=identity(pos_i,j);
            		}
            		else b[i*matrix_side+j]=0.;
        	}
    	}
    	makeBlockMatrix_Rows(a, matrix_side, block_side, total_pr, current_pr);
    	makeBlockMatrix_Rows(b, matrix_side, block_side, total_pr, current_pr);
    	
        MPI_printUpperLeftBlock(a, matrix_side, block_side, total_pr, current_pr, blocks_order_reversed, buf_string, buf_string_2);
        MPI_printUpperLeftBlock(b, matrix_side, block_side, total_pr, current_pr, blocks_order_reversed, buf_string, buf_string_2);

        MPI_Barrier(MPI_COMM_WORLD);
        time_w1=MPI_Wtime();

        gaussInvert(a, b, 
            matrix_side, block_side, total_pr, current_pr, 
            blocks_order_reversed, blocks_order, 
            buf_1, buf_2, buf_string, buf_string_2);
        
        MPI_Barrier(MPI_COMM_WORLD);
        time_w2=MPI_Wtime();
        
        MPI_printUpperLeftBlock(b, matrix_side, block_side, total_pr, current_pr, blocks_order_reversed, buf_string, buf_string_2);
        
        if(init_method==1){
            initMatrixByColumns(a, matrix_side, block_side, total_pr, current_pr);
        }
        if(init_method==2){
            readMatrixByColumns(name, a, matrix_side, block_side, total_pr, current_pr);
        }
        
        makeBlockMatrix_Columns(a, matrix_side, block_side, total_pr, current_pr);

        residual = MPI_getResidual(a, b, matrix_side, block_side, total_pr, current_pr, blocks_order);
        if(current_pr==0){
            
            printf("Residual: %.3le\n", residual);
            printf("Time elapsed: %.3le\n", time_w2-time_w1);
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);

      	delete[] a;
      	delete[] b;

        delete[] blocks_order_reversed;
        delete[] blocks_order;

        delete[] buf_1;
        delete[] buf_2;
        delete[] buf_string;
        delete[] buf_string_2;

    	MPI_Finalize();
    	return 0;
}
