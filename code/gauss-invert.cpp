#include "tools.h"
#include "io.h"
#include <mpi.h>
#include <stddef.h>

typedef struct {
    int label;
    int min_k;
    int rank;
    double minnorm;
} mainBlockInfo;

void searchMainBlock(void *inv, void *inoutv, int *len, MPI_Datatype *MPI_mainBlockInfo){
  mainBlockInfo *in = (mainBlockInfo *) inv;
  mainBlockInfo *inout = (mainBlockInfo *) inoutv;
  int i;
  for (i=0; i<*len; ++i){
    if(in->label){
      if(inout->label){
        if(in->minnorm<inout->minnorm){
          inout->minnorm = in->minnorm;
          inout->rank = in->rank;
          inout->min_k = in->min_k;
        }
      }
      else{
        inout->minnorm = in->minnorm;
        inout->rank = in->rank;
        inout->min_k = in->min_k;
      }
    }
      inout->label += in->label;
     // in++; 
     // inout++;
  }
}

int gaussInvert(double *a, double *b, int matrix_side, int block_side, 
  int total_pr, int current_pr, 
  int* blocks_order_reversed, int* blocks_order,
  double* buf_1, double* buf_2, double* buf_string, double* buf_string_2){
  MPI_Status status;

  int first_row, first_row_proc_id, last_row_c;
  int current_row, current_row_proc_id;

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

  int i, j, k, min_j, min_k_global;
  int res;

  int buf_size = 2 * block_string_size;

	double temp=-1.;

  struct MPI_Double_Int{
   	double minnorm;
   	int rank;
 	};

  mainBlockInfo in, out;

  MPI_Op MPI_searchMainBlock;
  MPI_Datatype MPI_mainBlockInfo;

  MPI_Datatype type[4] = { MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE };

  int blocklen[4] = { 1, 1, 1 };

  MPI_Aint displ[4] = { 0, 0, 0, 0 };

  displ[0] = offsetof(mainBlockInfo, label);
  displ[1] = offsetof(mainBlockInfo, min_k);
  displ[2] = offsetof(mainBlockInfo, rank);
  displ[3] = offsetof(mainBlockInfo, minnorm);

  MPI_Type_create_struct(4, blocklen, displ, type, &MPI_mainBlockInfo);

  MPI_Type_commit(&MPI_mainBlockInfo);

  MPI_Op_create(searchMainBlock, 1, &MPI_searchMainBlock);

 	in.rank = current_pr;
	in.minnorm = 1000000000.;
  in.label = 0;
  in.min_k = 0;

 	for (i=0; i<buf_size; i++){
  	buf_string[i] = 0.;
		buf_string_2[i] = 0.;
 	}
 	for (i=0; i<total_full_block_rows; i++){
		first_row = (i+total_pr-1-current_pr)/total_pr;
		first_row_proc_id = i%total_pr;
    min_j = 0;
    min_k_global = 0;

   	in.minnorm = 0.;
   	in.rank = current_pr;
    in.label = 0;
    in.min_k = i;
		temp = 0.;
	
	  for (j=first_row; j<current_pr_full_rows; j++){
   		for (k=i; k<total_full_block_rows; k++){
   	    res = simpleInvert(a + j*block_string_size + k*block_size, buf_1, buf_2, block_side);
        if (!res) {
   			  temp = matrixNorm(buf_1, block_side);
          if (in.label){
            if (temp<in.minnorm){
     		      in.minnorm = temp;
              min_j=j;
              in.min_k = k;
            }
          }
          else{
            in.label = 1;
            in.minnorm = temp;
            min_j=j;
            in.min_k = k;
          }
        }
      }
    }

#ifdef WO_PIVOT_SEARCH_ATALL
		if (current_pr==first_row_proc_id){
			min_j=first_row;
			in.minnorm = 1.;
			in.min_k=i;
			in.label=1;
		}
#endif

    MPI_Allreduce(&in, &out, 1, MPI_mainBlockInfo, MPI_searchMainBlock, MPI_COMM_WORLD);

		if (out.label==0){
			if (current_pr==0){
			  printf("Main block not found!\n\t -- Step %d\n", i);
			}
			fflush(stdout);
			MPI_Barrier(MPI_COMM_WORLD);
			return -1;
		}

#ifdef DEBUG_MODE
		if (current_pr == first_row_proc_id){
      printf("**\nOUT RANK %d\n", out.rank);
      printf("IN RANK %d\n**\n", in.rank);
     	fflush(stdout);
    }
#endif
   	min_k_global = out.min_k;
#ifdef W_FULL_PIVOT_SEARCH
  	for (j=0; j<max_block_rows_pp; j++){
			swapMatrix(a + j*block_string_size + i*block_size, a + j*block_string_size + min_k_global*block_size, block_size);
 		}
#endif
    temp = blocks_order[i];
    blocks_order[i]=blocks_order[min_k_global];
    blocks_order[min_k_global]=temp;

		//for debug purposes:
#ifdef WO_PIVOT_SEARCH_ATALL
		out.rank = first_row_proc_id;
		//don't forget about it
#endif
      /***************Multiply string by the main block*****************/
    if (current_pr==out.rank){
      simpleInvert(a + min_j*block_string_size + i*block_size, buf_1, buf_2, block_side);
      for (j=i+1; j<total_full_block_rows; j++){
        simpleMatrixMultiply(buf_1, a + min_j*block_string_size + j*block_size, buf_2, block_side, block_side, block_side);
        copyMatrix(buf_2, a + min_j*block_string_size + j*block_size, block_size);
      }
      for (j=0; j<total_full_block_rows; j++){
        simpleMatrixMultiply(buf_1, b + min_j*block_string_size + j*block_size, buf_2, block_side, block_side, block_side);
        copyMatrix(buf_2, b + min_j*block_string_size + j*block_size, block_size);
      }
      if(small_block_row_width){
        simpleMatrixMultiply(buf_1, b + min_j*block_string_size + total_full_block_rows*block_size, buf_2, block_side, block_side, small_block_row_width);
        copyMatrix(buf_2, b + min_j*block_string_size + total_full_block_rows*block_size, small_block_size);

        simpleMatrixMultiply(buf_1, a + min_j*block_string_size + total_full_block_rows*block_size, buf_2, block_side, block_side, small_block_row_width);
        copyMatrix(buf_2, a + min_j*block_string_size + total_full_block_rows*block_size, small_block_size);
      }
			for (j=0; j<block_string_size; j++){
		    buf_string[j]=a[min_j*block_string_size+j];
			}
			for (j=0; j<block_string_size; j++){
				buf_string[j+block_string_size]=b[min_j*block_string_size+j];
			}
    }

		MPI_Bcast(buf_string, buf_size, MPI_DOUBLE, out.rank, MPI_COMM_WORLD);

		if (out.rank!=first_row_proc_id){
			if (current_pr==first_row_proc_id){
				for (j=0; j<block_string_size; j++){
					buf_string_2[j]=a[first_row*block_string_size+j];
          a[first_row*block_string_size+j]=buf_string[j];
				}
				for (j=0; j<block_string_size; j++){
					buf_string_2[j+block_string_size]=b[first_row*block_string_size+j];
          b[first_row*block_string_size+j]=buf_string[j+block_string_size];
				}
				MPI_Send(buf_string_2, buf_size, MPI_DOUBLE, out.rank, 42, MPI_COMM_WORLD);
			}
			if(current_pr==out.rank){
				MPI_Recv(buf_string_2, buf_size, MPI_DOUBLE, first_row_proc_id, 42, MPI_COMM_WORLD, &status);
				for (j=0; j<block_string_size; j++){
					a[min_j*block_string_size+j]=buf_string_2[j];
				}
				for (j=0; j<block_string_size; j++){
					b[min_j*block_string_size+j]=buf_string_2[j+block_string_size];
				}
			}
		}
		else{
			if(current_pr==out.rank){
				if (min_j!=first_row){
					for(j=0; j<block_string_size; j++){
						temp=a[min_j*block_string_size+j];
						a[min_j*block_string_size+j]=a[first_row*block_string_size+j];
						a[first_row*block_string_size+j]=temp;
					}
					for(j=0; j<block_string_size; j++){
						temp=b[min_j*block_string_size+j];
						b[min_j*block_string_size+j]=b[first_row*block_string_size+j];
						b[first_row*block_string_size+j]=temp;
					}
    		}
			}
		}

#ifdef DEBUG_MODE
		if (current_pr == out.rank){
			printf("SUCCESS SEARCHING MAIN BLOCK STEP %d\n", i);
			fflush(stdout);
		}
#endif

    if (current_pr == first_row_proc_id){
    	first_row++;
  	}

  	for (j=first_row; j<max_block_rows_pp; j++){
      for (k = i+1; k<total_full_block_rows; k++){
		    simpleMatrixMultiply(a + j*block_string_size + i*block_size, buf_string + k*block_size, buf_2, block_side, block_side, block_side);
			  subtractFromMatrix(a + j*block_string_size + k*block_size, buf_2, block_size);
    	}
      for (k = 0; k<total_full_block_rows; k++){
				simpleMatrixMultiply(a + j*block_string_size + i*block_size, buf_string + block_string_size + k*block_size, buf_2, block_side, block_side, block_side);
				subtractFromMatrix(b + j*block_string_size + k*block_size, buf_2, block_size);
      }
      if(small_block_row_width){
        simpleMatrixMultiply(a + j*block_string_size + i*block_size, buf_string + total_full_block_rows*block_size, buf_2, block_side, block_side, small_block_row_width);
        subtractFromMatrix(a + j*block_string_size + total_full_block_rows*block_size, buf_2, small_block_size);

        simpleMatrixMultiply(a + j*block_string_size + i*block_size, buf_string + total_full_block_rows*block_size + block_string_size, buf_2, block_side, block_side, small_block_row_width);
        subtractFromMatrix(b + j*block_string_size + total_full_block_rows*block_size, buf_2, small_block_size);
    	}
   	}
	}

  if(small_block_row_width){
    if (current_pr==last_block_row_proc_id){
      simpleInvert(a + current_pr_full_rows*block_string_size + total_full_block_rows*block_size, buf_1, buf_2, small_block_row_width);
      for (k=0; k<total_full_block_rows; k++){
        simpleMatrixMultiply(buf_1, b + current_pr_full_rows*block_string_size + k*block_size, buf_2, 
          small_block_row_width, small_block_row_width, block_side);
        copyMatrix(buf_2, b + current_pr_full_rows*block_string_size + k*block_size, small_block_size);
    	}
      simpleMatrixMultiply(buf_1, b + current_pr_full_rows*block_string_size + total_full_block_rows*block_size, buf_2, 
      	small_block_row_width, small_block_row_width, small_block_row_width);
      copyMatrix(buf_2, b + current_pr_full_rows*block_string_size + total_full_block_rows*block_size, small_block_row_width*small_block_row_width);
    }
  }
#ifdef DEBUG_MODE
  if (current_pr==first_row_proc_id){
  	printf("SUCCESS IN DIRECT FLOW!!!\n");
  	fflush(stdout);
	}
#endif

  for(j=0; j<total_full_block_rows; j++){
    blocks_order_reversed[blocks_order[j]]=j;
  }

#ifdef W_REVERSE_FLOW
	if(small_block_row_width){
    if (current_pr==last_block_row_proc_id){
      for (j=0; j<block_string_size; j++){
        buf_string[j]=b[(last_block_row_in_current_pr-1)*block_string_size + j];
      }
    }
    MPI_Bcast(buf_string, block_string_size, MPI_DOUBLE, last_block_row_proc_id, MPI_COMM_WORLD);
    last_row_c = (current_pr==last_block_row_proc_id) ? (last_block_row_in_current_pr-1) : (last_block_row_in_current_pr);
    for (j=last_row_c-1;j>=0;j--){
      for (k=0; k<total_full_block_rows;k++){
        simpleMatrixMultiply(a+j*block_string_size+total_full_block_rows*block_size, buf_string + k*block_size, buf_2,block_side,small_block_row_width,block_side);
        subtractFromMatrix(b+j*block_string_size+k*block_size, buf_2, block_size);
      } 
      simpleMatrixMultiply(a+j*block_string_size+total_full_block_rows*block_size, buf_string + total_full_block_rows*block_size, buf_2,block_side,small_block_row_width,small_block_row_width);
      subtractFromMatrix(b+j*block_string_size+total_full_block_rows*block_size, buf_2, small_block_size);
    }
  }	
  for (i=total_full_block_rows-1; i>0; i--){//i-ю строчку вычитаем из всех
    //current_row = (i+total_pr-1-current_pr)/total_pr;//first row in cur_pr not upper than the subtracted string
    current_row_proc_id = i%total_pr;    
    current_row = i/total_pr;
    if(current_pr==current_row_proc_id){
    	for (j=0; j<block_string_size; j++){
       	buf_string[j]=b[current_row*block_string_size + j];
  		}
  	}
        
    if(current_pr<current_row_proc_id){
      current_row++;
    }
        
    MPI_Bcast(buf_string, block_string_size, MPI_DOUBLE, current_row_proc_id, MPI_COMM_WORLD);

    for(j=current_row-1; j>=0; j--){//из j-й строчки правой матрицы вычитается
      for (k=0; k<total_full_block_rows;k++){
        simpleMatrixMultiply(a+j*block_string_size+i*block_size, buf_string + k*block_size, buf_2, block_side, block_side, block_side);
      	subtractFromMatrix(b+j*block_string_size+k*block_size, buf_2, block_size);
      }
      if (small_block_row_width){
        simpleMatrixMultiply(a+j*block_string_size+i*block_size, buf_string + total_full_block_rows*block_size, buf_2,block_side,block_side,small_block_row_width);
        subtractFromMatrix(b+j*block_string_size+total_full_block_rows*block_size, buf_2, small_block_size);
    	}
    }
  }
#endif
#ifdef DEBUG_MODE
  if (current_pr==0){
    printf("SUCCESS IN REVERSE FLOW!!!\n");
    fflush(stdout);
    printf("Exit from gauss\n");
    fflush(stdout);
  }
#endif
  
  MPI_Type_free(&MPI_mainBlockInfo);
  MPI_Op_free(&MPI_searchMainBlock);
	return 0;
}

