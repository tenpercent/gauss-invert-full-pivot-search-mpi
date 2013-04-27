#include "tools.h"
#include "io.h"
#include <cmath>

void blocksMultiply(double *a, int p, int q, double *b, int i, int j, double *out, int k, int m, int n, int r){
  //out=a_{pq}*b_{ij}
  //k = n/m
  //r = m%n
  //m - размер блока
  int x, y, z, t;
  x = (p<k) ? m:r;
  y = (q<k) ? m:r;
  z = (j<k) ? m:r;
  t = (i<k) ? m:r;
  //int block_size_first = m*x;
  //int block_size_second = m*y;
  simpleMatrixMultiply(a + p*m*n + q*m*x, b + i*m*n + j*m*t, out, x, y, z);
}

void blocksMultiplyLast(double *a, double *b, int i, int j, double *out, int k, int m, int n, int r){
  //out=a_{pq}*b_{ij}
  //k = n/m
  //r = m%n
  //m - размер блока
  int z;
  z = (j<k) ? m:r;
  //int block_size_first = m*x;
  //int block_size_second = m*y;
  simpleMatrixMultiply(a, b + i*m*n + j*m*r, out, r, r, z);
}

void blocksSubtract(double *a, int i, int j, double *b, int k, int m, int n, int r){
  //a_{ij}-=b;
  int x, y, size;
  x = (i<k)? m:r;
  y = (j<k)? m:r;
  size = x*y;
  subtractFromMatrix(a + i*m*n + j*m*x, b, size);
}

void blocksCopy(double *a, int i, int j, double *b, int k, int m, int n, int r){
  //a_{ij}=b;
  int x, y, size;
  x = (i<k)? m:r;
  y = (j<k)? m:r;
  size = x*y;
  copyMatrix(b, a + i*m*n + j*m*x, size);
}
void blocksSwap(double *a, int i, int j, double *b, int p, int q, int k, int m, int n, int r){
  //a_{ij}=b;
  int x, y, z, size;
  x = (i<k)? m:r;
  y = (j<k)? m:r;
  z = (p<k)? m:r;
  size = x*y;
  swapMatrix(b + p*m*n + q*m*z, a + i*m*n + j*m*x, size);
}

void printMatrix(double* const a, int h, int l){
  int point=0;
  for (int j=0; j<2*l; j++){
    printf(" _");
  }
  printf("\n ");
  for (int i=0; i<h; i++){
    printf("|");
    for(int j=0; j<l; j++){
      printf("%8.4lf", a[point]);
      point++;
    }
    printf("|\n ");
  }
  for (int j=0; j<2*l; j++){
    printf("_ ");
  }
  printf("\n");
  //return 0;
}
int printBlockMatrix(double* const a, int n, int m){
  const int s = n/m;
  const int r = (n - s*m);
  const int string_size = m*n;
  const int block_size = m*m;
  const int small_block_size = m*r;
  if(r!=0){
    for (int i=0; i<s; i++){
      for (int j=0; j<s; j++){
	printMatrix(a+i*string_size+j*block_size, m, m);
      }
      printMatrix(a+i*string_size+s*block_size, m, r);
    }
    for (int j=0; j<s; j++){
      printMatrix(a+s*string_size+j*small_block_size, r, m);
    }
    printMatrix(a+s*string_size+s*small_block_size, r, r);
  }
  else{
    for (int i=0; i<s; i++){
      for (int j=0; j<s; j++){
	printMatrix(a+i*string_size+j*block_size, m, m);
      }
    }
  }
  return 0;
}
int makeBlockMatrix_Rows(double* const a, int matrix_side, int block_side, int total_pr, int current_pr){
  if ((matrix_side<1)||(block_side<1)){
    printf("incorrect size of a matrix or a block!\n\t --makeBlockMatrix");
    return -1;
  }
  int i, j, k, l;
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

  int count=0;
  int res = 0;
  double* const b = new double[matrix_size_current_pr];

  if (small_block_row_width!=0){
    for (i=0; i<max_block_rows_pp; i++){//цикл по блочным строчкам
      for (j=0; j<total_full_block_rows; j++){//цикл по блокам внутри строки
        for (k=0; k<block_side; k++){
            for (l=0; l<block_side; l++){
                b[count]=a[i*block_string_size+j*block_side+k*matrix_side+l];
                count++;
            }
        }
      }
      for (k=0; k<block_side; k++){
        for (l=0; l<small_block_row_width; l++){
            b[count]=a[i*block_string_size+total_full_block_rows*block_side+k*matrix_side+l];
            count++;
        }
      }
    }
    res = copyMatrix(b, a, matrix_size_current_pr);
    if (res!=0){
      printf("copyMatrix failed!\nLet me sing a song about it\n\t -- makeBlockMatrix\n");
    }
  }
  else{
    for (i=0; i<max_block_rows_pp; i++){//цикл по блочным строчкам
      for (j=0; j<total_full_block_rows; j++){//цикл по блокам внутри строки
        for (k=0; k<block_side; k++){
            for (l=0; l<block_side; l++){
                b[count]=a[i*block_string_size+j*block_side+k*matrix_side+l];
                count++;
            }
        }
      }
    }
    res = copyMatrix(b, a, matrix_size_current_pr);
    if (res!=0){
      printf("copyMatrix failed!\nLet me sing a song about it\n\t -- makeBlockMatrix\n");
    }
  }

  delete[] b;
  return 0;
}
int makeBlockMatrix_Columns(double* const a, int matrix_side, int block_side, int total_pr, int current_pr){
 if ((matrix_side<1)||(block_side<1)){
    printf("incorrect size of a matrix or a block!\n\t --makeBlockMatrix");
    return -1;
  }
  int i, j, k, l;
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

  int count=0;
  int res = 0;
  double* const b = new double[matrix_size_current_pr];

  if (small_block_row_width!=0){
    for (i=0; i<total_full_block_rows; i++){//цикл по блочным строчкам
      for (j=0; j<max_block_rows_pp; j++){//цикл по блокам внутри строки
        for (k=0; k<block_side; k++){
            for (l=0; l<block_side; l++){
                b[count]=a[i*short_block_string_size+j*block_side+k*max_rows_pp+l];
                count++;
            }
        }
      }
    }
    for(j=0; j<max_block_rows_pp;j++){
      for (k=0; k<small_block_row_width; k++){
        for (l=0; l<block_side; l++){
          b[count]=a[total_full_block_rows*short_block_string_size+j*block_side+k*max_rows_pp+l];
          count++;
        }
      }
    }
    res = copyMatrix(b, a, matrix_size_current_pr);
    if (res!=0){
      printf("copyMatrix failed!\nLet me sing a song about it\n\t -- makeBlockMatrix\n");
    }
  }
  else{
    for (i=0; i<total_full_block_rows; i++){//цикл по блочным строчкам
      for (j=0; j<max_block_rows_pp; j++){//цикл по блокам внутри строки
        for (k=0; k<block_side; k++){
            for (l=0; l<block_side; l++){
                b[count]=a[i*short_block_string_size+j*block_side+k*max_rows_pp+l];
                count++;
            }
        }
      }
    }
    res = copyMatrix(b, a, matrix_size_current_pr);
    if (res!=0){
      printf("copyMatrix failed!\nLet me sing a song about it\n\t -- makeBlockMatrix\n");
    }
  }

  delete[] b;
  return 0;
}
int makeOrdinaryMatrix(double* const a, int n, int m){
  if ((m<1)||(n<1)){
    printf("bad sizes \n\t --makeOrdinaryMatrix\n");
    return -1;
  }
  if ((n<1)||(m<1)){
    printf("incorrect size of a matrix or a block!\n\t --makeBlockMatrix");
    return -1;
  }
  const int s = n/m;
  const int r = (n - m*s);
  const int matrix_size = n*n;
  //int block_size = m*m;
  const int string_size = n*m;
  int count=0;
  int res = 0;
  double* const b = new double[matrix_size];
  if (r!=0){
    for (int i=0; i<s; i++){//цикл по блочным строчкам
      for (int j=0; j<s; j++){//цикл по блокам внутри строки
	for (int k=0; k<m; k++){
	  for (int l=0; l<m; l++){//for for for for
	    b[i*string_size+j*m+k*n+l]=a[count];
	    count++;
	  }
	}
      }
      for (int k=0; k<m; k++){
	for (int l=0; l<r; l++){
	  b[i*string_size+s*m+k*n+l]=a[count];
	  count++;
	}
      }
    }
    for (int j=0; j<s; j++){
      for (int k=0; k<r; k++){
	for (int l=0; l<m; l++){
	  b[s*string_size+j*m+k*n+l]=a[count];
	  count++;
	}
      }
    }
    for (int k=0; k<r; k++){
      for (int l=0; l<r; l++){
	b[s*string_size+s*m+k*n+l]=a[count];
	count++;
      }
    }
    res = copyMatrix(b, a, matrix_size);
    if (res!=0){
      printf("copyMatrix failed!\nLet me sing a song about it\n\t -- makeBlockMatrix\n");
    }
  }
  else{
    for (int i=0; i<s; i++){//цикл по блочным строчкам
      for (int j=0; j<s; j++){//цикл по блокам внутри строки
	for (int k=0; k<m; k++){
	  for (int l=0; l<m; l++){//for for for for
	    b[i*string_size+j*m+k*n+l]=a[count];
	    count++;
	  }
	}
      }
    }
    res = copyMatrix(b, a, matrix_size);
    if (res!=0){
      printf("copyMatrix failed!\nLet me sing a song about it\n\t -- makeBlockMatrix\n");
    }
  }
  delete[] b;
  return 0;
}
int copyMatrix(double* const in, double* const out, int size){
  if (!size){
    printf("tried to copy nothing!\n\t -- copyMatrix\n");
    return -1;
  }
  int i;
  //double tempi;
  for (i=0; i<size; i++){
    //tempi=in[i];
    //out[i]=tempi;// мало ли...
    out[i] = in[i];
  }
  return 0;
}

int zeroMatrix(double* const a, int size){
  if (!size){
    printf("Bad matrix size!\n\t--zeroMatrix\n");
    return -1;
  }
  int i;
  for (i=0; i<size; i++){
    a[i]=0.;
  }
  return 0;
}
int idMatrix(double* const a, int m){
  if (!m){
    printf("Bad matrix size!\n\t--idMatrix\n");
    return -1;
  }
  int i;
  //int res;
  int size = m*m;
  //res = zeroMatrix(a, size);
  //if (res!=0){
  //  printf("zeroMatrix failed!\n\t--idMatrix\n");
  //}
  for (i=0; i<size; i++){
    *(a+i)=0.;
  }
  for (i=0; i<size; i+=(m+1)){
    a[i]=1.;
  }
  // printf("idMatrix says:\n");
  // printMatrix(a, m, m);
  return 0;
}

int swapMatrix(double* const a, double* const b, int size){
  double tempo;
  int i;
  if (!size){
    printf("Trying to swap empty matrices, bad boy\n");
    return -1;
  }
  for (i=0; i<size; i++){
    tempo = a[i];
    a[i] = b[i];
    b[i] = tempo;
  }
  return 0;
}
int addToMatrix(double* const a, double* const b, int size){
  if (size<1){
    printf("Trying to add empty matrices, bad boy\n");//случай, когда матрицы разного размера, принципиально не отслеживается
    return -1;
  }

  for(int i=0; i<size; i++){
    a[i]+=b[i];
  }
  return 0;
}
int subtractFromMatrix(double* const a, double* const b, int size){
  //if (size<1){
  int i;
  if (!size){
    printf("Trying to subtract empty matrices, bad boy\n");
    return -1;
  }
  for (i=0; i<size; i++){
    a[i]-=b[i];
  }
  return 0;
}

double matrixNorm(double* const a, int m){
  //здесь мы будем как-нибудь считать норму m*m матрицы
  //if (m<1){
  if (!m){
    printf("Can't find a norm of empty matrix\n");
    return -1;
  }
  double temp = 0.;
  double ans = 0.;
  int i, j;
  for (i=0; i<m; i++){
    temp = 0.;
    for (j=0; j<m; j++){
      temp+=fabs(a[j*m+i]);
    }
    if (temp>ans){
      ans=temp;
    }
  }
  return ans;
}
/*
int simpleMatrixMultiply(double* const a, double* const b, double* const out, int p, int q, int r){
  //double* c = new double[p*r];
  //if ((p<1)||(q<1)||(r<1)){
  if (!(p&&q&&r)){
    printf("Bad matrix; can't multiply\n\t--simpleMAtrixMultiply\n");
    return -1;
  }
  for (int i=0; i<p; i++){
    for (int j=0; j<r; j++){
      out[i*r+j]=0.;
      for (int k=0; k<q; k++){
	out[i*r+j]+=(a[i*q+k]*b[k*r+j]);
      }
    }
  }
  return 0;
}*/

int blockMatrixMultiply(double* a, double* b, double* out, int n, int m){
  const int s = n/m;
  const int r = n - s*m;
  const int matrix_size = n*n;
  const int block_size = m*m;
  const int string_size = m*n;
  const int small_block_size = r*m;
  const int smallest_block_size = r*r;

  int i, j, k;
  int res = zeroMatrix(out, matrix_size);
  if (res!=0){
    printf("zeroMatrix failed!\n\t--blockmatrixmultiply\n");
  }
  double* const temp_block = new double[block_size];
  if (r>0){
    for (i=0; i<s; i++){
      for (j=0; j<s; j++){
	for (k=0; k<s; k++){
	  res = simpleMatrixMultiply(a+i*string_size+k*block_size, b+k*string_size+j*block_size, temp_block, m, m, m);
	  res = addToMatrix(out+i*string_size+j*block_size, temp_block, block_size);
	}
	res = simpleMatrixMultiply(a+i*string_size+s*block_size, b+s*string_size+j*small_block_size, temp_block, m, r, m);
	res = addToMatrix(out+i*string_size+j*block_size, temp_block, block_size);
      }
      for (k=0; k<s; k++){
	res = simpleMatrixMultiply(a+i*string_size+k*block_size, b+k*string_size+s*block_size, temp_block, m, m, r);
	res = addToMatrix(out+i*string_size+s*block_size, temp_block, small_block_size);
      }
      res = simpleMatrixMultiply(a+i*string_size+s*block_size, b+s*string_size+s*small_block_size, temp_block, m, r, r);
      res = addToMatrix(out+i*string_size+s*block_size, temp_block, small_block_size);
    }
    for (j=0; j<s; j++){
      for (k=0; k<s; k++){
	res = simpleMatrixMultiply(a+s*string_size+k*small_block_size, b+k*string_size+j*block_size, temp_block, r, m, m);
	res = addToMatrix(out+s*string_size+j*small_block_size, temp_block, small_block_size);
      }
      res = simpleMatrixMultiply(a+s*string_size+s*small_block_size, b+s*string_size+j*small_block_size, temp_block, r, r, m);
      res = addToMatrix(out+s*string_size+j*small_block_size, temp_block, small_block_size);
    }
    for (k=0; k<s; k++){
      res = simpleMatrixMultiply(a+s*string_size+k*small_block_size, b+k*string_size+s*block_size, temp_block, r, m, r);
      res = addToMatrix(out+s*string_size+s*small_block_size, temp_block, smallest_block_size);
    }
    res = simpleMatrixMultiply(a+s*string_size+s*small_block_size, b+s*string_size+s*small_block_size, temp_block, r, r, r);
    res = addToMatrix(out+s*string_size+s*small_block_size, temp_block, smallest_block_size);
  }
  else{
    for (i=0; i<s; i++){
      for (j=0; j<s; j++){
	for (k=0; k<s; k++){
	  res = simpleMatrixMultiply(a+i*string_size+k*block_size, b+k*string_size+j*block_size, temp_block, m, m, m);
	  res = addToMatrix(out+i*string_size+j*block_size, temp_block, block_size);
	}
      }
    }
  }
  delete[] temp_block;
  return 0;
}


double* simpleMatrixAdd(double* a, double* b, int p, int q){//про запас
  int size = p*q;
  double* c = new double[size];
  for (int i=0; i<size; i++){
    c[i]=a[i]+b[i];
  }
  return c;
}

int printUpperLeftBlock(double* a, int n, int m){
  const int q=7;
  double* const temp = new double[q*q];
  if (q<n){
    for (int i=0; i<q; i++){
      for (int j=0; j<q; j++){
	temp[i*q+j]=getIJ(a, n, m, i, j);
      }
    }
    printMatrix(temp, q, q);
  }
  else{
    for (int i=0; i<n; i++){
      for (int j=0; j<n; j++){
	temp[i*n+j]=getIJ(a, n, m, i, j);
      }
    }
    printMatrix(temp, n, n);
  }
  delete[] temp;
  return 0;
}

double getIJ(double* const a, int n, int m, int i, int j){//Fixed: matrix is distributed in another way
  int vert_block = i/m;
  int hor_block = j/m;
  int vert_block_pos = i%m;
  int hor_block_pos = j%m;

  int r=n%m;
  int s=n/m;

  if (r==0){
    return a[vert_block*m*n+hor_block*m*m+vert_block_pos*m+hor_block_pos];
  }
  else{
    if (vert_block!=s){
      if (hor_block!=s){
        return a[vert_block*m*n+hor_block*m*m+vert_block_pos*m+hor_block_pos];
      }
      else{
        return a[vert_block*m*n+hor_block*m*m+vert_block_pos*r+hor_block_pos];
      }
    }
    else{
      if (hor_block!=s){
	       return a[s*m*n+hor_block*m*m+vert_block_pos*m+hor_block_pos];
      }
      else{
	       return a[s*m*n+s*m*m+vert_block_pos*r+hor_block_pos];
      }
    }
  }
}
