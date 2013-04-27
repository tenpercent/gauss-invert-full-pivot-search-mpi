#include <cstdio>

void defineMatrixWithFunction(double* const, int);
int readMatrixFromFile(FILE*, double* const, int);

void printMatrix(double* const, int, int);
int printBlockMatrix(double* const, int, int);

int copyMatrix(double* const, double* const, int);
int swapMatrix(double* const, double* const, int);

int zeroMatrix(double* const, int);
int idMatrix(double* const, int);

int makeBlockMatrix_Rows(double* const, int, int, int, int);
int makeBlockMatrix_Columns(double* const, int, int, int, int);

int makeOrdinaryMatrix(double* const, int, int);

int addToMatrix(double* const, double* const, int);
int subtractFromMatrix(double* const, double* const, int);

double matrixNorm(double* const, int);

void matrInit (double*, int);

int simpleInvert(double* const, double* const, double* const, int);

int simpleMatrixMultiply(double* const, double* const, double* const, int, int, int);
int smartMatrixMultiply(double* const, double* const, double* const, int, int, int);
int blockMatrixMultiply(double*, double*, double*, int, int);//только для квадратных

int checkIfOkay(double*, double*, int);

int calculationDeviation(double* const, double* const, int, int);

int printUpperLeftBlock(double* cons, int, int);

double getIJ(double* const, int, int, int, int);

inline int max(int i, int j){
  return ((i>j)?i:j);
}

void blocksMultiply(double*,int,int,double*,int,int,double*,int,int,int,int);
void blocksMultiplyLast(double*,double*,int,int,double*,int,int,int,int);
void blocksSubtract(double*,int,int,double*,int,int,int,int);
void blocksCopy(double*,int,int,double*,int,int,int,int);
void blocksSwap(double*,int,int,double*,int,int,int,int,int,int);

int sync(int);
