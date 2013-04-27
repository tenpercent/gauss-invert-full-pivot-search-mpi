#include "tools.h"
#include <cmath>

int simpleInvert(double* const c, double* const b, double* const a, int m){
  //ищем обратный к блоку размера m*m
  //double eps = 2e-16;//для гильбертовой  
  const double eps = 1e-8;//обычно хватает
  //const double eps = (m<20)?1e-28:1e-8;
  
  if(m<=0){
    //printf("Tryin to invert empty matrix; bad boy\n\t--simpleInvert\n");
    return 1;
  }
  if(m==1){
    if(fabs(c[0])>eps){
      b[0]=1./c[0];
      return 0;
    }
    else{
      return 1;
    }
  }
  const int size = m*m;
  double temp;
  int i, j, k;
  int p=0;
  double tryMain=0.;
  idMatrix(b, m);
  copyMatrix(c, a, size);//забекапились
  for (i=0; i<m; i++){
    //цикл по подматрицам
    temp=0.;
    for (j=i; j<m; j++){
      tryMain=fabs(a[j*m+i]);
      if (tryMain>temp){
	temp=tryMain;
	p=j;
      }
    }
    if (temp<eps){
      return 1;
    }
    for (j=i; j<m; j++){
      //std::swap(a[i*m+j], a[p*m+j]);
      temp = a[i*m+j];
      a[i*m+j] = a[p*m+j];
      a[p*m+j] = temp;
    }
    for (j=0; j<m; j++){
      //std::swap(b[i*m+j], b[p*m+j]);
      temp = b[i*m+j];
      b[i*m+j] = b[p*m+j];
      b[p*m+j] = temp;
    }   
    temp=1./a[i*m+i];
    for (j=i+1; j<m; j++){
      a[i*m+j]*=temp;
    }
    for (j=0; j<m; j++){
      b[i*m+j]*=temp;
    }
    for (j=i+1; j<m; j++){//из j-й строчки вычитаем i-ю      
      temp=a[j*m+i];//если раскомментить, то ничего не будет.
      for (k=0; k<m; k++){
	b[j*m+k]-=b[i*m+k]*temp;
      }
      for (k=i+1; k<m; k++){
	a[j*m+k]-=a[i*m+k]*temp;
      }
    }
  
  /*else{
    return 1;
    }*/
  }
  for (i=m-1; i>0; i--){//строчку с номером i вычитаем из
    for (j=i-1; j>=0; j--){//первых (i-1) строчек
      temp = a[j*m+i];
      for(k=0; k<m; k++){
	b[j*m+k]-=b[i*m+k]*temp;
      }
    }
  }
  return 0;
}
