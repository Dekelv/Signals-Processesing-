#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <tgmath.h>

const double pi = 3.141592653589793238462643383279502884;

int *readSignal(int *len) {
  int *x;
  char c;
  scanf("%d:", len);
  x = calloc(*len, sizeof(int));
  do c = getchar(); while (c != '[');
  if (*len > 0) {
    scanf("%d", &x[0]);
    for (int i=1; i < *len; i++) scanf(",%d", &x[i]);
  }
  do c = getchar(); while (c != ']');
  return x;
}

void *safeMalloc(int sz) {
  void *p = calloc(sz, 1);
  if (p == NULL) {
    fprintf(stderr, "Fatal error: safeMalloc(%d) failed.\n", sz);
    exit(EXIT_FAILURE);
  }
  return p;
}

int *makeIntArray(int n) {
  /* allocates dynamic int array of size/length n */
  return safeMalloc(n*sizeof(int));
}

void printSignalc(int l1, int l2, double complex *x) {  
  printf("%d: [", l1);
  
  if (l1 > 0) {
    printf("%d", (int)round(creal(x[0]))/l2);
    for (int i=1; i < l1; i++) printf(",%d", (int)round(creal(x[i]))/l2);
  }
  printf("]\n");
}

double complex *FFT(double complex *a, double complex w, int N){
  double complex *y = calloc(N, sizeof(double complex));
  if(N==1){
    y[0]=a[0];
    return y;
  }

  double complex x=pow(w,0);  

  double complex *aEven, *aOdd;
  
  aEven = calloc(N/2, sizeof(double complex));
  aOdd = calloc(N/2, sizeof(double complex));

  //divide step
  for(int i=0;i<N/2;i++){
      aEven[i] = a[2*i]; 
      aOdd[i] = a[i*2+1];
  }
  
  double complex *yEven = FFT(aEven,pow(w,2),N/2);
  double complex *yOdd = FFT(aOdd,pow(w,2),N/2);

  for(int i=0;i<=N/2-1;i++){
    y[i] = yEven[i]+x*yOdd[i];
    y[i+N/2] = yEven[i]-x*yOdd[i];
    x = x*w;
  }
  return y;
}

int findPowOf2(int n){
  int i=1;
  while(i<n){
    i=i*2;
  }
  return i;
}

double complex *pad(int *s, int l1, int l2){
  
  double complex *x = calloc(l2, sizeof(double complex));

  for(int i=0;i<l2;i++){
      if(i<l1){
        x[i]=s[i];
      }else{
        x[i]=0;
      }
  }
  return x;
}

double complex *multi(double complex *x, double complex *h, int n) {
	double complex *y = calloc(n, sizeof(double complex));
	for(int i = 0; i < n; i++) {
		y[i] = x[i] * h[i];
	}
	return y;
}

int main(int argc, char *argv[]) {
  int lenF, *FIL;
  int N, *X;
  
  FIL = readSignal(&lenF);
  X = readSignal(&N);

  int lenY= N+lenF-1;

  int n2 = findPowOf2(lenY); 

  double complex *x=calloc(n2, sizeof(double complex));

  double complex *f=calloc(n2, sizeof(double complex));

  x = pad(X,N,n2);
  f = pad(FIL,lenF,n2);

  double p = 2*pi/n2;

  double complex omega = exp(I*p);

  double complex *fftx = FFT(x,omega,n2);

  double complex *fftf = FFT(f,omega,n2);

  double complex *conv = multi(fftx,fftf,n2);

  omega = exp(I*p*-1);

  double complex *ifft = FFT(conv,omega,n2);

  printSignalc(lenY, n2, ifft);
  free(FIL);
  free(X);
  free(conv);
  free(ifft);
  free(fftx);
  free(x);
  free(f);
  return 0;
}
