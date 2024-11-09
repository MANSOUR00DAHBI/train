//mat.h 
//strat _MAIN_H_
#ifndef _MAIN_H_
#define _MAIN_H_
#endif 
//End _MAIN_H_

//Start _IMPLEMENTATION_
#ifndef _IMPLEMENTATION_
#define _IMPLEMENTATION_
#define _MAIN_H_
#endif
//End _IMPLEMENTATION_
//Start MALLOC
#ifndef MALLOC
#define MALLOC malloc
#endif//End MALLOC
//Start ASSERT
#ifndef ASSERT
#define ASSERT assert
#endif //End ASSERT
// strat defintion   _MAIN_H_
#ifdef  _MAIN_H_

#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<time.h>
#include<math.h>
#include <assert.h>
#include <stdint.h>
#include <signal.h>
#endif


 #ifdef _IMPLEMENTATION_
typedef float Cell[3];

float sigmoidf(float x);
float rand_float();
typedef struct 
{
  size_t rows ;size_t cols ; size_t stride ;
  float *es ;	
}Mat;
typedef struct 
{
	Mat a0             ;
    Mat w1 , b1  , a1  ;
    Mat w2 , b2  , a2  ;
}Mat_xor;
Mat mat_alloc (size_t row , size_t col);
Mat mat_row(Mat m , size_t row);
Mat rand_mat();
void mat_file(Mat w1 , float x);
void mat_doc(Mat w1 , Mat w2 , Mat b );
void mat_sum(Mat w1 , Mat w2 );
void mat_sig(Mat w1);
void mat_copy(Mat  dst , Mat src);
void forward(Mat_xor);
float cost(Mat_xor m , Mat ti , Mat to );
void finite_diff(Mat_xor m , Mat_xor g , float eps , Mat ti , Mat to);
void learn_xor (Mat_xor m ,Mat_xor g , float rate );
Mat_xor xor_alloc();

void mat_print(Mat m,const char *name);
#define MAT_PRINT(m) mat_print((m),(#m));
#define MAT_ES(m,i,j)(m).es[(i) * (m).stride + (j)] 
 #endif // End  defintion   _MAIN_H_



