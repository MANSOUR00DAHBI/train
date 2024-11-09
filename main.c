// main.C 

#include "main.h"
#define LINE printf("------------------------------------------------------\n")
float dt[]={
	0,0,0,
	0,1,1,
    1,0,1,
    1,1,0,
};
int main(void){
    
    srand(time(0));
    float eps  = 1e-1;
    float rate = 1e-1;

    size_t stride = 3;
    size_t n = sizeof(dt)/sizeof(dt[0])/stride;
    Mat ti = { .rows = n, .cols = 2 , .stride = stride, .es = dt};
    Mat to = { .rows = n, .cols = 1 , .stride = stride, .es = dt + 2};
    Mat_xor m = xor_alloc();
    Mat_xor g = xor_alloc();
    mat_file(m.w2,10.8f);
    mat_file(m.w1,20.5f);
    MAT_PRINT(m.a0);
    MAT_PRINT(m.w1);
    MAT_PRINT(m.b1); 
    LINE;
    MAT_PRINT(m.w1);
    MAT_PRINT(m.b1);
    MAT_PRINT(m.a1);  
    LINE;
    for (int i = 0; i < 100; ++i)
    {
    finite_diff(m,g,eps,ti,to);
    learn_xor(m,g,rate);
    }
    printf("cost :%f\n",cost(m,ti,to) );
    MAT_PRINT(m.w2);
    MAT_PRINT(m.b2);
    MAT_PRINT(m.a2);


     //bit = 1 << bit;
   #define bites 1 << bitest ;
    size_t h = 4 << h;

	return 0;
  



}