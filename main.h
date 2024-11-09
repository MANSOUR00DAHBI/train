#include "Mat.h"
 #define _MAIN_H_




float sigmoidf(float x){ return 1.f / (1.f + expf(-x));}
float rand_float (){return (float) rand()/ RAND_MAX ;}
Mat mat_alloc (size_t row , size_t col){
	Mat m ;
	m.rows      =  row  ;
	m.cols      =  col ;
	m.stride    =  m.cols ;
	m.es        = MALLOC(sizeof(*m.es )*row*col);
	ASSERT(m.es != NULL);
	return m ;

}
Mat mat_row(Mat m , size_t row){
	return (Mat){
		.rows   = 1,
		.cols   = m.cols,
		.stride = m.stride,
		.es     = &MAT_ES(m,row,0)
	};
}
void mat_doc(Mat w1 , Mat w2 , Mat b )
{
	ASSERT(w2.es != NULL);
	ASSERT(w1.es != NULL);
	ASSERT(b.es != NULL);
	
	ASSERT(w2.cols == b.rows);
	ASSERT(w1.rows == w2.rows);
	ASSERT(w1.cols == b.cols);
	
	size_t count_cols = w2.cols ;
    for(size_t i = 0; i < w1.rows; ++i){
		for(size_t j = 0; j < w1.cols; ++j){
		    MAT_ES(w1,i,j) = 0;
			for(size_t k = 0 ; k < count_cols ;++k){
	            MAT_ES(w1,i,j) = MAT_ES(w2,i,k) * MAT_ES(b,k,j);
			}
		}
    }
}        
void mat_sum(Mat w1 , Mat b )
{
	ASSERT(w1.rows == b.rows);
	ASSERT(w1.cols == b.cols);

   for (size_t i = 0; i < w1.rows; ++i){
		for (size_t j = 0; j < w1.cols; ++j){
          MAT_ES(w1,i,j) += MAT_ES(b,i,j);
		}
	}
}
void mat_sig(Mat w1){


   for (size_t i = 0; i < w1.rows; ++i){
		for (size_t j = 0; j < w1.cols; ++j){
          MAT_ES(w1,i,j)  = sigmoidf(MAT_ES(w1,i,j));
		}
	}
}
void mat_file(Mat w1 , float x){


   for (size_t i = 0; i < w1.rows; ++i){
		for (size_t j = 0; j < w1.cols; ++j){
          MAT_ES(w1,i,j)  = x;
		}
	}
}
void mat_copy(Mat  dst , Mat src){
	ASSERT(dst.rows == src.rows);
	ASSERT(dst.cols == src.cols);
   for (size_t i = 0; i < dst.rows; ++i){
		for (size_t j = 0; j < dst.cols; ++j){
          MAT_ES(dst,i,j)  = MAT_ES(src,i,j);
		}
	}	
}

void forward(Mat_xor m){

    MAT_ES(m.a0,0,0)= 0;
    MAT_ES(m.a0,0,1)= 1;
    mat_doc(m.a1, m.a0 , m.w1 );
    mat_sum(m.a1,m.b1);
    mat_sig(m.a1);

    mat_doc(m.a2, m.a1 , m.w2 );
    mat_sum(m.a2,m.b2);
    mat_sig(m.a2);  

}

float cost(Mat_xor m , Mat ti , Mat to ){
	ASSERT(ti.rows == to.rows);
	ASSERT(to.cols == m.a2.cols);
	size_t n  = ti.rows;
	size_t q  = to.cols;
	float c ;
	for (size_t i = 0; i < n; ++i){
		Mat x = mat_row(ti ,i);
		Mat y = mat_row(to, i);
		mat_copy(m.a0, x);
		for(size_t j = 0 ; j < q ; ++j){
			float d = MAT_ES(m.a2,0,j) - MAT_ES(y,0,j);
			c      += d * d;
			c      /= n;
		}
	}
	return c ;
}
void finite_diff(Mat_xor m , Mat_xor g , float eps , Mat ti , Mat to){
	float c = cost (m ,ti,to);
	float saved ;
	for(size_t i = 0 ; i < m.w1.rows;++i){
		for(size_t j = 0 ; j < m.w1.cols;++j){
            saved             = MAT_ES(m.w1,i,j);
            MAT_ES(m.w1,i,j) += eps ;
            MAT_ES(g.w1,i,j)  = (cost(m,ti,to) - c )/ eps ;
            MAT_ES(m.w1,i,j)  = saved ; 
		}
	} 
	for(size_t i = 0 ; i < m.b1.rows;++i){
		for(size_t j = 0 ; j < m.b1.cols;++j){
            saved             = MAT_ES(m.b1,i,j);
            MAT_ES(m.b1,i,j) += eps ;
            MAT_ES(g.b1,i,j)  = (cost(m,ti,to) - c )/ eps ;
            MAT_ES(m.b1,i,j)  = saved ; 
		}
	} 
	for(size_t i = 0 ; i < m.a1.rows;++i){
		for(size_t j = 0 ; j < m.a1.cols;++j){
            saved             = MAT_ES(m.a1,i,j);
            MAT_ES(m.a1,i,j) += eps ;
            MAT_ES(g.a1,i,j)  = (cost(m,ti,to) - c )/ eps ;
            MAT_ES(m.a1,i,j)  = saved ; 
		}
	} 		
	for(size_t i = 0 ; i < m.w2.rows;++i){
		for(size_t j = 0 ; j < m.w2.cols;++j){
            saved             = MAT_ES(m.w2,i,j);
            MAT_ES(m.w2,i,j) += eps ;
            MAT_ES(g.w2,i,j)  = (cost(m,ti,to) - c )/ eps ;
            MAT_ES(m.w2,i,j)  = saved ; 
		}
	} 
	for(size_t i = 0 ; i < m.b2.rows;++i){
		for(size_t j = 0 ; j < m.b2.cols;++j){
            saved             = MAT_ES(m.b2,i,j);
            MAT_ES(m.b2,i,j) += eps ;
            MAT_ES(g.b2,i,j)  = (cost(m,ti,to) - c )/ eps ;
            MAT_ES(m.b2,i,j)  = saved ; 
		}
	} 
	for(size_t i = 0 ; i < m.a2.rows;++i){
		for(size_t j = 0 ; j < m.a2.cols;++j){
            saved             = MAT_ES(m.a2,i,j);
            MAT_ES(m.a2,i,j) += eps ;
            MAT_ES(g.a2,i,j)  = (cost(m,ti,to) - c )/ eps ;
            MAT_ES(m.a2,i,j)  = saved ; 
		}
	} 		


}
void learn_xor (Mat_xor m ,Mat_xor g , float rate ){
	
	for(size_t i = 0 ; i < m.w1.rows;++i){
		for(size_t j = 0 ; j < m.w1.cols;++j){
            MAT_ES(m.w1,i,j) -= rate * MAT_ES(g.w1,i,j) ;   
		}
	} 
	for(size_t i = 0 ; i < m.b1.rows;++i){
		for(size_t j = 0 ; j < m.b1.cols;++j){
            MAT_ES(m.b1,i,j) -= rate * MAT_ES(g.b1,i,j) ;   
		}
	}
	for(size_t i = 0 ; i < m.a1.rows;++i){
		for(size_t j = 0 ; j < m.a1.cols;++j){
            MAT_ES(m.a1,i,j) -= rate * MAT_ES(g.a1,i,j) ;   
		}
	}
	for(size_t i = 0 ; i < m.w2.rows;++i){
		for(size_t j = 0 ; j < m.w2.cols;++j){
            MAT_ES(m.w2,i,j) -= rate * MAT_ES(g.w2,i,j) ;   
		}
	}
	for(size_t i = 0 ; i < m.b2.rows;++i){
		for(size_t j = 0 ; j < m.b2.cols;++j){
            MAT_ES(m.b2,i,j) -= rate * MAT_ES(g.b2,i,j) ;   
		}
	}
	for(size_t i = 0 ; i < m.a2.rows;++i){
		for(size_t j = 0 ; j < m.a2.cols;++j){
            MAT_ES(m.a2,i,j) -= rate * MAT_ES(g.a2,i,j) ;   
		}
	}					
}

Mat_xor xor_alloc()
{

	Mat_xor m;
	m.a0    = mat_alloc(1,2);
    m.w1    = mat_alloc(2,2);
    m.b1    = mat_alloc(1,2);
    m.a1    = mat_alloc(1,2);

    m.w2    = mat_alloc(2,1);
    m.b2    = mat_alloc(1,1);
    m.a2    = mat_alloc(1,1);
    return m;
}
void mat_print (Mat m ,const char *name)
{
	printf("%s [\n",name );
	for (size_t i = 0; i < m.rows; ++i){
		for (size_t j = 0; j < m.cols; ++j){
         printf("    |%zu ^ %zu| %f ",i,j,MAT_ES(m,i,j) );
		}
		printf("\n");
	}
	printf("  ]\n");
}