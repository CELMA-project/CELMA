/************************************************************************/
/*
  If the USE_SSE2 flag is set, an implementation of the Arakawa procedure
  using SSE2 instructions is used. Otherwise, a normal implementation is
  used. 
 */

#ifdef USE_SSE2  /* If the USE_SSE2 flag is set */
#ifdef __SSE2__   /* If the compiler supports SSE2  */
void Util_3DArakawaNl(double ***resd,double ***ad,double ***bd,double **fac,int nz,int nx,int ny)
{
  const double twelfth = 1.0/12.0;
  int i,j,k;
  static double **a,**b,**res;
  __m128d mm_a=_mm_set1_pd(twelfth);
  __m128d mm_b;
  __m128d a_i0_j0,a_i1_j0,a_i2_j0,a_i0_j1,a_i1_j1,a_i2_j1,a_i0_j2,a_i1_j2,a_i2_j2;
  __m128d b_i0_j0,b_i1_j0,b_i2_j0,b_i0_j1,b_i1_j1,b_i2_j1,b_i0_j2,b_i1_j2,b_i2_j2;
  __m128d tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  __m128d r1, r2, r3, r4;
   
  for(k=0;k<nz;k++)
    {
      a=ad[k];
      b=bd[k];
      res=resd[k];
      
      for(i=0;i<ny;i++) 
	for(j=0;j<nx;j+=2) {
	  a_i0_j0 = _mm_load_pd(a[i-1]+j-1);
	  a_i0_j2 = _mm_load_pd(a[i-1]+j+1);
	  a_i0_j1 = _mm_shuffle_pd(a_i0_j0,a_i0_j2,_MM_SHUFFLE2(0,1));
	  a_i1_j0 = _mm_load_pd(a[i]+j-1);
	  a_i1_j2 = _mm_load_pd(a[i]+j+1);
	  a_i1_j1 = _mm_shuffle_pd(a_i1_j0,a_i1_j2,_MM_SHUFFLE2(0,1));
	  a_i2_j0 = _mm_load_pd(a[i+1]+j-1);
	  a_i2_j2 = _mm_load_pd(a[i+1]+j+1);
	  a_i2_j1 = _mm_shuffle_pd(a_i2_j0,a_i2_j2,_MM_SHUFFLE2(0,1));

	  b_i0_j0 = _mm_load_pd(b[i-1]+j-1);
	  b_i0_j2 = _mm_load_pd(b[i-1]+j+1);
	  /* b_i0_j1 = _mm_shuffle_pd(b_i0_j0,b_i0_j2,_MM_SHUFFLE2(0,1)); */
	  b_i0_j1 = _mm_loadu_pd(b[i-1]+j);
	  b_i1_j0 = _mm_load_pd(b[i]+j-1);
	  b_i1_j2 = _mm_load_pd(b[i]+j+1);
	  /* b_i1_j1 = _mm_shuffle_pd(b_i1_j0,b_i1_j2,_MM_SHUFFLE2(0,1)); */
	  b_i1_j1 = _mm_loadu_pd(b[i  ]+j);
	  b_i2_j0 = _mm_load_pd(b[i+1]+j-1);
	  b_i2_j2 = _mm_load_pd(b[i+1]+j+1);
	  /* b_i2_j1 = _mm_shuffle_pd(b_i2_j0,b_i2_j2,_MM_SHUFFLE2(0,1)); */
	  b_i2_j1 = _mm_loadu_pd(b[i+1]+j);

	  /*
	    If the GCC compiler is used, it is also possible to write the formula
	    using the normal arithmetic operations +, - and *
	    However, the PGI compiler does not understand this, but requires that
	    arithmetic operations on SSE values are implemented with calls to the
	    SSE2 intrinsic functions _mm_add_pd, _mm_sub_pd and _mm_mul_pd.
	  */
	  /*
          mm_b = (
                   ( a_i1_j0 + a_i2_j0 - a_i1_j2 - a_i2_j2 ) * ( b_i2_j1 - b_i1_j1 ) 
                 + ( a_i0_j0 + a_i1_j0 - a_i0_j2 - a_i1_j2 ) * ( b_i1_j1 - b_i0_j1 )
                 + ( a_i2_j1 + a_i2_j2 - a_i0_j1 - a_i0_j2 ) * ( b_i1_j2 - b_i1_j1 )
                 + ( a_i2_j0 + a_i2_j1 - a_i0_j0 - a_i0_j1 ) * ( b_i1_j1 - b_i1_j0 )
  	         + ( a_i2_j1 - a_i1_j2 ) * ( b_i2_j2 - b_i1_j1 )
	         + ( a_i1_j0 - a_i0_j1 ) * ( b_i1_j1 - b_i0_j0 )
	         + ( a_i1_j2 - a_i0_j1 ) * ( b_i0_j2 - b_i1_j1 )
	         + ( a_i2_j1 - a_i1_j0 ) * ( b_i1_j1 - b_i2_j0 )
                 );
	  */
	  /* Compute each row of the above formula into temporary variables */
	  tmp1 = 
	    _mm_mul_pd(
		       _mm_sub_pd(
				  _mm_add_pd(a_i1_j0, a_i2_j0),
				  _mm_add_pd(a_i1_j2, a_i2_j2)
				  ), 
		       _mm_sub_pd(b_i2_j1, b_i1_j1)
		       );
	  /* tmp1 = ( (a_i1_j0 + a_i2_j0) - (a_i1_j2 + a_i2_j2) ) * ( b_i2_j1 - b_i1_j1 ); */
	  tmp2 = _mm_mul_pd(
			    _mm_sub_pd(
				       _mm_add_pd(a_i0_j0, a_i1_j0),
				       _mm_add_pd(a_i0_j2, a_i1_j2)
				       ),
			    _mm_sub_pd(b_i1_j1, b_i0_j1)
			    );
	  /* tmp2 = ( (a_i0_j0 + a_i1_j0) - (a_i0_j2 + a_i1_j2) ) * ( b_i1_j1 - b_i0_j1 ); */
	  tmp3 = _mm_mul_pd(
			    _mm_sub_pd(
				       _mm_add_pd(a_i2_j1, a_i2_j2), 
				       _mm_add_pd(a_i0_j1, a_i0_j2)
				       ),
			    _mm_sub_pd(b_i1_j2, b_i1_j1)
			    );
	  /* tmp3 = ( (a_i2_j1 + a_i2_j2) - (a_i0_j1 + a_i0_j2) ) * ( b_i1_j2 - b_i1_j1 ); */
	  tmp4 = _mm_mul_pd(
			    _mm_sub_pd(
				       _mm_add_pd(a_i2_j0, a_i2_j1),
				       _mm_add_pd(a_i0_j0, a_i0_j1)
				       ),
			    _mm_sub_pd(b_i1_j1, b_i1_j0)
			    );
	  /* tmp4 = ( (a_i2_j0 + a_i2_j1) - (a_i0_j0 + a_i0_j1) ) * ( b_i1_j1 - b_i1_j0 ); */
	  
          tmp5 = _mm_mul_pd(
			    _mm_sub_pd(a_i2_j1, a_i1_j2),
			    _mm_sub_pd(b_i2_j2, b_i1_j1)
			    );
	  /* tmp5 = ( a_i2_j1 - a_i1_j2 ) * ( b_i2_j2 - b_i1_j1 ); */
	  
          tmp6 = _mm_mul_pd(
			    _mm_sub_pd(a_i1_j0, a_i0_j1),
			    _mm_sub_pd(b_i1_j1, b_i0_j0)
			    );
	  /* tmp6 = ( a_i1_j0 - a_i0_j1 ) * ( b_i1_j1 - b_i0_j0 ); */
	  
          tmp7 = _mm_mul_pd(
			    _mm_sub_pd(a_i1_j2, a_i0_j1), 
			    _mm_sub_pd(b_i0_j2, b_i1_j1)
			    );
	  /* tmp7 = ( a_i1_j2 - a_i0_j1 ) * ( b_i0_j2 - b_i1_j1 ); */
	  
          tmp8 = _mm_mul_pd(
			    _mm_sub_pd(a_i2_j1, a_i1_j0), 
			    _mm_sub_pd(b_i1_j1, b_i2_j0)
			    );
	  /* tmp8 = ( a_i2_j1 - a_i1_j0 ) * ( b_i1_j1 - b_i2_j0 ); */
	  
          /* Add all the temporary variables into mm_b */
	  r1 = _mm_add_pd(tmp1, tmp2);
	  r2 = _mm_add_pd(tmp3, tmp4);
	  r3 = _mm_add_pd(tmp5, tmp6);
	  r4 = _mm_add_pd(tmp7, tmp8);
          mm_b = _mm_add_pd(_mm_add_pd(r1, r2), _mm_add_pd(r3, r4));
	  
          /* Store the result into res[i][j] */
          _mm_storeu_pd(
		       res[i]+j, 
		       _mm_mul_pd(_mm_mul_pd(mm_a, mm_b), _mm_loadu_pd(fac+k*nx+j))
		       );
          /* _mm_storeu_pd(res[i]+j, mm_a * mm_b * _mm_loadu_pd(fac+j) ); */

	} /* End for j */
    }  /* End for k */
}
#endif

#else   /* The USE_SSE2 flag is not on, use the normal version of the procedure */

void Util_3DArakawaNl(double ***resd,double ***ad,double ***bd,double **fac,int nz,int nx,int ny)
{
    register  int i,j,k;
    double **a,**b,**res;
    const double  twelfth = 1.0/12.0;
    
    for(k=0;k<nz;k++)
    {
        a=ad[k];
        b=bd[k];
        res=resd[k];
  	
        for(i=0;i<ny;i++) for(j=0;j<nx;j++)  {   
            res[i][j] =
                  (a[i  ][j-1] + a[i+1][j-1]  - a[i  ][j+1]  -  a[i+1][j+1])*(b[i+1][j  ] - b[i  ][j  ])
                + (a[i-1][j-1] + a[i  ][j-1]  - a[i-1][j+1]  -  a[i  ][j+1])*(b[i  ][j  ] - b[i-1][j  ])
                + (a[i+1][j  ] + a[i+1][j+1]  - a[i-1][j  ]  -  a[i-1][j+1])*(b[i  ][j+1] - b[i  ][j  ]) 
                + (a[i+1][j-1] + a[i+1][j  ]  - a[i-1][j-1]  -  a[i-1][j  ])*(b[i  ][j  ] - b[i  ][j-1])
                + (a[i+1][j  ] - a[i  ][j+1])*(b[i+1][j+1] - b[i  ][j  ])
                + (a[i  ][j-1] - a[i-1][j  ])*(b[i  ][j  ] - b[i-1][j-1])
                + (a[i  ][j+1] - a[i-1][j  ])*(b[i-1][j+1] - b[i  ][j  ]) 
                + (a[i+1][j  ] - a[i  ][j-1])*(b[i  ][j  ] - b[i+1][j-1]);
		res[i][j] *=  twelfth*fac[k][j];
        }  
    }

    BUGREPORT;

}

#endif






/***********************************************************************/
/* Implements the convective nonlinearity in the classical form as
   proposed by arakawa in Jour. Comp. Phys No. 1
   to conserve energy and vorticity

   Calculates -{a,b}   */


void Util_ArakawaNl(double **res,double **a,double **b,double *fac,int nx,int ny)
{
    register int ip,ir;

    for(ip=0;ip<ny;ip++) for(ir=0;ir<nx;ir++)   
        res[ip][ir] =
            -1./12.*fac[ir]*(
                ((a[ip][ir+1]-a[ip][ir-1])*(b[ip+1][ir]-b[ip-1][ir])
	             -(a[ip+1][ir]-a[ip-1][ir])*(b[ip][ir+1]-b[ip][ir-1]))
                + 
                ( a[ip  ][ir+1]*( b[ip+1][ir+1] - b[ip-1][ir+1])
                  - a[ip  ][ir-1]*( b[ip+1][ir-1] - b[ip-1][ir-1])
                  - a[ip+1][ir  ]*( b[ip+1][ir+1] - b[ip+1][ir-1])
                  + a[ip-1][ir  ]*( b[ip-1][ir+1] - b[ip-1][ir-1]))
                +
                (a[ip+1][ir+1]*( b[ip+1][ir  ] - b[ip  ][ir+1])
                 - a[ip-1][ir-1]*( b[ip  ][ir-1] - b[ip-1][ir  ])
                 - a[ip+1][ir-1]*( b[ip+1][ir  ] - b[ip  ][ir-1])
                 + a[ip-1][ir+1]*( b[ip  ][ir+1] - b[ip-1][ir  ]))
                );	

}
/************************************************************************/

void Util_Arakawa45(double **res,double **a,double **b,double *fac,int nx,int ny)
{
#ifdef _OPENMP
    static  int i,j;
#pragma omp threadprivate (i,j)
#else
    register int i,j;
#endif


/* Equation 45 of Arakawa  */
/* Evaluates -{a,b}        */
/* As in                   */
/* d_t w = -{phi,w}        */

    BUGREPORT;

    for(i=0;i<ny;i++) for(j=0;j<nx;j++)   
        res[i][j] =
            (a[i  ][j-1] + a[i+1][j-1]  - a[i  ][j+1]  -  a[i+1][j+1])*(b[i+1][j  ] - b[i  ][j  ])
            + 	(a[i-1][j-1] + a[i  ][j-1]  - a[i-1][j+1]  -  a[i  ][j+1])*(b[i  ][j  ] - b[i-1][j  ])
            +   (a[i+1][j  ] + a[i+1][j+1]  - a[i-1][j  ]  -  a[i-1][j+1])*(b[i  ][j+1] - b[i  ][j  ]) 
            + 	(a[i+1][j-1] + a[i+1][j  ]  - a[i-1][j-1]  -  a[i-1][j  ])*(b[i  ][j  ] - b[i  ][j-1]);

    for(i=0;i<ny;i++) for(j=0;j<nx;j++)   
        res[i][j] +=
            (a[i+1][j  ] - a[i  ][j+1])*(b[i+1][j+1] - b[i  ][j  ])
            + (a[i  ][j-1] - a[i-1][j  ])*(b[i  ][j  ] - b[i-1][j-1])
            + (a[i  ][j+1] - a[i-1][j  ])*(b[i-1][j+1] - b[i  ][j  ]) 
            + (a[i+1][j  ] - a[i  ][j-1])*(b[i  ][j  ] - b[i+1][j-1]);


    BUGREPORT;
    for(i=0;i<ny;i++) for(j=0;j<nx;j++)   
        res[i][j] *= 1./12.*fac[j];


}

/************************************************************************/


void Util_Arakawa46(double **res,double **a,double **b,double *fac,int nx,int ny)
{
    register int i,j;

/* Equation 46 of Arakawa */
/* Evaluates -{a,b}        */
    for(i=0;i<ny;i++) for(j=0;j<nx;j++)   
        res[i][j] = 1./12.*fac[j]*(
            (a[i  ][j-1] + a[i+1][j-1]  - a[i  ][j+1]  -  a[i+1][j+1])*(b[i+1][j  ] + b[i  ][j  ])
            - 	(a[i-1][j-1] + a[i  ][j-1]  - a[i-1][j+1]  -  a[i  ][j+1])*(b[i  ][j  ] + b[i-1][j  ]) 
            + 	(a[i+1][j  ] + a[i+1][j+1]  - a[i-1][j  ]  -  a[i-1][j+1])*(b[i  ][j+1] + b[i  ][j  ]) 
            - 	(a[i+1][j-1] + a[i+1][j  ]  - a[i-1][j-1]  -  a[i-1][j  ])*(b[i  ][j  ] + b[i  ][j-1]) 
            +  (a[i+1][j  ] - a[i  ][j+1])*(b[i+1][j+1] + b[i  ][j  ])
            - (a[i  ][j-1] - a[i-1][j  ])*(b[i  ][j  ] + b[i-1][j-1])
            + (a[i  ][j+1] - a[i-1][j  ])*(b[i-1][j+1] + b[i  ][j  ]) 
            - (a[i+1][j  ] - a[i  ][j-1])*(b[i  ][j  ] + b[i+1][j-1]) );
    
}

/************************************************************************/

void Util_3DSet2Zero(double ***res,double ***a,double ***b,double **fac,int nz,int nx,int ny)
{
register int i,j,k;

for(i=0;i<nz;i++)
 for(j=0;j<ny;j++)
   for(k=0;k<nx;k++)
     res[i][j][k] = 0.;
}

