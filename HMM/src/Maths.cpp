/*
 * Maths.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: Marcin Bogusz
 */

#include "Maths.hpp"

namespace EBC
{

Maths::Maths()
{
	srand((unsigned int)time(0));
	z_rndu = rand();

}

double* Maths::matrixMultiply(double *matA, double *matB, int size)
{
	//FIXME - make sure the resulting matrix is zeroed before multiplication!!!
	double* matResult = new double[size*size];
	int i;

	for (int i=0; i<size; i++)
	{
		for (int j=0; j<size; j++)
		{
			matResult[i*size + j]=0;
			for (int k=0; k<size; k++)
			{
				matResult[i*size + j] += matA[i*size+ k] * matB[k *size + j];
			}
		}
	}
	return matResult;
}

void Maths::matrixByDiagonalMultiply(double *matA, double *matDiag, int size)
{
	for (int i=0; i<size; i++)
		for (int j=0; j<size; j++)
			matA[i*size+j] *=  matDiag[j];
}

void Maths::vectorMultiply(double* vector, double factor, int size)
{
	for(int i=0; i<size; i++)
		vector[i] *= factor;
}

void Maths::expLambdaT(double* lambda, double t, int size)
{
	for(int i=0; i<size; i++)
		lambda[i] = exp(lambda[i]*t);
}

double Maths::getRandom(double lo=0.0, double hi=1.0)
{
	double divider = RAND_MAX/hi;
	return lo +  (rand()/divider);
}

void Maths::revMatrixExponentiation(double* Q, double* P)
{

}

int Maths::eigenRealSym(double A[], int n, double Root[], double work[])
{
/* This finds the eigen solution of a real symmetrical matrix A[n*n].  In return,
   A has the right vectors and Root has the eigenvalues.
   work[n] is the working space.
   The matrix is first reduced to a tridiagonal matrix using HouseholderRealSym(),
   and then using the QL algorithm with implicit shifts.

   Adapted from routine tqli in Numerical Recipes in C, with reference to LAPACK
   Ziheng Yang, 23 May 2001
*/
   int status=0;
   HouseholderRealSym(A, n, Root, work);
   status = EigenTridagQLImplicit(Root, work, n, A);
   EigenSort(Root, A, n);

   return(status);
}

int Maths::eigenQREV (double Q[], double pi[], int n,
               double Root[], double U[], double V[], double spacesqrtpi[])
{
/*
   This finds the eigen solution of the rate matrix Q for a time-reversible
   Markov process, using the algorithm for a real symmetric matrix.
   Rate matrix Q = S * diag{pi} = U * diag{Root} * V,
   where S is symmetrical, all elements of pi are positive, and U*V = I.
   space[n] is for storing sqrt(pi).

   [U 0] [Q_0 0] [U^-1 0]    [Root  0]
   [0 I] [0   0] [0    I]  = [0     0]

   Ziheng Yang, 25 December 2001 (ref is CME/eigenQ.pdf)
*/
   int i,j, inew, jnew, nnew, status;
   double *pi_sqrt=spacesqrtpi, small=1e-100;

   for(j=0,nnew=0; j<n; j++)
      if(pi[j]>small)
         pi_sqrt[nnew++] = sqrt(pi[j]);

   /* store in U the symmetrical matrix S = sqrt(D) * Q * sqrt(-D) */

   if(nnew==n) {
      for(i=0; i<n; i++)
         for(j=0,U[i*n+i] = Q[i*n+i]; j<i; j++)
            U[i*n+j] = U[j*n+i] = (Q[i*n+j] * pi_sqrt[i]/pi_sqrt[j]);

      status=eigenRealSym(U, n, Root, V);
      for(i=0; i<n; i++) for(j=0; j<n; j++)  V[i*n+j] = U[j*n+i] * pi_sqrt[j];
      for(i=0; i<n; i++) for(j=0; j<n; j++)  U[i*n+j] /= pi_sqrt[i];
   }
   else {
      for(i=0,inew=0; i<n; i++) {
         if(pi[i]>small) {
            for(j=0,jnew=0; j<i; j++)
               if(pi[j]>small) {
                  U[inew*nnew+jnew] = U[jnew*nnew+inew]
                                    = Q[i*n+j] * pi_sqrt[inew]/pi_sqrt[jnew];
                  jnew++;
               }
            U[inew*nnew+inew] = Q[i*n+i];
            inew++;
         }
      }

      status=eigenRealSym(U, nnew, Root, V);

      for(i=n-1,inew=nnew-1; i>=0; i--)   /* construct Root */
         Root[i] = (pi[i]>small ? Root[inew--] : 0);
      for(i=n-1,inew=nnew-1; i>=0; i--) {  /* construct V */
         if(pi[i]>small) {
            for(j=n-1,jnew=nnew-1; j>=0; j--)
               if(pi[j]>small) {
                  V[i*n+j] = U[jnew*nnew+inew]*pi_sqrt[jnew];
                  jnew--;
               }
               else
                  V[i*n+j] = (i==j);
            inew--;
         }
         else
            for(j=0; j<n; j++)  V[i*n+j] = (i==j);
      }
      for(i=n-1,inew=nnew-1; i>=0; i--) {  /* construct U */
         if(pi[i]>small) {
            for(j=n-1,jnew=nnew-1;j>=0;j--)
               if(pi[j]>small) {
                  U[i*n+j] = U[inew*nnew+jnew]/pi_sqrt[inew];
                  jnew--;
               }
               else
                  U[i*n+j] = (i==j);
            inew--;
         }
         else
            for(j=0;j<n;j++)
               U[i*n+j] = (i==j);
      }
   }
}

void Maths::HouseholderRealSym(double a[], int n, double d[], double e[])
{
/* This uses HouseholderRealSym transformation to reduce a real symmetrical matrix
   a[n*n] into a tridiagonal matrix represented by d and e.
   d[] is the diagonal (eigends), and e[] the off-diagonal.
*/
   int m,k,j,i;
   double scale,hh,h,g,f;

   for (i=n-1;i>=1;i--) {
      m=i-1;
      h=scale=0;
      if (m > 0) {
         for (k=0;k<=m;k++)
            scale += fabs(a[i*n+k]);
         if (scale == 0)
            e[i]=a[i*n+m];
         else {
            for (k=0;k<=m;k++) {
               a[i*n+k] /= scale;
               h += a[i*n+k]*a[i*n+k];
            }
            f=a[i*n+m];
            g=(f >= 0 ? -sqrt(h) : sqrt(h));
            e[i]=scale*g;
            h -= f*g;
            a[i*n+m]=f-g;
            f=0;
            for (j=0;j<=m;j++) {
               a[j*n+i]=a[i*n+j]/h;
               g=0;
               for (k=0;k<=j;k++)
                  g += a[j*n+k]*a[i*n+k];
               for (k=j+1;k<=m;k++)
                  g += a[k*n+j]*a[i*n+k];
               e[j]=g/h;
               f += e[j]*a[i*n+j];
            }
            hh=f/(h*2);
            for (j=0;j<=m;j++) {
               f=a[i*n+j];
               e[j]=g=e[j]-hh*f;
               for (k=0;k<=j;k++)
                  a[j*n+k] -= (f*e[k]+g*a[i*n+k]);
            }
         }
      }
      else
         e[i]=a[i*n+m];
      d[i]=h;
   }
   d[0]=e[0]=0;

   /* Get eigenvectors */
   for (i=0;i<n;i++) {
      m=i-1;
      if (d[i]) {
         for (j=0;j<=m;j++) {
            g=0;
            for (k=0;k<=m;k++)
               g += a[i*n+k]*a[k*n+j];
            for (k=0;k<=m;k++)
               a[k*n+j] -= g*a[k*n+i];
         }
      }
      d[i]=a[i*n+i];
      a[i*n+i]=1;
      for (j=0;j<=m;j++) a[j*n+i]=a[i*n+j]=0;
   }
}

void Maths::EigenSort(double d[], double U[], int n)
{
/* this sorts the eigen values d[] and rearrange the (right) eigen vectors U[]
*/
   int k,j,i;
   double p;

   for (i=0; i<n-1; i++) {
      p = d[k=i];
      for (j=i+1; j<n; j++)
         if (d[j] >= p) p = d[k=j];
      if (k != i) {
         d[k] = d[i];
         d[i] = p;
         for (j=0;j<n;j++) {
            p = U[j*n+i];
            U[j*n+i] = U[j*n+k];
            U[j*n+k] = p;
         }
      }
   }
}

int Maths::EigenTridagQLImplicit(double d[], double e[], int n, double z[])
{
/* This finds the eigen solution of a tridiagonal matrix represented by d and e.
   d[] is the diagonal (eigenvalues), e[] is the off-diagonal
   z[n*n]: as input should have the identity matrix to get the eigen solution of the
   tridiagonal matrix, or the output from HouseholderRealSym() to get the
   eigen solution to the original real symmetric matrix.
   z[n*n]: has the orthogonal matrix as output

   Adapted from routine tqli in Numerical Recipes in C, with reference to
   LAPACK fortran code.
   Ziheng Yang, May 2001
*/
   int m,j,iter,niter=30, status=0, i,k;
   double s,r,p,g,f,dd,c,b, aa,bb;

   for (i=1;i<n;i++) e[i-1]=e[i];  e[n-1]=0;
   for (j=0;j<n;j++) {
      iter=0;
      do {
         for (m=j;m<n-1;m++) {
            dd=fabs(d[m])+fabs(d[m+1]);
            if (fabs(e[m])+dd == dd) break;  /* ??? */
         }
         if (m != j) {
            if (iter++ == niter) {
               status=-1;
               break;
            }
            g=(d[j+1]-d[j])/(2*e[j]);

            /* r=pythag(g,1); */

            if((aa=fabs(g))>1)  r=aa*sqrt(1+1/(g*g));
            else                r=sqrt(1+g*g);

            g=d[m]-d[j]+e[j]/(g+((g) >= 0.0 ? fabs(r) : -fabs(r)));
            s=c=1;
            p=0;
            for (i=m-1;i>=j;i--) {
               f=s*e[i];
               b=c*e[i];

               /*  r=pythag(f,g);  */
               aa=fabs(f); bb=fabs(g);
               if(aa>bb)       { bb/=aa;  r = aa*sqrt(1+bb*bb); }
               else if(bb==0)             r = 0;
               else            { aa/=bb;  r = bb*sqrt(1+aa*aa); }

               e[i+1]=r;
               if (r == 0) {
                  d[i+1] -= p;
                  e[m]=0;
                  break;
               }
               s=f/r;
               c=g/r;
               g=d[i+1]-p;
               r=(d[i]-g)*s+2*c*b;
               d[i+1]=g+(p=s*r);
               g=c*r-b;
               for (k=0;k<n;k++) {
                  f=z[k*n+i+1];
                  z[k*n+i+1]=s*z[k*n+i]+c*f;
                  z[k*n+i]=c*z[k*n+i]-s*f;
               }
            }
            if (r == 0 && i >= j) continue;
            d[j]-=p; e[j]=g; e[m]=0;
         }
      } while (m != j);
   }
   return(status);
}

/*From PAML by Ziheng Yang*/
double Maths::rndu (void)
{
/* 32-bit integer assumed.
   From Ripley (1987) p. 46 or table 2.4 line 2.
   This may return 0 or 1, which can be a problem.
*/
   z_rndu = z_rndu*69069 + 1;
   if(z_rndu==0 || z_rndu==4294967295)  z_rndu = 13;
   return z_rndu/4294967295.0;
}



} /* namespace EBC */
