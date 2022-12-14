#include "KyMUtil.h"


#define KY_JACOBI_ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
        a[k][l]=h+s*(g-h*tau);
#define KY_JACOBI_MAX_ROTATIONS  20

// Jacobi iteration for the solution of eigenvectors/eigenvalues of a nxn
// real symmetric matrix. Square nxn matrix a; size of matrix in n;
// output eigenvalues in w; and output eigenvectors in v. Resulting
// eigenvalues/vectors are sorted in decreasing order; eigenvectors are
// normalized.
int KyMUtil::Jacobi(double **a, int n, double *w, double **v)
{
	int i, j, k, iq, ip, numPos;
	double tresh, theta, tau, t, sm, s, h, g, c, tmp;
	double bspace[4], zspace[4];
	double *b = bspace;
	double *z = zspace;

	// only allocate memory if the matrix is large
	if (n > 4) {
		b = new double[n];
		z = new double[n]; 
	}

	// initialize
	for (ip = 0; ip < n; ip++) {
		for (iq = 0; iq < n; iq++) {
			v[ip][iq] = 0.0;
		}
		v[ip][ip] = 1.0;
	}
	for (ip = 0; ip < n; ip++) {
		b[ip] = w[ip] = a[ip][ip];
		z[ip] = 0.0;
	}

	// begin rotation sequence
	for (i = 0; i < KY_JACOBI_MAX_ROTATIONS; i++) {
		sm = 0.0;
		for (ip = 0; ip < n-1; ip++) 
			for (iq = ip+1; iq < n; iq++)
				sm += fabs(a[ip][iq]);
		
		if (sm == 0.0)
			break;
		if (i < 3)                                // first 3 sweeps
		    tresh = 0.2*sm/(n*n);
		else
			tresh = 0.0;
	
		for (ip = 0; ip < n-1; ip++) {
			for (iq = ip+1; iq < n; iq++) {
				g = 100.0*fabs(a[ip][iq]);
				// after 4 sweeps
				if (i > 3 && (fabs(w[ip])+g) == fabs(w[ip])
					&& (fabs(w[iq])+g) == fabs(w[iq]))
					a[ip][iq] = 0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h = w[iq] - w[ip];
					if ( (fabs(h)+g) == fabs(h))
						t = (a[ip][iq]) / h;
					else {
						theta = 0.5*h / (a[ip][iq]);
						t = 1.0 / (fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0)
							t = -t;
					}
					c = 1.0 / sqrt(1+t*t);
					s = t*c;
					tau = s/(1.0+c);
					h = t*a[ip][iq];
	  			    z[ip] -= h;
				    z[iq] += h;
				    w[ip] -= h;
				    w[iq] += h;
				    a[ip][iq]=0.0;

	  			    // ip already shifted left by 1 unit
					for (j = 0; j <= ip-1; j++) {
						KY_JACOBI_ROTATE(a, j, ip, j, iq)
					}
					 
				    // ip and iq already shifted left by 1 unit
					for (j = ip+1; j <= iq-1; j++) {
						KY_JACOBI_ROTATE(a, ip, j, j, iq)
					}
			 
					// iq already shifted left by 1 unit
					for (j = iq+1; j < n; j++) {
						KY_JACOBI_ROTATE(a, ip, j, iq, j)
					}
			
					for (j = 0; j < n; j++) {
						KY_JACOBI_ROTATE(v, j, ip, j, iq)
					}
				}
			}
		}
	

		for (ip = 0; ip < n; ip++) {
			b[ip] += z[ip];
			w[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	
	// this is NEVER called
	if (i >= KY_JACOBI_MAX_ROTATIONS ) {
		if (n > 4) {
			delete [] b;
			delete [] z;
		}
		return 0;
	}
	

	// sort eigenfunctions                 these changes do not affect accuracy 
	for (j = 0; j < n-1; j++) {                 // boundary incorrect
		k = j;
		tmp = w[k];
		for (i = j+1; i < n; i++) {               // boundary incorrect, shifted already
			if (w[i] >= tmp) {                  // why exchage if same?
				k = i;
				tmp = w[k];
			}
		}
		if (k != j) {
			w[k] = w[j];
			w[j] = tmp;
			for (i=0; i<n; i++) {
				tmp = v[i][j];
				v[i][j] = v[i][k];
				v[i][k] = tmp;
			}
		}
	}
	// insure eigenvector consistency (i.e., Jacobi can compute vectors that
	// are negative of one another (.707,.707,0) and (-.707,-.707,0). This can
	// reek havoc in hyperstreamline/other stuff. We will select the most
	// positive eigenvector.
	int ceil_half_n = (n >> 1) + (n & 1);
	for (j = 0; j < n; j++) {
		for (numPos = 0, i = 0; i < n; i++)
			if ( v[i][j] >= 0.0 )
				numPos++;
		//    if ( numPos < ceil(double(n)/double(2.0)) )
		if ( numPos < ceil_half_n)
			for(i=0; i<n; i++)
				v[i][j] *= -1.0;
	}

	if (n > 4) {
		delete [] b;
		delete [] z;
	}
	return 1;
}
