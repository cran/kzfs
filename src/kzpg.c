// Ming Luo 2016-2-21

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

SEXP kzpg(SEXP pg, SEXP c, SEXP m)
{
    double pct = REAL(c)[0];
    int N1 = LENGTH(pg);
    int M1 = INTEGER_VALUE(m);
    double *pgp = REAL(pg);

	double sumpg = 0;
	double ccx = 0; 
	
	int ml, mr, lx, rx; 
	int tmpl = 0;
	int tmpr = 0;
	SEXP dl, dr, W, SP;
	
	PROTECT(dl = allocVector(INTSXP, M1+2));
	PROTECT(dr = allocVector(INTSXP, M1+2));
	PROTECT(W = allocVector(INTSXP, N1));
	PROTECT(SP = allocVector(REALSXP, N1));
	
	for (int i = 0; i < N1-1; i++) {
        sumpg += pow(pgp[i+1] - pgp[i], 2);
		INTEGER(W)[i] = 1;
    }
	INTEGER(W)[N1-1] = 1;
	ccx = pct * sumpg;
	
    for (int i = 0; i < N1; i++) {
		for (int k = 2; k <= M1; k++) {
			sumpg = 0;
            lx = MAX(0, (i - k + 1));
            rx = MIN(N1-1, (i + k - 1));
			for (int j = lx; j < rx; j++) {
				sumpg += pow(pgp[j+1] - pgp[j], 2);
			}
			if (sumpg <= ccx) {
                INTEGER(W)[i] = k;
            } else {
				break;
			}
        }
    }
	
    for (int i = 2; i < N1; i++) {
		if ((INTEGER(W)[i-1] > MAX(INTEGER(W)[i-2],INTEGER(W)[i])) & 
			(pgp[i-1] > MAX(pgp[i-2],pgp[i]))) {
				INTEGER(W)[i-1] = MAX(INTEGER(W)[i-1]-2,1);
			}
    }

	for (int i = 0; i < N1; i++) {
 		int ldl = MIN(M1, i); 
		for (int j = 0; j < ldl; j++) {
			tmpl = j + 1;
			INTEGER(dl)[j] = INTEGER(W)[i - j] - INTEGER(W)[i - j - 1];
			if (INTEGER(dl)[j]>=0) { break; }
			if ( j == (ldl - 1) ) { tmpl++; }
		}
        ml = MAX(INTEGER(W)[i], MIN(tmpl, ldl));
		
		int ldr = MIN(M1, N1-i-1);
		for (int j = 0; j < ldr; j++) {
			tmpr = j + 1;
			INTEGER(dr)[j] = INTEGER(W)[i+j+1] - INTEGER(W)[i+j];
			if (INTEGER(dr)[j]<=0) { break; }
			if ( j == (ldr - 1) ) { tmpr++; }
		}
        mr = MAX(INTEGER(W)[i], MIN(tmpr, ldr));
		
        lx = MAX(0, (i - ml + 1));
        rx = MIN(N1-1, (i + mr - 1));
		sumpg = 0;
		double ct = 0;
		for (int q = lx; q <= rx; q++) {
			ct ++;
			sumpg += pgp[q];
		}
        REAL(SP)[i] = sumpg/ct;
    }
	UNPROTECT(4);
    return SP;
}


SEXP kzp2w(SEXP pg, SEXP c, SEXP m)
{
    double pct = REAL(c)[0];
	int lx, rx, N1, N2, M1;
	SEXP dl, dr, W, dim;
	
	dim = GET_DIM(pg);
    N1 = INTEGER(dim)[0];
    N2 = INTEGER(dim)[1];
    M1 = INTEGER_VALUE(m);
    double *pgp = REAL(pg);

	double sumpg = 0;
	double ccx = 0; 
	
	PROTECT(dl = allocVector(INTSXP, M1+2));
	PROTECT(dr = allocVector(INTSXP, M1+2));
	PROTECT(W = allocMatrix(INTSXP, N1, N2));
	
	for (int j = 0; j < N2; j++) {
		for (int i = 0; i < N1-1; i++) {
			sumpg += pow(pgp[(i+1)+j*N1] - pgp[i+j*N1], 2);
			INTEGER(W)[i+j*N1] = 1;
		}
		INTEGER(W)[(N1-1)+j*N1] = 1;
	}
	ccx = pct * sumpg;
	
    for (int i = 0; i < N1; i++) {
		for (int j = 0; j < N2; j++) {
			for (int k = M1; k >= 2; k--) {
				sumpg = 0;
				lx = MAX(0, (i - k + 1));
				rx = MIN(N1-1, (i + k - 1));
				for (int l = lx; l < rx; l++) {
					sumpg += pow(pgp[(l+1)+j*N1] - pgp[l+j*N1], 2);
				}
				if (sumpg <= ccx) {
					INTEGER(W)[i+j*N1] = k;
					break;
				}
			}
		}
	}
	UNPROTECT(3);
    return W;
}
