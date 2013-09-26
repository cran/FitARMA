/* GetSimARMA(z,a,beta,par) computes the remainder of the simulation of an ARMA given
*  initial values in z and a as well as ARMA parameters in
*  beta = (phi(1),...,phi(p),theta(1),...,theta(q)).  This is used by the R function
*  SimulateGaussianARMA to speed up the computation.
*
*CAUTION: you must dimension z and a exactly right. That is z and a are vectors
*         whose length must be z[n] and a[n+q], where n=par[1] and q=par[3]. Also
*         note that par is passed as a double vector -- this seems more reliable
*         than an integer vector.
*
*  AIM, June 25, 2009
*/

#include <R.h>

void GetSimARMA(double *z, double *a, double *beta, double *par)

{
/* Input parameters:
*      z vector of length n.  z[1:p] contains intial values
*      a vector of length n+q containing Gaussian white noise
*      beta vector of length p+q containing the phi and theta parameters
*      par vector of length 3 containing n, p, q
* NB:since passing int or int * doesn't seem to work, par is passed as double
*
* Output parameters
*   z vector of length containing the simulated time series
*/
	int i, j, q, model, n, p;
			
	n = (int)(par[0]);
	p = (int)(par[1]);
	q = (int)(par[2]);
		
	if (p>0 && q>0)
		model= 3; /* mixed case */
	else if (p>0 && q==0)
		model = 2; /* ar */
	else if (p==0 && q>0) /* ma */
		model = 1;
	else
		model = 0; /* white noise */

	switch (model) {
	case 3: /* mixed */
		for (i=p; i<n; i++)
		{
			z[i] = a[i+q];
			for (j=0; j<p; j++)
				z[i] = z[i]+beta[j]*z[i-1-j];
			for (j=0; j<q; j++)
				z[i] = z[i]-beta[j+p]*a[i-1-j+q];
		}
		break;
	case 2:  /* pure ar */
		for (i=p; i<n; i++)
		{
			z[i] = a[i];
			for (j=0; j<p; j++)
				z[i] = z[i]+beta[j]*z[i-1-j];
		}
		break;
	case 1: /* pure ma */
		for (i=0; i<n; i++)
		{
			z[i] = a[i+q];
			for (j=0; j<q; j++)
				z[i] = z[i]-beta[j]*a[i-1-j+q];
		}
		break;
	case 0:  /* just white moise */
		for (i=0; i<n; i++)
		{
			z[i] = a[i];
		}
		break;
	}
	return;
	}
