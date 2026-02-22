#include <iostream>

// definiciones de los métodos "operadores"
#include <Eigen/Sparse>
#include "../include/parametros.h"
#include "../include/estado.h"


// Operadores rhs de las PDEs (2D)
double Estado::F(int i, int j){
	double Fij;

	Fij =	 n_2d_(i*(N+1) + j) * (  1.0 
								   - n_2d_(i*(N+1) + j) 
								  )

			- p  * n_2d_(i*(N+1) + j) 
				 * m_2d_(i*(N+1) + j)
			
			+ e0 * n_2d_(i*(N+1) + j) * (  d2mdx2_2d_(i*(N+1) + j)
								 		 + d2mdy2_2d_(i*(N+1) + j) 
									 	)
			
			+ e0 * (  dndx_2d_(i*(N+1) + j)
			     	* dmdx_2d_(i*(N+1) + j)
			     	+ dndy_2d_(i*(N+1) + j) 
			    	* dmdy_2d_(i*(N+1) + j) 
				   )
			
			+ D1_2d_(i*(N+1) + j) * (  dndx_2d_(i*(N+1) + j) 
								  	 * dndx_2d_(i*(N+1) + j) 
								  	 + dndy_2d_(i*(N+1) + j)
								  	 * dndy_2d_(i*(N+1) + j) 
								 	)
 
			+ D_2d_(i*(N+1) + j)  * (  d2ndx2_2d_(i*(N+1) + j) 
								 	 + d2ndy2_2d_(i*(N+1) + j) 
									);	
	//Fij = 1.0;	
return Fij;

}

////////////////////////////////////////////////////////////////

double Estado::G(int i, int j){
	double Gij;

	Gij =	  p * n_2d_(i*(N+1) + j) 
		  	    * m_2d_(i*(N+1) + j)

			-mu * m_2d_(i*(N+1) + j)

			-chi0 * ( m_2d_(i*(N+1) + j) * (
											  d2ndx2_2d_(i*(N+1) + j)
											+ d2ndy2_2d_(i*(N+1) + j)
										   )
						+ (   dndx_2d_(i*(N+1) + j)
							* dmdx_2d_(i*(N+1) + j)

							+ dndy_2d_(i*(N+1) + j)
							* dmdy_2d_(i*(N+1) + j)
						  )
					)
			+D0 * (   d2mdx2_2d_(i*(N+1) + j)
					+ d2mdy2_2d_(i*(N+1) + j)
				  );
	//Gij = 2.0;
return Gij;

}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


// Elementos del Jacobiano 2D
double Estado::dFdn(int i, int j, int l, int m){
	
	double dFdn_ijlm;
	int iprev, inext, jprev, jnext;

	iprev = (N+i)%(N+1);
	inext = (i+1)%(N+1);

	jprev = (N+j)%(N+1);
	jnext = (j+1)%(N+1);

	// cent
	if( (l==i) && (m==j) ){
		dFdn_ijlm =

		(1.0 - 2.0*n_2d_(i*(N+1) + j) - p*m_2d_(i*(N+1) + j) )
		
		+ e0 * (  d2mdx2_2d_(i*(N+1) + j)
				+ d2mdy2_2d_(i*(N+1) + j)
			   )

		+2.0 * (  (dndx_2d_(i*(N+1) + j))*(dndx_2d_(i*(N+1) + j))
				+ (dndy_2d_(i*(N+1) + j))*(dndy_2d_(i*(N+1) + j))
			   )

		+D1_2d_(i*(N+1) + j) * (  d2ndx2_2d_(i*(N+1) + j) 
								+ d2ndy2_2d_(i*(N+1) + j)
							   )

		-(4.0/(dx*dx)) * D_2d_(i*(N+1) + j);

		//dFdn_ijlm = (3.0 - 1.0)*(-2.0/dt);
	return dFdn_ijlm;
	}	


	// left
	if( (l==i) && (m==jprev) ){
		dFdn_ijlm =
		
		-(e0/(2.0*dx)) * dmdy_2d_(i*(N+1) + j)

		-(1.0/dx) * D1_2d_(i*(N+1) + j) * dndy_2d_(i*(N+1) + j)

		+(1.0/(dx*dx)) * D_2d_(i*(N+1) + j);

		//dFdn_ijlm = 2.0*(-2.0/dt);
	return dFdn_ijlm;
	}


	//right
	if( (l==i) && (m==jnext) ){
		dFdn_ijlm =
	
		(e0/(2.0*dx)) * dmdy_2d_(i*(N+1) + j)

		+ (1.0/dx) * D1_2d_(i*(N+1) + j) * dndy_2d_(i*(N+1) + j)

		+ (1.0/(dx*dx)) * D_2d_(i*(N+1) + j);

		//dFdn_ijlm = 4.0*(-2.0/dt);
	return dFdn_ijlm;
	}

	// up
	if( (l==iprev) && (m==j) ){
		dFdn_ijlm =

		-(e0/(2.0*dx)) * dmdx_2d_(i*(N+1) + j)

		- (1.0/dx) * D1_2d_(i*(N+1) + j) * dndx_2d_(i*(N+1) + j)

		+ (1.0/(dx*dx)) * D_2d_(i*(N+1) + j);
	
		//dFdn_ijlm = 1.0*(-2.0/dt);
	return dFdn_ijlm;
	}

	// down
	if( (l==inext) && (m==j) ){
		dFdn_ijlm = 
	
		(e0/(2.0*dx)) * dmdx_2d_(i*(N+1) + j)

		+ (1.0/dx) * D1_2d_(i*(N+1) + j) * dndx_2d_(i*(N+1) + j)

		+ (1.0/(dx*dx)) * D_2d_(i*(N+1) + j);

		//dFdn_ijlm = 5.0*(-2.0/dt);
	return dFdn_ijlm;
	}
}

////////////////////////////////////////////////////////////////

double Estado::dFdm(int i, int j, int l, int m){
	
	double dFdm_ijlm;
	int iprev, inext, jprev, jnext;

	iprev = (N+i)%(N+1);
	inext = (i+1)%(N+1);

	jprev = (N+j)%(N+1);
	jnext = (j+1)%(N+1);

	// cent
	if( (l==i) && (m==j) ){
		dFdm_ijlm =

		- p * n_2d_(i*(N+1) + j)

		- (4.0*e0/(dx*dx)) * n_2d_(i*(N+1) + j);

		//dFdm_ijlm = 3.0*(-2.0/dt);
	return dFdm_ijlm;
	}	


	// left
	if( (l==i) && (m==jprev) ){
		dFdm_ijlm =
		
		  (e0/(dx*dx)) * n_2d_(i*(N+1) + j)

		- (e0/(2.0*dx)) * dndy_2d_(i*(N+1) + j);

		//dFdm_ijlm = 2.0*(-2.0/dt);
	return dFdm_ijlm;
	}


	//right
	if( (l==i) && (m==jnext) ){
		dFdm_ijlm =
	
 		  (e0/(dx*dx)) * n_2d_(i*(N+1) + j)

		+ (e0/(2.0*dx)) * dndy_2d_(i*(N+1) + j);

		//dFdm_ijlm = 4.0*(-2.0/dt);
	return dFdm_ijlm;
	}

	// up
	if( (l==iprev) && (m==j) ){
		dFdm_ijlm =

		  (e0/(dx*dx)) * n_2d_(i*(N+1) + j)
 
		- (e0/(2.0*dx)) * dndx_2d_(i*(N+1) + j);
	
		//dFdm_ijlm = 1.0*(-2.0/dt);
	return dFdm_ijlm;
	}

	// down
	if( (l==inext) && (m==j) ){
		dFdm_ijlm = 
	
		(e0/(dx*dx)) * n_2d_(i*(N+1) + j)

		+ (e0/(2.0*dx)) * dndx_2d_(i*(N+1) + j);

		//dFdm_ijlm = 5.0*(-2.0/dt);
	return dFdm_ijlm;
	}
}

////////////////////////////////////////////////////////////////

double Estado::dGdn(int i, int j, int l, int m){
	
	double dGdn_ijlm;
	int iprev, inext, jprev, jnext;

	iprev = (N+i)%(N+1);
	inext = (i+1)%(N+1);

	jprev = (N+j)%(N+1);
	jnext = (j+1)%(N+1);

	// cent
	if( (l==i) && (m==j) ){
		dGdn_ijlm =

		  p * m_2d_(i*(N+1) + j)
		
		+ (4.0*chi0/(dx*dx)) * m_2d_(i*(N+1) + j);

		//dGdn_ijlm = 3.0*(-2.0/dt);
	return dGdn_ijlm;
	}	


	// left
	if( (l==i) && (m==jprev) ){
		dGdn_ijlm =
		
		- (chi0/(dx*dx)) * m_2d_(i*(N+1) + j)

		+ (chi0/(2.0*dx)) * dmdy_2d_(i*(N+1) + j);

		//dGdn_ijlm = 2.0*(-2.0/dt);
	return dGdn_ijlm;
	}


	//right
	if( (l==i) && (m==jnext) ){
		dGdn_ijlm =
	
		- (chi0/(dx*dx)) * m_2d_(i*(N+1) + j)

		- (chi0/(2.0*dx)) * dmdy_2d_(i*(N+1) + j);

		//dGdn_ijlm = 4.0*(-2.0/dt);
	return dGdn_ijlm;
	}

	// up
	if( (l==iprev) && (m==j) ){
		dGdn_ijlm =

		- (chi0/(dx*dx)) * m_2d_(i*(N+1) + j)

		+ (chi0/(2.0*dx)) * dmdx_2d_(i*(N+1) + j);
	
		//dGdn_ijlm = 1.0*(-2.0/dt);
	return dGdn_ijlm;
	}

	// down
	if( (l==inext) && (m==j) ){
		dGdn_ijlm = 
	
		- (chi0/(dx*dx)) * m_2d_(i*(N+1) + j)

		- (chi0/(2.0*dx)) * dmdx_2d_(i*(N+1) + j);

		//dGdn_ijlm = 5.0*(-2.0/dt);
	return dGdn_ijlm;
	}
}

////////////////////////////////////////////////////////////////

double Estado::dGdm(int i, int j, int l, int m){
	
	double dGdm_ijlm;
	int iprev, inext, jprev, jnext;

	iprev = (N+i)%(N+1);
	inext = (i+1)%(N+1);

	jprev = (N+j)%(N+1);
	jnext = (j+1)%(N+1);

	// cent
	if( (l==i) && (m==j) ){
		dGdm_ijlm =

		  p * n_2d_(i*(N+1) + j) - mu

		- chi0 * (  d2ndx2_2d_(i*(N+1) + j)
				  + d2ndy2_2d_(i*(N+1) + j)	
				 )

		- (4.0/(dx*dx)) * D0;			  

		//dGdm_ijlm = (3.0 - 1.0)*(-2.0/dt);
	return dGdm_ijlm;
	}	


	// left
	if( (l==i) && (m==jprev) ){
		dGdm_ijlm =
	
		(chi0/(2.0*dx)) * dndy_2d_(i*(N+1) + j)
		
		+ (D0/(dx*dx));		

		//dGdm_ijlm = 2.0*(-2.0/dt);
	return dGdm_ijlm;
	}


	//right
	if( (l==i) && (m==jnext) ){
		dGdm_ijlm =
	
		- (chi0/(2.0*dx)) * dndy_2d_(i*(N+1) + j)

		+ (D0/(dx*dx));

		//dGdm_ijlm = 4.0*(-2.0/dt);
	return dGdm_ijlm;
	}

	// up
	if( (l==iprev) && (m==j) ){
		dGdm_ijlm =

		  (chi0/(2.0*dx)) * dndx_2d_(i*(N+1) + j)
		
		+ (D0/(dx*dx));
	
		//dGdm_ijlm = 1.0*(-2.0/dt);
	return dGdm_ijlm;
	}

	// down
	if( (l==inext) && (m==j) ){
		dGdm_ijlm = 
	
		- (chi0/(2.0*dx)) * dndx_2d_(i*(N+1) + j)

		+ (D0/(dx*dx));

		//dGdm_ijlm = 5.0*(-2.0/dt);
	return dGdm_ijlm;
	}
}

