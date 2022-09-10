#ifndef __MY_ALLOC_H__
#define __MY_ALLOC_H__

#include <iostream>
#include <complex>
using namespace std;

//==========================================================================//
int ALLOC_1D( double*   &p_1Darray, const int& N );
int ALLOC_2D( double**  &p_2Darray, const int& Nx, const int& Ny );
int ALLOC_3D( double*** &p_3Darray, 
	      const int& Nx, 
	      const int& Ny,
	      const int& Nz );

int ALLOC_4D( double****&p_4Darray, 
	      const int& Ni, 
	      const int& Nj,
	      const int& Nk,
	      const int& Nm
	      );

int ALLOC_5D( double*****&p_5Darray, 
	      const int& Ni, 
	      const int& Nj,
	      const int& Nk,
	      const int& Nm,
	      const int& Nn
	      );

int ALLOC_6D( double******&p_6Darray, 
	      const int& Ni, 
	      const int& Nj,
	      const int& Nk,
	      const int& Nm,
	      const int& Nn,
	      const int& Np
	      );


int ALLOC_7D( double *******&p_7Darray, 
	      const int& Ni, 
	      const int& Nj,
	      const int& Nk,
	      const int& Nm,
	      const int& Nn,
	      const int& Np,
	      const int& Nq
	      );

//==========================================================================//

int ALLOC_1D_int( int*   &p_1Darray, const int& N );
int ALLOC_2D_int( int**  &p_2Darray, const int& Nx, const int& Ny );
int ALLOC_3D_int( int*** &p_3Darray, 
		  const int& Nx, 
		  const int& Ny,
		  const int& Nz );

int ALLOC_4D_int( int****&p_4Darray, 
		  const int& Ni, 
		  const int& Nj,
		  const int& Nk,
		  const int& Nm
		  );

int ALLOC_5D_int( int*****&p_5Darray, 
		  const int& Ni, 
		  const int& Nj,
		  const int& Nk,
		  const int& Nm,
		  const int& Nn
		  );

int ALLOC_6D_int( int******&p_6Darray, 
		  const int& Ni, 
		  const int& Nj,
		  const int& Nk,
		  const int& Nm,
		  const int& Nn,
		  const int& Np
		  );


int ALLOC_7D_int( int *******&p_7Darray, 
		  const int& Ni, 
		  const int& Nj,
		  const int& Nk,
		  const int& Nm,
		  const int& Nn,
		  const int& Np,
		  const int& Nq
		  );

//==========================================================================//
int ALLOC_1D_region( double*& p_1Darray, 
		     const int& i_start,  
		     const int& i_end );
int ALLOC_2D_region( double**& p_2Darray, 
		     const int& i_start,  const int& i_end, 
		     const int& j_start,  const int& j_end );
int ALLOC_3D_region( double***& p_3Darray, 
		     const int& i_start,  const int& i_end, 
		     const int& j_start,  const int& j_end, 
		     const int& k_start,  const int& k_end  );
int ALLOC_4D_region( double****& p_4Darray, 
		     const int& i_start,  const int& i_end, 
		     const int& j_start,  const int& j_end, 
		     const int& k_start,  const int& k_end,
		     const int& m_start,  const int& m_end
		     );

int ALLOC_5D_region( double*****& p_5Darray, 
		     const int& i_start,  const int& i_end, 
		     const int& j_start,  const int& j_end, 
		     const int& k_start,  const int& k_end,  
		     const int& m_start,  const int& m_end,
		     const int& n_start,  const int& n_end  
		     );
int ALLOC_6D_region( double******& p_6Darray, 
		     const int& i_start,  const int& i_end, 
		     const int& j_start,  const int& j_end, 
		     const int& k_start,  const int& k_end,  
		     const int& m_start,  const int& m_end,
		     const int& n_start,  const int& n_end,  
		     const int& p_start,  const int& p_end  
		     );

//==========================================================================// 
int ALLOC_1D_int_region( int*& p_1Darray,
			 const int& i_start,
			 const int& i_end );
int ALLOC_2D_int_region( int**& p_2Darray, 
			 const int& i_start,  const int& i_end, 
			 const int& j_start,  const int& j_end );
int ALLOC_3D_int_region( int***& p_3Darray, 
			 const int& i_start,  const int& i_end, 
			 const int& j_start,  const int& j_end, 
			 const int& k_start,  const int& k_end  );
int ALLOC_4D_int_region( int****& p_4Darray, 
			 const int& i_start,  const int& i_end, 
			 const int& j_start,  const int& j_end, 
			 const int& k_start,  const int& k_end,
			 const int& m_start,  const int& m_end
			 );
//==========================================================================//
int ALLOC_1D_FREE( double*& p_1Darray );
int ALLOC_2D_FREE( double**& p_2Darray, const int& Nx );
int ALLOC_3D_FREE( double***& p_3Darray, const int& Nx, const int& Ny );
int ALLOC_4D_FREE( double****& p_4Darray, 
		   const int& Ni, const int& Nj, const int& Nk );
int ALLOC_5D_FREE( double*****& p_5Darray, 
		   const int& Ni, const int& Nj, const int& Nk, const int& Nm);
int ALLOC_6D_FREE( double ******& p_6Darray, 
		   const int& Ni, const int& Nj, const int& Nk,const int& Nm,
		   const int& Nn );
int ALLOC_7D_FREE( double *******& p_7Darray, 
		   const int& Ni, const int& Nj, const int& Nk,const int& Nm,
		   const int& Nn, const int& Np );
//==========================================================================//
int ALLOC_1D_int_FREE( int*& p_1Darray );
int ALLOC_2D_int_FREE( int**& p_2Darray, const int& Nx );
int ALLOC_3D_int_FREE( int***& p_3Darray, const int& Nx, const int& Ny );
int ALLOC_4D_int_FREE( int****& p_4Darray, 
		       const int& Ni, const int& Nj, const int& Nk );
int ALLOC_5D_int_FREE( int*****& p_5Darray, 
		       const int& Ni, const int& Nj, const int& Nk, 
		       const int& Nm  );
int ALLOC_6D_int_FREE( int ******& p_6Darray, 
		       const int& Ni, const int& Nj, const int& Nk, 
		       const int& Nm, const int& Nn );
int ALLOC_7D_int_FREE( int *******& p_7Darray, 
		       const int& Ni, const int& Nj, const int& Nk,
		       const int& Nm, const int& Nn, const int& Np  );
//==========================================================================//
int ALLOC_1D_region_FREE( double*& p_1Darray, 
			  const int& i_start,
			  const int& i_end   );
int ALLOC_2D_region_FREE( double**& p_2Darray, 
			  const int& i_start,  const int& i_end, 
			  const int& j_start,  const int& j_end );
int ALLOC_3D_region_FREE( double***& p_3Darray, 
			  const int& i_start,  const int& i_end, 
			  const int& j_start,  const int& j_end, 
			  const int& k_start,  const int& k_end  );
int ALLOC_4D_region_FREE( double****& p_4Darray, 
			  const int& i_start,  const int& i_end, 
			  const int& j_start,  const int& j_end, 
			  const int& k_start,  const int& k_end,
			  const int& m_start,  const int& m_end
			  );
int ALLOC_5D_region_FREE( double*****& p_5Darray, 
			  const int& i_start,  const int& i_end, 
			  const int& j_start,  const int& j_end, 
			  const int& k_start,  const int& k_end, 
			  const int& m_start,  const int& m_end,   
			  const int& n_start,  const int& n_end  
			  );

//==========================================================================//
//==========================================================================//
// << COMPLEX >> 
//==========================================================================//
//==========================================================================//
int ALLOC_1D_complex( complex<double>* &p_1DC_array, const int& N );
//==========================================================================//
int ALLOC_1D_complex_FREE( complex<double>* &p_1DC_array ); 
//==========================================================================//
int ALLOC_2D_complex( complex<double>** &p_2DC_array, 
		      const int &Nx, 
		      const int &Ny 
		      );
//==========================================================================//
int ALLOC_2D_complex_FREE( complex<double>** &p_1DC_array,
			   const int &Nx, 
			   const int &Ny 
			   );
//==========================================================================//
int ALLOC_3D_complex( complex<double>***& p_2DC_array, 
		      const int& Nx, 
		      const int& Ny, 
		      const int& Nz 
		      );
//==========================================================================//
int ALLOC_3D_complex_FREE( complex<double>**& p_3DC_array,
			   const int& Nx, 
			   const int& Ny 
			   ) ;
//==========================================================================//
#endif // __MY_ALLOC_H__
