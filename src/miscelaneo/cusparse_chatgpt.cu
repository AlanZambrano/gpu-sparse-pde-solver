#include "common.h" // para usar los "CHECK"
#include <iostream>
#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <cuda.h>

//#define EIGEN_NO_CUDA // para quitar los warnings de Eigen
#include <Eigen/Sparse>



// g++ -I ./path/to/eigen/ main.cpp -o main -std=c++11
int main(int arg, char **argv){
	
	int N = 4;
// Definir la matriz sparse M utilizando Eigen
    //Eigen::SparseMatrix<double, RowMajor> M(3, 3); NO SIRVE EL SEG ARGUMENTO. Para que no falle, dejar su valor por default (ColMajor)
    Eigen::SparseMatrix<double> Mat(N,N);
	
	Mat.insert(0, 0) = 1.0;
    Mat.insert(1, 1) = 1.0;
	Mat.insert(2, 2) = 1.0;
	Mat.insert(3, 3) = 1.0;
	Mat.insert(4, 4) = 1.0;
	Mat.makeCompressed();
	
   // rhs
	/*
    Eigen::SparseMatrix<double> b(N,1);
    b.insert(0, 0) = 1.0;
    b.insert(1, 0) = 2.0;
	b.insert(2, 0) = 3.0;
	b.insert(3, 0) = 4.0;
	b.makeCompressed();
	*/
//	std::cout << b << std::endl;
	
	Eigen::VectorXd b(N);
	b << 1.0, 2.0, 3.0, 4.0;
	//std::cout << b << std::endl;
	double *h_b = b.data();
	
	for(int i=0; i<N; i++){
		std::cout << h_b[i] << std::endl;
	}	

	/*
	for(int nz=0; nz<b.nonZeros(); nz++){
		std::cout << "val_b: " << h_b[nz] << std::endl;
	}
	*/

/*
	std::cout << "M:" << std::endl;
	std::cout << Mat  << std::endl;
*/
	// trasponer a forma RowMajor
	Mat =Mat.transpose();
	Mat.makeCompressed();	

	// triada sparse en host
	double *h_csrValMat    = Mat.valuePtr();	  // valores/datos
	int    *h_csrRowPtrMat = Mat.outerIndexPtr(); // acumulador
	int    *h_csrColIndMat = Mat.innerIndexPtr(); // índice de columna

	int nnz  = Mat.nonZeros();
	int rows = Mat.rows();
	//std::cout << "h_nnz  = " << h_nnz  << std::endl;
	//std::cout << "h_rows = " << h_rows << std::endl;
	/*
	for(int nz=0; nz<h_nnz; nz++){
		std::cout << "val: "	<< h_csrValMat[nz]			<< std::endl;
		std::cout << "colInd: " << h_csrColIndMat[nz]<<"\n" << std::endl;
	}	
	for(int n=0; n<N+1; n++){
		std::cout << "acc: " << h_csrRowPtrMat[n] << std::endl;	
	}
	*/


	// tríada sparse y datos de M en device
	double *d_csrValMat;		
	int    *d_csrRowPtrMat; 
	int    *d_csrColIndMat;

	CHECK(cudaMalloc((void **)&d_csrValMat,    sizeof(double)*   nnz  ));
    CHECK(cudaMalloc((void **)&d_csrRowPtrMat, sizeof(int)   * (N + 1)));
    CHECK(cudaMalloc((void **)&d_csrColIndMat, sizeof(int)   *   nnz  ));

	CHECK(cudaMemcpy(d_csrValMat,    h_csrValMat,    sizeof(double)*nnz, cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(d_csrRowPtrMat, h_csrRowPtrMat, sizeof(int)*(N+1), cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(d_csrColIndMat, h_csrColIndMat, sizeof(int)*nnz, cudaMemcpyHostToDevice));


	// input/output de cusolver
	double *d_b;		// rhs
	double *h_x, *d_x;	// solución	
	double tol  = 1e-5;
    int reorder = 0;	// no-reorder during solving
	int singularity = 0;// cusolver le asigna -1 si M es invertible

	h_x = (double*)malloc( N * sizeof(double) );
	if(h_x == NULL){
		std::cout << "h_x bad alloc" << std::endl;
		exit(1);
	}
	
    CHECK(cudaMalloc((void **)&d_b, sizeof(double) * N ));
    CHECK(cudaMalloc((void **)&d_x, sizeof(double) * N ));

	CHECK(cudaMemcpy(d_b, h_b, sizeof(double)*N, cudaMemcpyHostToDevice));

	/*
	// Create the cuSPARSE handle
    CHECK_CUSPARSE(cusparseCreate(&handle));
	*/

	// create cusolver status y handle
	cusolverStatus_t 
		cusolver_status = CUSOLVER_STATUS_SUCCESS ;
	cusolverSpHandle_t 
		handle;
    cusolver_status = cusolverSpCreate(&handle);	


	// Descriptor
	cusparseMatDescr_t descr;
    CHECK_CUSPARSE(
		cusparseCreateMatDescr(&descr));
    CHECK_CUSPARSE(
		cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL));
    CHECK_CUSPARSE(
		cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO));

	// solver
	cusolver_status = cusolverSpDcsrlsvlu(
		handle, rows, nnz, descr, 
		d_csrValMat, d_csrRowPtrMat, d_csrColIndMat,
		d_b,
		tol, reorder,
		d_x, &singularity); 
	
	if(cusolver_status != CUSOLVER_STATUS_SUCCESS){
		std::cout << "cusolver error: " << cusolver_status << std::endl; 
	};

	CHECK(cudaMemcpy(h_x, d_x, sizeof(double) * N, cudaMemcpyDeviceToHost ));


	/*
	// en algún lado:
	cudaDeviceSynchronize();
	*/

	// Free and Destroy
	cusolverSpDestroy(handle);
	cusparseDestroyMatDescr(descr);
	
	//cudaFree( d_A );




//	1) cudaMalloc para alocar en device M, H y x

//	2) cudaMemcpy (... HostToDevice)
	
//	3) (hace falta crear el descriptor o se setea por default?) 		declaro, creo y seteo el descriptor de la matriz. (no es simétrica, ni antisimétrica debido al término convectivo)

/*	 
//Kolton:
	float *A, *dA;
	int *dANnzPerRow;
	float *dCsrValA;
    int *dCsrRowPtrA;
    int *dCsrColIndA;
    int totalANnz;

	cusparseHandle_t handle = 0;
    cusparseMatDescr_t Adescr = 0;

// Create the cuSPARSE handle
    CHECK_CUSPARSE(cusparseCreate(&handle));

 // Allocate device memory for vectors and the dense form of the matrix A
    CHECK(cudaMalloc((void **)&dA, sizeof(float) * M * N));

// Construct a descriptor of the matrix A
    CHECK_CUSPARSE(cusparseCreateMatDescr(&Adescr));

    CHECK_CUSPARSE(cusparseSetMatType(Adescr, CUSPARSE_MATRIX_TYPE_GENERAL));

    CHECK_CUSPARSE(cusparseSetMatIndexBase(Adescr, CUSPARSE_INDEX_BASE_ZERO));


// Transfer the input vectors and dense matrix A to the device
    CHECK(cudaMemcpy(dA, A, sizeof(float) * M * N, cudaMemcpyHostToDevice));

// (opcional?)Compute the number of non-zero elements in A
    CHECK_CUSPARSE(cusparseSnnz(handle, CUSPARSE_DIRECTION_ROW, M, N, Adescr, dA, M, dANnzPerRow, &totalANnz));

// Allocate device memory to store the sparse CSR representation of A
    CHECK(cudaMalloc((void **)&dCsrValA, sizeof(float) * totalANnz));
    CHECK(cudaMalloc((void **)&dCsrRowPtrA, sizeof(int) * (M + 1)));
    CHECK(cudaMalloc((void **)&dCsrColIndA, sizeof(int) * totalANnz));
*/



	/* //ChatGPT:
	cusparseMatDescr_t descr;
	cusparseCreateMatDescr(&descr);
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

	cusolverSpHandle_t handle;
	cusolverSpCreate(&handle);

	int singularity;
	cusolverSpDcsrlsvlu(handle, rows, nnz, descr, d_M, csrRowPtr.data(), csrColInd.data(), d_b, tolerance, 0, d_x, &singularity);

	cusolverSpDestroy(handle);
	cusparseDestroyMatDescr(descr);
	*/	

//	4) cuSolver "csrlsvlu" (basado en LU con piv. parcial)


return 0;
}











