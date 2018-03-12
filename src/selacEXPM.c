#ifdef USE_AVX2

#include <vector>
#include <immintrin.h>

#inlcude "Rcpp.h"
#include "R.h"


//I dont know how to enable AVX2 (2 is important, because AVX only supports 128 byte) in R compilation
//have to learn how to write configure files
void matmult_AVX(std::vector<double> &A, std::vector<double> &B, std::vector<double> &C, unsigned N)
{
    // use jb and kb to find optimal loop tilling
	size_t jb = std::min(512u, N); 
	size_t kb = std::min(24u, N);

//#pragma omp parallel for num_threads(4) schedule(dynamic, 16)
	for (int jj = 0; jj < N; jj += jb)
	{
		for (size_t kk = 0; kk < N; kk += kb)
		{
			for (size_t i = 0; i < N; i += 2) {
				for (size_t j = jj; j < jj + jb; j += 8) {
					__m256d sumA_1, sumB_1, sumA_2, sumB_2;
					if (kk == 0) {
						sumA_1 = sumB_1 = sumA_2 = sumB_2 = _mm256_setzero_pd();
					}
					else {
						sumA_1 = _mm256_load_pd(&C[i * N + j]);
						sumB_1 = _mm256_load_pd(&C[i * N + j + 4]);
						sumA_2 = _mm256_load_pd(&C[(i + 1) * N + j]);
						sumB_2 = _mm256_load_pd(&C[(i + 1) * N + j + 4]);
					}
					size_t limit = std::min(N, (unsigned)(kk + kb));
					for (size_t k = kk; k < limit; k++) {
						__m256d bc_mat1_1 = _mm256_set1_pd(A[i * N + k]);
						__m256d vecA_mat2 = _mm256_loadu_pd(&B[k * N + j]);
						__m256d vecB_mat2 = _mm256_loadu_pd(&B[k * N + j + 4]);
						sumA_1 = _mm256_add_pd(sumA_1, _mm256_mul_pd(bc_mat1_1, vecA_mat2));
						sumB_1 = _mm256_add_pd(sumB_1, _mm256_mul_pd(bc_mat1_1, vecB_mat2));
						__m256d bc_mat1_2 = _mm256_set1_pd(A[(i + 1) * N + k]);
						sumA_2 = _mm256_add_pd(sumA_2, _mm256_mul_pd(bc_mat1_2, vecA_mat2));
						sumB_2 = _mm256_add_pd(sumB_2, _mm256_mul_pd(bc_mat1_2, vecB_mat2));
					}
					_mm256_storeu_pd(&C[i * N + j], sumA_1);
					_mm256_storeu_pd(&C[i * N + j + 4], sumB_1);
					_mm256_storeu_pd(&C[(i + 1) * N + j], sumA_2);
					_mm256_storeu_pd(&C[(i + 1) * N + j + 4], sumB_2);
				}
			}
		}
	}
}

std::vector<double>& MatrixTimesMatrix(std::vector<double> &m1, std::vector<double> &m2)
{
	std::vector<double> out(m1.size(), 0.0);
    matmult_AVX(m1, m2, out, std::sqrt(m1.size()));
	return out;
}

void add(std::vector<double> &out, const std::vector<double> &addition)
{
	unsigned n2 = out.size();
	for (unsigned i = 0; i < n2; i += 4)
	{
		out[i] = out[i] + addition[i];
		out[i + 1] = out[i + 1] + addition[i + 1];
		out[i + 2] = out[i + 2] + addition[i + 2];
		out[i + 3] = out[i + 3] + addition[i + 3];
	}
}

void sub(std::vector<double> &out, const std::vector<double> &substrat)
{
	unsigned n2 = out.size();
	for (unsigned i = 0; i < n2; i += 4)
	{
		out[i] = out[i] - substrat[i];
		out[i + 1] = out[i + 1] - substrat[i + 1];
		out[i + 2] = out[i + 2] - substrat[i + 2];
		out[i + 3] = out[i + 3] - substrat[i + 3];
	}
}

void negativeMatrix(std::vector<double> &m)
{
	unsigned n2 = m.size();
	for (unsigned i = 0; i < n2; i += 4)
	{
		m[i] = -m[i];
		m[i + 1] = -m[i + 1];
		m[i + 2] = -m[i + 2];
		m[i + 3] = -m[i + 3];
	}
}

std::vector<double> I(unsigned n, double s)
{
	std::vector<double> IdentityMatrix(n*n, 0.0);
	for (unsigned i = 0; i < n; i++)
	{
		IdentityMatrix[i * n + i] = s;
	}
	return IdentityMatrix;
}

void MatrixTimesScalar(std::vector<double> &m, double s)
{
	unsigned n2 = m.size();
	for (unsigned i = 0; i < n2; i++)
	{
		m[i] = m[i] * s;
	}
}

void MatrixPlusIndentity(std::vector<double> &m)
{
	unsigned n = sqrt(m.size());
	for (unsigned i = 0; i < n; i++)
	{
		m[i * n + i] += 1.0;
	}
}

// invert matrix without pivot element (less stable but faster, so far unecessary, maybe with HMM matrix)
// add AVX to improve speed
void invert(const std::vector<double>& Aprime, std::vector<double>& R)
{
	unsigned n = sqrt(Aprime.size());
	std::vector<double> L = I(n, 1.0);
	std::vector<double> U(Aprime.size(), 0.0);
	std::vector<double> P = I(n, 1.0);
	// LU Decomposition
	for (unsigned i = 0; i < n; i++)
	{
		for (unsigned j = 0; j < i + 1; j++)
		{
			double s = 0.0;
			for (unsigned k = 0; k < j; k++)
				s += L[j * n + k] * U[k * n + i];
			U[j * n + i] = Aprime[j * n + i] - s;
		}
		for (unsigned j = i; j < n; j++)
		{
			double s = 0.0;
			for (unsigned k = 0; k < i; k++)
				s += L[j * n + k] * U[k * n + i];
			L[j * n + i] = (Aprime[j * n + i] - s) / U[i * n + i];
	
		}
	}

	// Invert Matrix
	std::vector<double> Z(Aprime.size(), 0.0);
	for (unsigned i = 0; i < n; i++) 
	{
		// Find Z (L^-1) with Forward Substitution
		for (unsigned j = 0; j < n; j++) 
		{
			Z[j * n + i] = (i == j ? 1.0 : 0.0);
			for (unsigned k = 0; k < n; k++) 
			{
				if (k != j) 
					Z[j * n + i] -= (L[j * n + k] * Z[k * n + i]);
			}
		}
		// Find X (A^-1) with Backward Substitution
		for (int j = n - 1; j >= 0; j--) 
		{
			R[j * n + i] = Z[j * n + i];
			for (unsigned k = 0; k < n; k++) 
			{
				if (k != j) 
					R[j * n + i] -= (U[j * n + k] * R[k * n + i]);
			}
			R[j * n + i] /= U[j * n + j];
		}
	}
}


// [[Rcpp::export]]
std::vector<double>& expm(std::vector<double>& H, double t, const unsigned p) {
	unsigned n = sqrt(H.size());

	// Calcuate Pade coefficients
	double* c = new double[p + 1u];
	c[0] = 1;
	for (unsigned i = 0, j = 1; i < p; ++i, ++j)
		c[j] = c[i] * ((p - i) / (j * (2.0 * p - i)));

	// Calcuate the infinty norm of H, which is defined as the largest row sum of a matrix
	double norm = 0.0;
	bool all_H_are_zero = true;
	for (unsigned i = 0; i < n; ++i) {
		double temp = 0.0;
		for (unsigned j = 0; j < n; j++) {
			temp += std::abs(H[i * n + j]);
			if (H[i * n + j] != 0.0)
				all_H_are_zero = false;			
		}
		norm = t * std::max<double>(norm, temp);
	}
	// If norm = 0, and all H elements are not NaN or infinity but zero,
	// then U should be identity.
	if (norm == 0.0) {
		if (all_H_are_zero == true || t == 0.0) return I(n, 1.0);
		// Some error happens, H has elements which are NaN or infinity.
		std::cerr << "Null input error in the template expm_pad.\n";
	}


	// Scaling, seek s such that || H*2^(-s) || < 1/2, and set scale = 2^(-s)
	int s = 0;
	double scale = 1.0;
	if (norm > 0.5) {
		s = std::max<int>( 0, static_cast<int>((log(norm) / log(2.0) + 2.0)) );
		scale /= std::pow(2.0, s);
		MatrixTimesScalar(H, (scale * t)); 
	} // else is H = H * 1;

	// Horner evaluation of the irreducible fraction
	// Initialise P (numerator) and Q (denominator)
	std::vector<double> H2 = MatrixTimesMatrix(H, H);
	std::vector<double> Q = I(n, c[p]);
	std::vector<double> P = I(n, c[p - 1]);

	unsigned odd = 1u;
	for (unsigned k = p - 1; k > 0; --k) {
		if (odd == 1) 
		{
			Q = MatrixTimesMatrix(Q, H2);
			add(Q, I(n, c[k - 1]));
		} else {
			P = MatrixTimesMatrix(P, H2);
			add(P, I(n, c[k - 1]));
		}
		odd = 1u - odd;
	}
	(odd == 1) ? (Q = MatrixTimesMatrix(Q, H)) : (P = MatrixTimesMatrix(P, H));
	sub(Q, P);

	delete c;

	//Matrix Inversion
	invert(Q, H2);

	H = MatrixTimesMatrix(H2, P);
	MatrixTimesScalar(H, 2.0);
	MatrixPlusIndentity(H);
	if (odd == 1) negativeMatrix(H);

	// Squaring
	// Spending all my time here, why is this different from H^s, which can be done much faster.
	for (unsigned i = 0; i < s; ++i)
		H = MatrixTimesMatrix(H, H);

	return H;
}
#endif
