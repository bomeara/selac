#define USE_OPENCL
#ifdef USE_OPENCL
#include <vector>

#include "CL/cl.hpp"
#include <clblast.h>

#inlcude "Rcpp.h"
#include "R.h"

//I dont know how to enable OpenCL in R compilation have to learn how to write configure files

// OpenCL platform/device settings
const cl::Platform platform_id = 0;
const cl::Platform device_id = 0;
std::vector<cl::Platform> platforms = std::vector<cl::Platform>();
int error_val = cl::Platform::get(&platforms);
cl::Platform platform = platforms[platform_id];
// Initializes the OpenCL device
cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[0])(), 0 };
cl::Context context(CL_DEVICE_TYPE_CPU, properties);
std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
cl::Device device = devices[device_id];
auto queue = cl::CommandQueue(context, device);
auto event = cl_event{ nullptr };


//produce helper functions like add, sub, and so on
const char *kernelAdd =                                       "\n" \
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable                    \n" \
"__kernel void add(  __global double *a,                       \n" \
"                       __global double *b,                       \n" \
"                       const unsigned int n)                    \n" \
"{                                                               \n" \
"    //Get our global thread ID                                  \n" \
"    int id = get_global_id(0);                                  \n" \
"                                                                \n" \
"    //Make sure we do not go out of bounds                      \n" \
"    if (id < n)                                                 \n" \
"        a[id] = a[id] + b[id];                                  \n" \
"}                                                               \n" \
                                                                "\n" ;
// Create the compute program from the source buffer
// Build kernel from source string
cl::Program::Sources source(1, std::make_pair(kernelAdd,strlen(kernelAdd)));
cl::Program programADD = cl::Program(context, source);
programADD.build(devices);
// Create kernel object
cl::Kernel sub(programADD, "add", &err);

const char *kernelSub =                                       "\n" \
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable                    \n" \
"__kernel void sub(  __global double *a,                       \n" \
"                       __global double *b,                       \n" \
"                       __global double *c,                       \n" \
"                       const unsigned int n)                    \n" \
"{                                                               \n" \
"    //Get our global thread ID                                  \n" \
"    int id = get_global_id(0);                                  \n" \
"                                                                \n" \
"    //Make sure we do not go out of bounds                      \n" \
"    if (id < n)                                                 \n" \
"        c[id] = a[id] - b[id];                                  \n" \
"}                                                               \n" \
                                                                "\n" ;                                                               
// Create the compute program from the source buffer
// Build kernel from source string
cl::Program::Sources source(1, std::make_pair(kernelSub,strlen(kernelSub)));
cl::Program programSUB = cl::Program(context, source);
programSUB.build(devices);
// Create kernel object
cl::Kernel sub(programSUB, "sub", &err);

const char *kernelNeg =                                       "\n" \
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable                    \n" \
"__kernel void neg(  __global double *a,                       \n" \
"                       const unsigned int n)                    \n" \
"{                                                               \n" \
"    //Get our global thread ID                                  \n" \
"    int id = get_global_id(0);                                  \n" \
"                                                                \n" \
"    //Make sure we do not go out of bounds                      \n" \
"    if (id < n)                                                 \n" \
"        a[id] = -a[id];                                  \n" \
"}                                                               \n" \
                                                                "\n" ;                                                               
// Create the compute program from the source buffer
// Build kernel from source string
cl::Program::Sources source(1, std::make_pair(kernelNeg,strlen(kernelNeg)));
cl::Program programNEG = cl::Program(context, source);
programNEG.build(devices);
// Create kernel object
cl::Kernel neg(programNEG, "neg", &err);


const char *kernelMtS =                                         "\n" \
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable                   \n" \
"__kernel void mts(  __global double *a,                         \n" \
"                       const double s                           \n" \
"                       const unsigned int n)                    \n" \
"{                                                               \n" \
"    //Get our global thread ID                                  \n" \
"    int id = get_global_id(0);                                  \n" \
"                                                                \n" \
"    //Make sure we do not go out of bounds                      \n" \
"    if (id < n)                                                 \n" \
"        a[id] = a[id] * s;                                      \n" \
"}                                                               \n" \
                                                                "\n" ;                                                               
// Create the compute program from the source buffer
// Build kernel from source string
cl::Program::Sources source(1, std::make_pair(kernelMtS,strlen(kernelMtS)));
cl::Program programMTS = cl::Program(context, source);
programMTS.build(devices);
// Create kernel object
cl::Kernel mts(programMTS, "mts", &err);


// Number of work items in each local work group
cl::NDRange localSize(64);
// Number of total work items - localSize must be devisor
cl::NDRange globalSize((int)(ceil(n/(float)64)*64));

void add(std::vector<double>& m1, const std::vector<double>& m2)
{
    unsigned int n  = m.size();
    auto device_a = cl::Buffer(context, CL_MEM_READ_ONLY, n*sizeof(double));
    auto device_b = cl::Buffer(context, CL_MEM_READ_ONLY, n*sizeof(double));
    queue.enqueueWriteBuffer(device_a, CL_TRUE, 0, n*sizeof(double), m1.data());
    queue.enqueueWriteBuffer(device_b, CL_TRUE, 0, n*sizeof(double), m2.data());

    // Set the arguments to our compute kernel
	add.setArgs(0, device_a);
	add.setArgs(1, device_b);
	add.setArgs(2, n);

    // Execute the kernel over the entire range of the data set
	queue.enqueueNDRangeKernel(add, cl::NullRange, globalSize, localSize, NULL, &event);

    queue.enqueueReadBuffer(device_a, CL_TRUE, 0, n*sizeof(double), &m1[0], NULL, NULL);           
}

void sub(std::vector<double>& m1, const std::vector<double>& m2)
{
    unsigned int n  = m.size();
    auto device_a = cl::Buffer(context, CL_MEM_READ_ONLY, n*sizeof(double));
    auto device_b = cl::Buffer(context, CL_MEM_READ_ONLY, n*sizeof(double));
    queue.enqueueWriteBuffer(device_a, CL_TRUE, 0, n*sizeof(double), m1.data());
    queue.enqueueWriteBuffer(device_b, CL_TRUE, 0, n*sizeof(double), m2.data());

    // Set the arguments to our compute kernel
	sub.setArgs(0, device_a);
	sub.setArgs(1, device_b);
	sub.setArgs(2, n);

    // Execute the kernel over the entire range of the data set
	queue.enqueueNDRangeKernel(sub, cl::NullRange, globalSize, localSize, NULL, &event);

    queue.enqueueReadBuffer(device_a, CL_TRUE, 0, n*sizeof(double), &m1[0], NULL, NULL);           
}

void negate(std::vector<double> m)
{
    unsigned int n  = m.size();
    auto device_a = cl::Buffer(context, CL_MEM_READ_ONLY, n*sizeof(double));
    queue.enqueueWriteBuffer(device_a, CL_TRUE, 0, n*sizeof(double), m.data());

    // Set the arguments to our compute kernel
	neg.setArgs(0, device_a);
	neg.setArgs(1, n);

    // Execute the kernel over the entire range of the data set
	queue.enqueueNDRangeKernel(neg, cl::NullRange, globalSize, localSize, NULL, &event);

    queue.enqueueReadBuffer(device_a, CL_TRUE, 0, n*sizeof(double), &m[0], NULL, NULL);           
}

void MatrixTimesScalar(std::vector<double> m, const double s)
{
    unsigned int n  = m.size();
    auto device_a = cl::Buffer(context, CL_MEM_READ_ONLY, n*sizeof(double));
    queue.enqueueWriteBuffer(device_a, CL_TRUE, 0, n*sizeof(double), m.data());

    // Set the arguments to our compute kernel
	mts.setArgs(0, device_a);
	mts.setArgs(1, s);
	mts.setArgs(2, n);

    // Execute the kernel over the entire range of the data set
	queue.enqueueNDRangeKernel(mts, cl::NullRange, globalSize, localSize, NULL, &event);

    queue.enqueueReadBuffer(device_a, CL_TRUE, 0, n*sizeof(double), &m[0], NULL, NULL);           
}

std::vector<double> MatrixTimesMatrix(std::vector<double> m1, std::vector<double> m2)
{
	std::vector<double> out(m1.size(), 0.0);
	unsigned n = sqrt(m1.size());
	// Copy the matrices to the device
	auto device_a = cl::Buffer(context, CL_MEM_READ_ONLY, m1.size()*sizeof(double));
	auto device_b = cl::Buffer(context, CL_MEM_READ_ONLY, m2.size()*sizeof(double));
	auto device_c = cl::Buffer(context, CL_MEM_WRITE_ONLY, out.size()*sizeof(double));
	queue.enqueueWriteBuffer(device_a, CL_TRUE, 0, m1.size()*sizeof(double), m1.data());
	queue.enqueueWriteBuffer(device_b, CL_TRUE, 0, m2.size()*sizeof(double), m2.data());
	queue.enqueueWriteBuffer(device_c, CL_TRUE, 0, out.size()*sizeof(double), out.data());
	auto queue_plain = queue();
	auto status = clblast::Gemm(clblast::Layout::kRowMajor, clblast::Transpose::kNo, 
		clblast::Transpose::kNo, n, n, n, 1.0, device_a(), 0, n, device_b(), 0, n,
		0.0, device_c(), 0, n, &queue_plain, &event);
	queue.enqueueReadBuffer(device_c, CL_TRUE, 0, n*sizeof(double), &out[0], NULL, NULL);
	return out;
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

void MatrixPlusIndentity(std::vector<double> &m)
{
	unsigned n = sqrt(m.size());
	for (unsigned i = 0; i < n; i++)
	{
		m[i * n + i] += 1.0;
	}
}

// invert matrix without pivot element (less stable but faster, so far unecessary, maybe with HMM matrix)
// add AVX or opencl to improve speed
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
	if (odd == 1) negate(H);

	// Squaring
	// Spending all my time here, why is this different from H^s, which can be done much faster.
	for (unsigned i = 0; i < s; ++i)
		H = MatrixTimesMatrix(H, H);

	return H;
}
#endif