/* Copyright (c) 2021, National University of Defense Technology. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef YHAMG_H
#define YHAMG_H

#include <mpi.h>
#include <stdint.h>
#include <cmath>

namespace YHAMG
{

struct complex
{
	float real;
	float imag;

	complex(float re = 0.0, float im = 0.0)
		: real(re),
		imag(im)
	{
	}

	complex(const complex &a)
		: real(a.real),
		imag(a.imag)
	{
	}

	complex operator+() const
	{
		return *this;
	}

	complex operator-() const
	{
		return complex(-real, -imag);
	}

	complex operator+(const float& a) const
	{
		return complex(real + a, imag);
	}

	complex operator-(const float& a) const
	{
		return complex(real - a, imag);
	}

	complex operator*(const float& a) const
	{
		return complex(real * a, imag * a);
	}

	complex operator/(const float& a) const
	{
		return complex(real / a, imag / a);
	}

	complex operator+(const complex& a) const
	{
		return complex(real + a.real, imag + a.imag);
	}

	complex operator-(const complex& a) const
	{
		return complex(real - a.real, imag - a.imag);
	}

	complex operator*(const complex& a) const
	{
		return complex(real * a.real - imag * a.imag, real * a.imag + imag * a.real);
	}

	complex operator/(const complex& a) const
	{
		float denom = 1.0 / (a.real * a.real + a.imag * a.imag);
		return complex((real * a.real + imag * a.imag) * denom, (-real * a.imag + imag * a.real) * denom);
	}

	complex& operator+=(const float& a)
	{
		real += a;
		return *this;
	}

	complex& operator-=(const float& a)
	{
		real -= a;
		return *this;
	}

	complex& operator*=(const float& a)
	{
		real *= a;
		imag *= a;
		return *this;
	}

	complex& operator/=(const float& a)
	{
		real /= a;
		imag /= a;
		return *this;
	}
	
	complex& operator+=(const complex& a)
	{
		real += a.real;
		imag += a.imag;
		return *this;
	}

	complex& operator-=(const complex& a)
	{
		real -= a.real;
		imag -= a.imag;
		return *this;
	}

	complex& operator*=(const complex& a)
	{
		float temp = real * a.real - imag * a.imag;
		imag = real * a.imag + imag * a.real;
		real = temp;
		return *this;
	}

	complex& operator/=(const complex& a)
	{
		float denom = 1.0 / (a.real * a.real + a.imag * a.imag);
		float temp = (real * a.real + imag * a.imag) * denom;
		imag = (-real * a.imag + imag * a.real) * denom;
		real = temp;
		return *this;
	}

	bool operator==(const float& a) const
	{
		return real == a && imag == 0;
	}

	bool operator!=(const float& a) const
	{
		return real != a || imag != 0;
	}

	bool operator==(const complex& a) const
	{
		return real == a.real && imag == a.imag;
	}

	bool operator!=(const complex& a) const
	{
		return real != a.real || imag != a.imag;
	}

	complex& operator=(const float& a)
	{
		real = a;
		imag = 0.0;
		return *this;
	}

	complex& operator=(const complex& a)
	{
		real = a.real;
		imag = a.imag;
		return *this;
	}
};

static inline float creal(const complex a)
{
	return a.real;
}

static inline float cimag(const complex a)
{
	return a.imag;
}

static inline float cabs(const complex a)
{
	return sqrt(a.real * a.real + a.imag * a.imag);
}

static inline complex cconj(const complex a)
{
	return complex(a.real, -a.imag);
}

static inline complex csqrt(const complex a)
{
	if (a.real > 0)
	{
		float temp = sqrt(0.5 * (cabs(a) + a.real));
		return complex(temp, 0.5 * a.imag / temp);
	}
	else
	{
		float temp = sqrt(0.5 * (cabs(a) - a.real));
		if (a.imag >= 0)
			return complex(0.5 * a.imag / temp, temp);
		return complex(-0.5 * a.imag / temp, -temp);
	}
}

struct zomplex
{
	double real;
	double imag;

	zomplex(double re = 0.0, double im = 0.0)
		: real(re),
		imag(im)
	{
	}

	zomplex(const zomplex &a)
		: real(a.real),
		imag(a.imag)
	{
	}

	zomplex operator+() const
	{
		return *this;
	}

	zomplex operator-() const
	{
		return zomplex(-real, -imag);
	}

	zomplex operator+(const double& a) const
	{
		return zomplex(real + a, imag);
	}

	zomplex operator-(const double& a) const
	{
		return zomplex(real - a, imag);
	}

	zomplex operator*(const double& a) const
	{
		return zomplex(real * a, imag * a);
	}

	zomplex operator/(const double& a) const
	{
		return zomplex(real / a, imag / a);
	}

	zomplex operator+(const zomplex& a) const
	{
		return zomplex(real + a.real, imag + a.imag);
	}

	zomplex operator-(const zomplex& a) const
	{
		return zomplex(real - a.real, imag - a.imag);
	}

	zomplex operator*(const zomplex& a) const
	{
		return zomplex(real * a.real - imag * a.imag, real * a.imag + imag * a.real);
	}

	zomplex operator/(const zomplex& a) const
	{
		double denom = 1.0 / (a.real * a.real + a.imag * a.imag);
		return zomplex((real * a.real + imag * a.imag) * denom, (-real * a.imag + imag * a.real) * denom);
	}

	zomplex& operator+=(const double& a)
	{
		real += a;
		return *this;
	}

	zomplex& operator-=(const double& a)
	{
		real -= a;
		return *this;
	}

	zomplex& operator*=(const double& a)
	{
		real *= a;
		imag *= a;
		return *this;
	}

	zomplex& operator/=(const double& a)
	{
		real /= a;
		imag /= a;
		return *this;
	}
	
	zomplex& operator+=(const zomplex& a)
	{
		real += a.real;
		imag += a.imag;
		return *this;
	}

	zomplex& operator-=(const zomplex& a)
	{
		real -= a.real;
		imag -= a.imag;
		return *this;
	}

	zomplex& operator*=(const zomplex& a)
	{
		double temp = real * a.real - imag * a.imag;
		imag = real * a.imag + imag * a.real;
		real = temp;
		return *this;
	}

	zomplex& operator/=(const zomplex& a)
	{
		double denom = 1.0 / (a.real * a.real + a.imag * a.imag);
		double temp = (real * a.real + imag * a.imag) * denom;
		imag = (-real * a.imag + imag * a.real) * denom;
		real = temp;
		return *this;
	}

	bool operator==(const double& a) const
	{
		return real == a && imag == 0;
	}

	bool operator!=(const double& a) const
	{
		return real != a || imag != 0;
	}

	bool operator==(const zomplex& a) const
	{
		return real == a.real && imag == a.imag;
	}

	bool operator!=(const zomplex& a) const
	{
		return real != a.real || imag != a.imag;
	}

	zomplex& operator=(const double& a)
	{
		real = a;
		imag = 0.0;
		return *this;
	}

	zomplex& operator=(const zomplex& a)
	{
		real = a.real;
		imag = a.imag;
		return *this;
	}
};

static inline double zreal(const zomplex a)
{
	return a.real;
}

static inline double zimag(const zomplex a)
{
	return a.imag;
}

static inline double zabs(const zomplex a)
{
	return sqrt(a.real * a.real + a.imag * a.imag);
}

static inline zomplex zconj(const zomplex a)
{
	return zomplex(a.real, -a.imag);
}

static inline zomplex zsqrt(const zomplex a)
{
	if (a.real > 0)
	{
		double temp = sqrt(0.5 * (zabs(a) + a.real));
		return zomplex(temp, 0.5 * a.imag / temp);
	}
	else
	{
		double temp = sqrt(0.5 * (zabs(a) - a.real));
		if (a.imag >= 0)
			return zomplex(0.5 * a.imag / temp, temp);
		return zomplex(-0.5 * a.imag / temp, -temp);
	}
}

struct global
{
	unsigned int indl;
	unsigned int from;
	
	global(unsigned int _indl, unsigned int _from = 0)
	: indl(_indl),
	from(_from)
	{
	}

	global(uint64_t indg = 0)
	: indl(indg & 0xffffffff),
	from(indg >> 32)
	{
	}
	
	operator uint64_t() const
	{
		return (uint64_t)from << 32 | (uint64_t)indl;
	}
	
	bool operator<(const global& a) const
	{
		return from == a.from ? indl < a.indl : from < a.from;
	}
};

static inline unsigned int local(global indg)
{
	return indg.indl;
}

static inline unsigned int owner(global indg)
{
	return indg.from;
}

struct Vector
{
	int ref;
	int size;
	double* values;

	Vector();
	Vector(int n);
	Vector(int n, double* values, int ref);
	Vector(const Vector& x);
	Vector(Vector&& x); 
	~Vector();
	Vector& operator=(double a);
	Vector& operator=(const Vector& x);
	Vector& operator=(Vector&& x);
	double& operator[](int i) const;

	void Free();
	void Resize(int n);
	void Fill(double a) const;
	void FillRandom() const;
	void Copy(const Vector& x) const;
	void Scale(double a) const;
	void AddScaled(double a, const Vector& x) const;
	void Add2Scaled(double a, const Vector& x, double b, const Vector& y) const;
	void Refer(const Vector& x);
};

void VecRead(const char* filename, Vector& x);
void VecWrite(const char* filename, const Vector& x);
void VecAXPBY(double alpha, const Vector& x, double beta, const Vector& y);
void VecAXPBYPCZ(double alpha, const Vector& x, double beta, const Vector& y, double gamma, const Vector& z);
double VecDot(const Vector& x, const Vector& y);
void VecElemMul(const Vector& x, const Vector& y);
void VecRecip(const Vector& x);

struct ParVector
{
	MPI_Comm comm;
	Vector local;

	ParVector(MPI_Comm comm = MPI_COMM_SELF);
	ParVector(MPI_Comm comm, int local_size);
	ParVector(MPI_Comm comm, int local_size, double* local_values, int local_ref);
	ParVector(const Vector& x);
	ParVector(const ParVector& x);
	ParVector(ParVector&& x);
	ParVector& operator=(double a);
	ParVector& operator=(const ParVector& x);
	ParVector& operator=(ParVector&& x);
	double& operator[](int i) const;

	void Free();
	void Resize(int n);
	void Refer(const ParVector& x);
	void Fill(double a) const;
	void FillRandom() const;
	void Copy(const ParVector& x) const;
	void Scale(double a) const;
	void AddScaled(double a, const ParVector& x) const;
	void Add2Scaled(double a, const ParVector& x, double b, const ParVector& y) const;
};

void ParVecAXPBY(double alpha, const ParVector& x, double beta, const ParVector& y);
void ParVecAXPBYPCZ(double alpha, const ParVector& x, double beta, const ParVector& y, double gamma, const ParVector& z);
double ParVecDot(const ParVector& x, const ParVector& y);
void ParVecElemMul(const ParVector& x, const ParVector& y);
void ParVecRecip(const ParVector& x);

struct MultiVector
{
	int ref;
	int nvec;
	int size;
	double* values;

	MultiVector();
	MultiVector(int m, int n);
	MultiVector(int m, int n, double* values, int ref);
	MultiVector(const MultiVector& X);
	MultiVector(MultiVector&& X);
	~MultiVector();
	MultiVector& operator=(const MultiVector& X);
	MultiVector& operator=(MultiVector&& X);

	void Free();
	void Allocate(int m, int n);
	void Refer(const MultiVector& X);
	const Vector operator()(int j) const;
};

struct ParMultiVector
{
	MPI_Comm comm;
	MultiVector local;

	ParMultiVector(MPI_Comm comm = MPI_COMM_SELF);
	ParMultiVector(MPI_Comm comm, int nvec, int local_size);
	ParMultiVector(MPI_Comm comm, int nvec, int local_size, double* local_values, int local_ref);
	ParMultiVector(const MultiVector& X);
	ParMultiVector(const ParMultiVector& X);
	ParMultiVector(ParMultiVector&& X);
	ParMultiVector& operator=(const ParMultiVector& X);
	ParMultiVector& operator=(ParMultiVector&& X);

	void Free();
	void Allocate(int m, int n);
	void Refer(const ParMultiVector& X);
	const ParVector operator()(int j) const;
};

class Operator
{
public:
	virtual ~Operator();
	virtual int InSize() const = 0;
	virtual int OutSize() const = 0;
	virtual void Apply(const Vector& x, const Vector& y) const = 0;
	virtual void Apply(const MultiVector& X, const MultiVector& Y) const;
};

class ParOperator
{
public:
	MPI_Comm comm;

	ParOperator(MPI_Comm comm = MPI_COMM_SELF);
	virtual ~ParOperator();
	virtual int InSize() const = 0;
	virtual int OutSize() const = 0;
	virtual void Apply(const ParVector& x, const ParVector& y) const = 0;
	virtual void Apply(const ParMultiVector& X, const ParMultiVector& Y) const;
};

class ParOperatorCompound : public ParOperator
{
private:
	const ParOperator* A;
	const ParOperator* B;
	const ParOperator* C;
	double alpha;
	double beta;
	double gamma;
	ParVector p;
	ParVector q;

public:	
	ParOperatorCompound(double beta, const ParOperator& C, double gamma = 0.0);	
	ParOperatorCompound(double alpha, const ParOperator& B, double beta, const ParOperator& C, double gamma = 0.0);
	ParOperatorCompound(double alpha, const ParOperator& A, const ParOperator& B, double gamma = 0.0);
	ParOperatorCompound(double alpha, const ParOperator& A, const ParOperator& B, double beta, const ParOperator& C, double gamma = 0.0);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& x, const ParVector& y) const;
	void Apply(const ParMultiVector& X, const ParMultiVector& Y) const;
};

class ParOperatorPoly : public ParOperator
{
private:
	int order;
	const double* coeff;
	const ParOperator* A;
	ParVector p;

public:
	ParOperatorPoly(int order, const double* coeff, const ParOperator& A);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& x, const ParVector& y) const;
	void Apply(const ParMultiVector& X, const ParMultiVector& Y) const;
};

class ParOperatorChebyshev : public ParOperator
{
private:
	int order;
	double lambda;
	double eigen_ratio;
	const ParOperator* A;
	ParVector z;
	ParVector r;

public:
	ParOperatorChebyshev(int order, double eigen_ratio, const ParOperator& A);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& x, const ParVector& y) const;
	void Apply(const ParMultiVector& X, const ParMultiVector& Y) const;
};

class ParOperatorIdent : public ParOperator
{
private:
	int size[2];

public:
	ParOperatorIdent(MPI_Comm comm, int n);
	ParOperatorIdent(MPI_Comm comm, int n, int m);	
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& x, const ParVector& y) const;
};

class ParOperatorZero : public ParOperator
{
private:
	int size[2];

public:
	ParOperatorZero(MPI_Comm comm, int n);
	ParOperatorZero(MPI_Comm comm, int n, int m);	
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& x, const ParVector& y) const;
};

class ParSolver
{
public:
	virtual ~ParSolver();
	virtual void operator()(const ParOperator& A, const ParVector& b, const ParVector& x, int& iter, double& relres) const;
	virtual void operator()(const ParOperator& A, const ParOperator& P, const ParVector& b, const ParVector& x, int& iter, double& relres) const = 0;
};

class ParSolverCG : public ParSolver
{
public:
	int    MaxIters;
	double Tolerance;
	int    PrintStats;

	using ParSolver::operator();
	ParSolverCG(int max_iters = 500, double tolerance = 1.0e-08, int print_stats = 0);
	void operator()(const ParOperator& A, const ParOperator& P, const ParVector& b, const ParVector& x, int& iter, double& relres) const;
};

class ParSolverGMRES : public ParSolver
{
public:
	int    Restart;
	int    MaxIters;
	double Tolerance;
	int    PrintStats;

	using ParSolver::operator();
	ParSolverGMRES(int restart = 10, int max_iters = 500, double tolerance = 1.0e-08, int print_stats = 0);
	void operator()(const ParOperator& A, const ParOperator& P, const ParVector& b, const ParVector& x, int& iter, double& relres) const;
};

class ParSolverBCGS : public ParSolver
{
public:
	int    MaxIters;
	double Tolerance;
	int    PrintStats;

	using ParSolver::operator();
	ParSolverBCGS(int max_iters = 500, double tolerance = 1.0e-08, int print_stats = 0);
	void operator()(const ParOperator& A, const ParOperator& P, const ParVector& b, const ParVector& x, int& iter, double& relres) const;
};

class ParSolverPipeCG : public ParSolver
{
public:
	int    MaxIters;
	double Tolerance;
	int    PrintStats;

	using ParSolver::operator();
	ParSolverPipeCG(int max_iters = 500, double tolerance = 1.0e-08, int print_stats = 0);
	void operator()(const ParOperator& A, const ParOperator& P, const ParVector& b, const ParVector& x, int& iter, double& relres) const;
};

class ParSolverPipeGMRES : public ParSolver
{
public:
	int    Restart;
	double Shift;
	int    MaxIters;
	double Tolerance;
	int    PrintStats;

	using ParSolver::operator();
	ParSolverPipeGMRES(int restart = 10, double shift = 1.0, int max_iters = 500, double tolerance = 1.0e-08, int print_stats = 0);
	void operator()(const ParOperator& A, const ParOperator& P, const ParVector& b, const ParVector& x, int& iter, double& relres) const;
};

class ParSolverPipeBCGS : public ParSolver
{
public:
	int    MaxIters;
	double Tolerance;
	int    PrintStats;

	using ParSolver::operator();
	ParSolverPipeBCGS(int max_iters = 500, double tolerance = 1.0e-08, int print_stats = 0);
	void operator()(const ParOperator& A, const ParOperator& P, const ParVector& b, const ParVector& x, int& iter, double& relres) const;
};

class ParEigenSolver
{
public:
	virtual ~ParEigenSolver();
	virtual void operator()(const ParOperator& A, const ParVector& x, double& lambda, int& iter, double& res) const;
	virtual void operator()(const ParOperator& A, const ParOperator& T, const ParVector& x, double& lambda, int& iter, double& res) const;
	virtual void operator()(const ParOperator& A, const ParOperator& B, const ParOperator& T, const ParVector& x, double& lambda, int& iter, double& res) const;
	virtual void operator()(const ParOperator& A, const ParMultiVector& Y, const ParVector& x, double& lambda, int& iter, double& res) const;
	virtual void operator()(const ParOperator& A, const ParOperator& T, const ParMultiVector& Y, const ParVector& x, double& lambda, int& iter, double& res) const;
	virtual void operator()(const ParOperator& A, const ParOperator& B, const ParOperator& T, const ParMultiVector& Y, const ParVector& x, double& lambda, int& iter, double& res) const;
	virtual void operator()(const ParOperator& A, const ParMultiVector& X, const Vector& Lambda, int& iter, const Vector& Res) const;
	virtual void operator()(const ParOperator& A, const ParOperator& T, const ParMultiVector& X, const Vector& Lambda, int& iter, const Vector& Res) const;
	virtual void operator()(const ParOperator& A, const ParOperator& B, const ParOperator& T, const ParMultiVector& X, const Vector& Lambda, int& iter, const Vector& Res) const;
	virtual void operator()(const ParOperator& A, const ParMultiVector& Y, const ParMultiVector& X, const Vector& Lambda, int& iter, const Vector& Res) const;
	virtual void operator()(const ParOperator& A, const ParOperator& T, const ParMultiVector& Y, const ParMultiVector& X, const Vector& Lambda, int& iter, const Vector& Res) const;
	virtual void operator()(const ParOperator& A, const ParOperator& B, const ParOperator& T, const ParMultiVector& Y, const ParMultiVector& X, const Vector& Lambda, int& iter, const Vector& Res) const = 0;
};

class ParEigenSolverLOBPCG : public ParEigenSolver
{
public:
	int    MaxIters;
	double Tolerance;
	int    PrintStats;

	using ParEigenSolver::operator();
	ParEigenSolverLOBPCG(int max_iters = 500, double tolerance = 1.0e-08, int print_stats = 0);
	void operator()(const ParOperator& A, const ParOperator& B, const ParOperator& T, const ParMultiVector& Y, const ParMultiVector& X, const Vector& Lambda, int& iter, const Vector& Res) const;
};

struct COOMatrix : public Operator
{
	int ref;
	int size[2];
	int nnz;
	
	int* rowind;
	int* colind;
	double* values;

	COOMatrix();
	COOMatrix(int n, int m, int nnz, int* rowind, int* colind, double* values, int ref);
	COOMatrix(const COOMatrix& A);
	COOMatrix(COOMatrix&& A);
	~COOMatrix();
	COOMatrix& operator=(const COOMatrix& A);
	COOMatrix& operator=(COOMatrix&& A);
	
	void Free();
	void Refer(const COOMatrix& A);
	int InSize() const;
	int OutSize() const;
	void Apply(const Vector& x, const Vector& y) const;
};

void COORead(const char* filename, COOMatrix& A);
void COOWrite(const char* filename, const COOMatrix& A);

struct CSRMatrix : public Operator
{
	int ref;
	int size[2];

	int* rowptr;
	int* colind;
	double* values;

	CSRMatrix();
	CSRMatrix(int n, int m, int* rowptr, int* colind, double* values, int ref);
	CSRMatrix(const COOMatrix& A);
	CSRMatrix(const CSRMatrix& A);
	CSRMatrix(CSRMatrix&& A);
	~CSRMatrix();
	CSRMatrix& operator=(const COOMatrix& A);
	CSRMatrix& operator=(const CSRMatrix& A);
	CSRMatrix& operator=(CSRMatrix&& A);

	void Free();
	void Refer(const CSRMatrix& A);
	int InSize() const;
	int OutSize() const;
	void Apply(const Vector& x, const Vector& y) const;
	void Apply(const MultiVector& X, const MultiVector& Y) const;
};

void CSRRead(const char* filename, CSRMatrix& A);
void CSRWrite(const char* filename, const CSRMatrix& A);
void CSRDiag(const CSRMatrix& A, const Vector& D);
void CSREliminZeros(const CSRMatrix& A);
void CSRScale(double alpha, const CSRMatrix& A);
void CSRScaleRows(const Vector& x, const CSRMatrix& A);
void CSRScaleCols(const Vector& x, const CSRMatrix& A);
void CSRMatAdd(const CSRMatrix& A, const CSRMatrix& B, CSRMatrix& C);
void CSRMatMul(const CSRMatrix& A, const CSRMatrix& B, CSRMatrix& C);
void CSRMatVec(double alpha, const CSRMatrix& A, const Vector& x, double beta, const Vector& y);
void CSRMatMultiVec(double alpha, const CSRMatrix& A, const MultiVector& X, double beta, const MultiVector& Y);
void CSRTrans(const CSRMatrix& A, CSRMatrix& B);

class CSRCreator
{
private:
	int size[2];

public:
	CSRCreator(int n, int m);
	void operator()(CSRMatrix& A) const;
};

class CSRMutator
{
private:
	CSRMatrix A;
	bool* marks;
	void* data;

public:
	CSRMutator(const CSRMatrix& A);
	~CSRMutator();
	double* Find(int row, int col) const;
	double* Insert(int row, int col) const;
	void Erase(int row, int col) const;
	void operator()(CSRMatrix& A);
};

class CSRAccessor
{
private:
	CSRMatrix A;

public:
	CSRAccessor(const CSRMatrix& A);
	double* Find(int row, int col) const;
};

struct ParCSRMatrix : public ParOperator
{
	CSRMatrix local;
	CSRMatrix exter;

	int nnb;
	int* nbrank;

	int* recvptr;
	int* recvind;

	int* sendptr;
	int* sendind;
	double* sendbuf;
	
	Vector recvx;

	ParCSRMatrix(MPI_Comm comm = MPI_COMM_SELF);
	ParCSRMatrix(MPI_Comm comm, int local_rows, int local_cols, int exter_cols, 
	int* local_rowptr, int* local_colind, double* local_values, int local_ref, 
	int* exter_rowptr, int* exter_colind, double* exter_values, int exter_ref,
	int nnb, int* nbrank, int* recvptr, int* recvind);
	ParCSRMatrix(const CSRMatrix& A);
	ParCSRMatrix(const ParCSRMatrix& A);
	ParCSRMatrix(ParCSRMatrix&& A);
	~ParCSRMatrix();
	ParCSRMatrix& operator=(const ParCSRMatrix& A);
	ParCSRMatrix& operator=(ParCSRMatrix&& A);
	
	void Free();
	void Refer(const ParCSRMatrix& A);
	void ExchangeHalo(const ParVector& x) const;
	void SetupHalo();
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& x, const ParVector& y) const;
	void Apply(const ParMultiVector& X, const ParMultiVector& Y) const;
};

void ParCSRDiag(const ParCSRMatrix& A, const ParVector& D);
void ParCSREliminZeros(const ParCSRMatrix& A);
void ParCSRScale(double alpha, const ParCSRMatrix& A);
void ParCSRScaleRows(const ParVector& x, const ParCSRMatrix& A);
void ParCSRScaleCols(const ParVector& x, const ParCSRMatrix& A);
void ParCSRMatAdd(const ParCSRMatrix& A, const ParCSRMatrix& B, ParCSRMatrix& C);
void ParCSRMatMul(const ParCSRMatrix& A, const ParCSRMatrix& B, ParCSRMatrix& C);
void ParCSRMatVec(double alpha, const ParCSRMatrix& A, const ParVector& x, double beta, const ParVector& y);
void ParCSRMatMultiVec(double alpha, const ParCSRMatrix& A, const ParMultiVector& X, double beta, const ParMultiVector& Y);
void ParCSRTrans(const ParCSRMatrix& A, ParCSRMatrix& B);

class ParCSRCreator
{
private:
	MPI_Comm comm;
	int size[2];

public:
	ParCSRCreator(MPI_Comm comm, int n, int m);
	void operator()(ParCSRMatrix& A) const;
};

class ParCSRMutator
{
private:
	MPI_Comm comm;
	CSRMatrix locref;
	CSRMatrix extref;
	int* recvind;
	int* recvfrom;
	bool* marks;
	void* data;

public:
	ParCSRMutator(const ParCSRMatrix& A);
	~ParCSRMutator();
	double* Find(int row, global col) const;
	double* Insert(int row, global col) const;
	void Erase(int row, global col) const;
	void operator()(ParCSRMatrix& A);
};

class ParCSRAccessor
{
private:
	MPI_Comm comm;
	CSRMatrix locref;
	CSRMatrix extref;
	int* recvind;
	int* recvfrom;

public:
	ParCSRAccessor(const ParCSRMatrix& A);
	~ParCSRAccessor();
	double* Find(int row, global col) const;
};

class ParCSRGenerator
{
public:
	MPI_Comm comm;

	ParCSRGenerator(MPI_Comm comm = MPI_COMM_SELF);
	virtual ~ParCSRGenerator();
	virtual void operator()(ParCSRMatrix& A) const = 0;
};

class ParCSRGeneratorIJ : public ParCSRGenerator
{
public:
	int size[2];
	const int* rowptr;
	const global* colind;
	const double* values;
	int pattern_symmetry;

	ParCSRGeneratorIJ(MPI_Comm comm = MPI_COMM_SELF);
	ParCSRGeneratorIJ(MPI_Comm comm, int n, int m, const int* rowptr, const global* colind, const double* values, int pattern_symmetry = 0);
	void operator()(ParCSRMatrix& A) const;
};

class ParCSRGenerator7pt : public ParCSRGenerator
{
public:
	int nx;
	int ny;
	int nz;

	int Px;
	int Py;
	int Pz;

	double cx;
	double cy;
	double cz;

	ParCSRGenerator7pt(MPI_Comm comm = MPI_COMM_SELF);
	ParCSRGenerator7pt(MPI_Comm comm, int nx, int ny, int nz, int Px, int Py, int Pz, double cx, double cy, double cz);
	void operator()(ParCSRMatrix& A) const;
};

class ParCSRGenerator27pt : public ParCSRGenerator
{
public:
	int nx;
	int ny;
	int nz;

	int Px;
	int Py;
	int Pz;

	ParCSRGenerator27pt(MPI_Comm comm = MPI_COMM_SELF);
	ParCSRGenerator27pt(MPI_Comm comm, int nx, int ny, int nz, int Px, int Py, int Pz);
	void operator()(ParCSRMatrix& A) const;
};

class ParCSRPrecond : public ParOperator
{
public:
	virtual ~ParCSRPrecond();
	virtual void Setup(const ParCSRMatrix& A, int REUSE = 0) = 0;
};

class ParCSRPrecondJacobi : public ParCSRPrecond
{
private:
	Vector D_recip;

public:
	void Free();
	void Setup(const ParCSRMatrix& A, int REUSE = 0);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& b, const ParVector& x) const;
};

class ParCSRPrecondSOR : public ParCSRPrecond
{
private:
	Vector D_local_recip;
	CSRMatrix L_local;
	CSRMatrix L_exter;
	CSRMatrix U_local;
	CSRMatrix U_exter;
	Vector recvx;

	bool* interior;

	int label;
	int nnb;
	int* nblab;
	int* nbrank;

	int* recvptr;
	int* sendptr;
	int* sendind;
	double* sendbuf;

public:
	int RelaxationType;
	double RelaxationFactor;

	ParCSRPrecondSOR(int relaxation_type = 2, double relaxation_factor = 1.0);
	~ParCSRPrecondSOR();

	void Free();
	void Setup(const ParCSRMatrix& A, int REUSE = 0);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& b, const ParVector& x) const;
};

class ParCSRPrecondILU0 : public ParCSRPrecond
{
private:
	Vector D_local_recip;
	CSRMatrix L_local;
	CSRMatrix L_exter;
	CSRMatrix U_local;
	CSRMatrix U_exter;
	Vector recvx;

	bool* interior;

	int label;
	int nnb;
	int* nblab;
	int* nbrank;

	int* recvptr;
	int* sendptr;
	int* sendind;
	double* sendbuf;

public:
	ParCSRPrecondILU0();
	~ParCSRPrecondILU0();
	
	void Free();
	void Setup(const ParCSRMatrix& A, int REUSE = 0);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& b, const ParVector& x) const;
};

class ParCSRPrecondBJSOR : public ParCSRPrecond
{
private:
	int nthd;
	Vector D_recip;
	CSRMatrix L;
	CSRMatrix U;

public:
	int RelaxationType;
	double RelaxationFactor;

	ParCSRPrecondBJSOR(int relaxation_type = 2, double relaxation_factor = 1.0);

	void Free();
	void Setup(const ParCSRMatrix& A, int REUSE = 0);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& b, const ParVector& x) const;
};

class ParCSRPrecondBJILU : public ParCSRPrecond
{
private:
	int nthd;
	Vector* D_recip;
	CSRMatrix* L;
	CSRMatrix* U;

public:
	int MaxFillins;
	double DropTolerance;

	ParCSRPrecondBJILU(int maxfil = 0, double droptol = 0.01);
	~ParCSRPrecondBJILU();

	void Free();
	int InSize() const;
	int OutSize() const;
	void Setup(const ParCSRMatrix& A, int REUSE = 0);
	void Apply(const ParVector& b, const ParVector& x) const;
};

class ParCSRPrecondAMG : public ParCSRPrecond
{
private:
	class Smoother;

	int nlev;
	ParCSRMatrix* A;
	ParCSRMatrix* P;
	ParCSRMatrix* R;
	Smoother* S;

	ParVector* r;
	ParVector* e;
	ParVector* g;

	void Strength(const ParCSRMatrix& A, bool* locstg, bool* extstg) const;
	void Coarsening(const ParCSRMatrix& A, const bool* locstg, const bool* extstg, int* cfmap, int AGGRESSIVE = 0) const;
	void Interpolation(const ParCSRMatrix& A, const bool* locstg, const bool* extstg, int* cfmap, ParCSRMatrix& P, int AGGRESSIVE = 0) const;
	void RAP(const ParCSRMatrix& R, const ParCSRMatrix& A, const ParCSRMatrix& P, ParCSRMatrix& C) const;
	void Sparsification(const ParCSRMatrix& A, ParCSRMatrix& B) const;
	void SetupSmoother(const ParCSRMatrix& A, Smoother& S, int REUSE = 0) const;
	void PreSmooth(int level, const ParVector& b, const ParVector& x, bool x0zero = 0) const;
	void PreSmooth_Restriction(int level, const ParVector& b, const ParVector& x, const ParVector& g, bool x0zero = 0) const;
	void PostSmooth(int level, const ParVector& b, const ParVector& x) const;
	void CoarseSmooth(int level, const ParVector& b, const ParVector& x, bool x0zero = 0) const;
	void Restriction(int level, const ParVector& b, const ParVector& x, const ParVector& g) const;
	void Prolongation(int level, const ParVector& e, const ParVector& x) const;
	void V_Cycle(int level, const ParVector& b, const ParVector& x, bool x0zero = 0) const;
	void W_Cycle(int level, const ParVector& b, const ParVector& x, bool x0zero = 0) const;
	void F_Cycle(int level, const ParVector& b, const ParVector& x) const;

public:
	int    MaxLevels;
	int    CoarseSize;
	double StrengthThreshold;
	int    AggressiveLevels;
	int    CoarsenType;
	int    InterpType;
	int    InterpMinElements;
	int    InterpMaxElements;
	double TruncationFactor;
	double SparsificationThreshold;
	int    CycleType;
	int    PreSweeps;
	int    PostSweeps;
	int    CoarseSweeps;
	int    SmoothType;
	double JacobiFactor;
	int    SORType;
	double SORFactor;
	int    ILUMaxFillins;
	double ILUDropTolerance;
	double ChebyshevEigenRatio;
	int    PrintStats;

	ParCSRPrecondAMG(
	int    max_levels = 20,
	int    coarse_size = 8,
	double strength_threshold = 0.25,
	int    aggressive_levels = 0,
	int    coarsen_type = 0,
	int    interp_type = 0,
	int    interp_min_elements = 4,
	int    interp_max_elements = 4,	
	double truncation_factor = 0.1,
	double sparsification_threshold = 0.0,
	int    cycle_type = 0,
	int    pre_sweeps = 1,
	int    post_sweeps = 1,
	int    coarse_sweeps = 2,
	int    smooth_type = 1,
	int    print_stats = 0);
	~ParCSRPrecondAMG();
	
	void Free();
	void Setup(const ParCSRMatrix& A, int REUSE = 0);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& b, const ParVector& x) const;
};

struct BSRMatrix : public Operator
{
	bool ref;
	int size[2];
	int bsize;

	int* rowptr;
	int* colind;
	double* values;

	BSRMatrix();
	BSRMatrix(int n, int m, int bsize, int* rowptr, int* colind, double* values, int ref);
	BSRMatrix(const BSRMatrix& A);
	BSRMatrix(int bsize, const CSRMatrix& A);
	BSRMatrix(BSRMatrix&& A);
	~BSRMatrix();
	BSRMatrix& operator=(const BSRMatrix& A);
	BSRMatrix& operator=(BSRMatrix&& A);

	void Free();
	void Refer(const BSRMatrix& A);
	int InSize() const;
	int OutSize() const;
	void Apply(const Vector& x, const Vector& y) const;
};

void CSRToBSR(int bsize, const CSRMatrix& A, BSRMatrix& B);
void BSRDiag(const BSRMatrix& A, const Vector& D);
void BSREliminZeros(const BSRMatrix& A);
void BSRScale(double alpha, const BSRMatrix& A);
void BSRScaleRows(const Vector& x, const BSRMatrix& A);
void BSRScaleCols(const Vector& x, const BSRMatrix& A);
void BSRMatAdd(const BSRMatrix& A, const BSRMatrix& B, BSRMatrix& C);
void BSRMatMul(const BSRMatrix& A, const BSRMatrix& B, BSRMatrix& C);
void BSRMatVec(double alpha, const BSRMatrix& A, const Vector& x, double beta, const Vector& y);
void BSRTrans(const BSRMatrix& A, BSRMatrix& B);

static inline void VecBlockFill(int bsize, double a, double* x)
{
	for (int i = 0; i < bsize; ++i)
		x[i] = a;
}

static inline void VecBlockScale(int bsize, double a, double* x)
{
	for (int i = 0; i < bsize; ++i)
		x[i] *= a;
}

static inline void VecBlockCopy(int bsize, const double* x, double* y)
{
	for (int i = 0; i < bsize; ++i)
		y[i] = x[i];
}

static inline void VecBlockAdd(int bsize, const double* x, double* y)
{
	for (int i = 0; i < bsize; ++i)
		y[i] += x[i];
}

static inline void VecBlockSub(int bsize, const double* x, double* y)
{
	for (int i = 0; i < bsize; ++i)
		y[i] -= x[i];
}

static inline void VecBlockScaleAdd(int bsize, double a, const double* x, double* y)
{
	for (int i = 0; i < bsize; ++i)
		y[i] += a * x[i];
}

static inline void BSRBlockFill(int bsize, double a, double* A)
{
	int bnnz = bsize * bsize;
	for (int i = 0; i < bnnz; ++i)
		A[i] = a;
}

static inline void BSRBlockScale(int bsize, double a, double* A)
{
	int bnnz = bsize * bsize;
	for (int i = 0; i < bnnz; ++i)
		A[i] *= a;
}

static inline void BSRBlockCopy(int bsize, const double* A, double* B)
{
	int bnnz = bsize * bsize;
	for (int i = 0; i < bnnz; ++i)
		B[i] = A[i];
}

static inline void BSRBlockTranspose(int bsize, const double* A, double* B)
{
	for (int i = 0; i < bsize; ++i)
		for (int j = 0; j < bsize; ++j)
			B[i + j * bsize] = A[i * bsize + j];
}

static inline void BSRBlockAdd(int bsize, const double* A, double* B)
{
	int bnnz = bsize * bsize;
	for (int i = 0; i < bnnz; ++i)
		B[i] += A[i];
}

static inline void BSRBlockSub(int bsize, const double* A, double* B)
{
	int bnnz = bsize * bsize;
	for (int i = 0; i < bnnz; ++i)
		B[i] -= A[i];
}

static inline void BSRBlockScaleAdd(int bsize, double a, const double* A, double* B)
{
	int bnnz = bsize * bsize;
	for (int i = 0; i < bnnz; ++i)
		B[i] += a * A[i];
}

static inline void BSRBlockMatVec(int bsize, const double* A, const double* x, double* y)
{
	for (int i = 0; i < bsize; ++i)
	{
		double temp = 0.0;
		for (int j = 0; j < bsize; ++j)
			temp += A[i + j * bsize] * x[j];
		y[i] = temp;
	}
}

static inline void BSRBlockMatVecAdd(int bsize, const double* A, const double* x, double* y)
{
	for (int i = 0; i < bsize; ++i)
	{
		double temp = y[i];
		for (int j = 0; j < bsize; ++j)
			temp += A[i + j * bsize] * x[j];
		y[i] = temp;
	}
}

static inline void BSRBlockMatVecSub(int bsize, const double* A, const double* x, double* y)
{
	for (int i = 0; i < bsize; ++i)
	{
		double temp = y[i];
		for (int j = 0; j < bsize; ++j)
			temp -= A[i + j * bsize] * x[j];
		y[i] = temp;
	}
}

static inline void BSRBlockMatMul(int bsize, const double* A, const double* B, double* C)
{
	for (int j = 0; j < bsize; ++j)
	{
		for (int i = 0; i < bsize; ++i)
		{
			double temp = 0.0;
			for (int k = 0; k < bsize; ++k)
				temp += A[i + k * bsize] * B[k + j * bsize];
			C[i + j * bsize] = temp;
		}
	}
}

static inline void BSRBlockMatMulAdd(int bsize, const double* A, const double* B, double* C)
{
	for (int j = 0; j < bsize; ++j)
	{
		for (int i = 0; i < bsize; ++i)
		{
			double temp = C[i + j * bsize];
			for (int k = 0; k < bsize; ++k)
				temp += A[i + k * bsize] * B[k + j * bsize];
			C[i + j * bsize] = temp;
		}
	}
}

static inline void BSRBlockMatMulSub(int bsize, const double* A, const double* B, double* C)
{
	for (int j = 0; j < bsize; ++j)
	{
		for (int i = 0; i < bsize; ++i)
		{
			double temp = C[i + j * bsize];
			for (int k = 0; k < bsize; ++k)
				temp -= A[i + k * bsize] * B[k + j * bsize];
			C[i + j * bsize] = temp;
		}
	}
}

static inline void BSRBlockLUFactorize(int bsize, double* LU)
{
	for (int k = 0; k < bsize; ++k)
	{
		double pivot = 1.0 / LU[k + k * bsize];
		LU[k + k * bsize] = pivot;
		for (int i = k + 1; i < bsize; ++i)
		{
			double factor = LU[i + k * bsize] * pivot;
			LU[i + k * bsize] = factor;
			for (int j = k + 1; j < bsize; ++j)
				LU[i + j * bsize] -= factor * LU[k + j * bsize];
		}
	}
}

static inline void BSRBlockLUVecSolve(int bsize, const double* LU, double* x)
{
	for (int i = 1; i < bsize; ++i)
	{
		double temp = x[i];
		for (int k = 0; k < i; ++k)
			temp -= x[k] * LU[i + k * bsize];
		x[i] = temp;
	}

	for (int i = bsize - 1; i >= 0; --i)
	{
		double temp = x[i];
		for (int k = bsize - 1; k > i; --k)
			temp -= x[k] * LU[i + k * bsize];
		x[i] = temp * LU[i + i * bsize];
	}
}

static inline void BSRBlockLUMatSolve(int bsize, const double* LU, double* X)
{
	for (int j = 0; j < bsize; ++j)
	{
		for (int i = 1; i < bsize; ++i)
		{
			double temp = X[i + j * bsize];
			for (int k = 0; k < i; ++k)
				temp -= X[k + j * bsize] * LU[i + k * bsize];
			X[i + j * bsize] = temp;
		}
	}
	
	for (int j = bsize - 1; j >= 0; --j)
	{
		for (int i = bsize - 1; i >= 0; --i)
		{
			double temp = X[i + j * bsize];
			for (int k = bsize - 1; k > i; --k)
				temp -= X[k + j * bsize] * LU[i + k * bsize];
			X[i + j * bsize] = temp * LU[i + i * bsize];
		}
	}
}

static inline void BSRBlockMatLUSolve(int bsize, const double* LU, double* X)
{
	for (int j = 0; j < bsize; ++j)
	{
		X[j] *= LU[0];
		for (int i = 1; i < bsize; ++i)
		{
			double temp = X[j + i * bsize];
			for (int k = 0; k < i; ++k)
				temp -= X[j + k * bsize] * LU[k + i * bsize];
			X[j + i * bsize] = temp * LU[i + i * bsize];
		}
	}

	for (int j = bsize - 1; j >= 0; --j)
	{
		for (int i = bsize - 2; i >= 0; --i)
		{
			double temp = X[j + i * bsize];
			for (int k = bsize - 1; k > i; --k)
				temp -= X[j + k * bsize] * LU[k + i * bsize];
			X[j + i * bsize] = temp;
		}
	}
}

static inline double BSRBlockSum(int bsize, const double* A)
{
	int bnnz = bsize * bsize;
	double result = 0.0;
	for (int i = 0; i < bnnz; ++i)
		result += A[i];
	return result;
}

static inline void BSRBlockDiagonal(int bsize, const double* A, double* x)
{
	for (int i = 0; i < bsize; ++i)
		x[i] = A[i + i * bsize];
}

static inline bool BSRBlockAllZero(int bsize, const double* A)
{
	int bnnz = bsize * bsize;
	for (int i = 0; i < bnnz; ++i)
		if (A[i] != 0.0) return 0;
	return 1;
}

static inline void BSRBlockScaleRows(int bsize, const double* x, double* A)
{
	for (int j = 0; j < bsize; ++j)
		for (int i = 0; i < bsize; ++i)
			A[i + j * bsize] *= x[i];
}

static inline void BSRBlockScaleCols(int bsize, const double* x, double* A)
{
	for (int j = 0; j < bsize; ++j)
		for (int i = 0; i < bsize; ++i)
			A[i + j * bsize] *= x[j];
}

class BSRCreator
{
private:
	int size[2];
	int bsize;

public:
	BSRCreator(int n, int m, int b);
	void operator()(BSRMatrix& A) const;
};

class BSRMutator
{
private:
	BSRMatrix A;
	bool* marks;
	void* data;

public:
	BSRMutator(const BSRMatrix& A);
	~BSRMutator();
	double* Find(int row, int col) const;
	double* Insert(int row, int col) const;
	void Erase(int row, int col) const;
	void operator()(BSRMatrix& A);
};

class BSRAccessor
{
private:
	BSRMatrix A;

public:
	BSRAccessor(const BSRMatrix& A);
	double* Find(int row, int col) const;
};

struct ParBSRMatrix : public ParOperator
{
	BSRMatrix local;
	BSRMatrix exter;

	int nnb;
	int* nbrank;

	int* recvptr;
	int* recvind;

	int* sendptr;
	int* sendind;
	double* sendbuf;
	
	Vector recvx;

	ParBSRMatrix(MPI_Comm comm = MPI_COMM_SELF);
	ParBSRMatrix(MPI_Comm comm, int local_rows, int local_cols, int exter_cols, int bsize, 
	int* local_rowptr, int* local_colind, double* local_values, int local_ref, 
	int* exter_rowptr, int* exter_colind, double* exter_values, int exter_ref,
	int nnb, int* nbrank, int* recvptr, int* recvind);
	ParBSRMatrix(const BSRMatrix& A);
	ParBSRMatrix(const ParBSRMatrix& A);
	ParBSRMatrix(int bsize, const ParCSRMatrix& A);
	ParBSRMatrix(ParBSRMatrix&& A);
	~ParBSRMatrix();
	ParBSRMatrix& operator=(const ParBSRMatrix& A);
	ParBSRMatrix& operator=(ParBSRMatrix&& A);
	
	void Free();
	void Refer(const ParBSRMatrix& A);
	void ExchangeHalo(const ParVector& x) const;
	void SetupHalo();
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& x, const ParVector& y) const;
};

void ParCSRToBSR(int bsize, const ParCSRMatrix& A, ParBSRMatrix& B);
void ParBSRDiag(const ParBSRMatrix& A, const ParVector& D);
void ParBSREliminZeros(const ParBSRMatrix& A);
void ParBSRScale(double alpha, const ParBSRMatrix& A);
void ParBSRScaleRows(const ParVector& x, const ParBSRMatrix& A);
void ParBSRScaleCols(const ParVector& x, const ParBSRMatrix& A);
void ParBSRMatAdd(const ParBSRMatrix& A, const ParBSRMatrix& B, ParBSRMatrix& C);
void ParBSRMatMul(const ParBSRMatrix& A, const ParBSRMatrix& B, ParBSRMatrix& C);
void ParBSRMatVec(double alpha, const ParBSRMatrix& A, const ParVector& x, double beta, const ParVector& y);
void ParBSRTrans(const ParBSRMatrix& A, ParBSRMatrix& B);

class ParBSRCreator
{
private:
	MPI_Comm comm;
	int size[2];
	int bsize;

public:
	ParBSRCreator(MPI_Comm comm, int n, int m, int b);
	void operator()(ParBSRMatrix& A) const;
};

class ParBSRMutator
{
private:
	MPI_Comm comm;
	BSRMatrix locref;
	BSRMatrix extref;
	int* recvind;
	int* recvfrom;
	bool* marks;
	void* data;

public:
	ParBSRMutator(const ParBSRMatrix& A);
	~ParBSRMutator();
	double* Find(int row, global col) const;
	double* Insert(int row, global col) const;
	void Erase(int row, global col) const;
	void operator()(ParBSRMatrix& A);
};

class ParBSRAccessor
{
private:
	MPI_Comm comm;
	BSRMatrix locref;
	BSRMatrix extref;
	int* recvind;
	int* recvfrom;

public:
	ParBSRAccessor(const ParBSRMatrix& A);
	~ParBSRAccessor();
	double* Find(int row, global col) const;
};

class ParBSRGenerator
{
public:
	MPI_Comm comm;

	ParBSRGenerator(MPI_Comm comm = MPI_COMM_SELF);
	virtual ~ParBSRGenerator();
	virtual void operator()(ParBSRMatrix& A) const = 0;
};

class ParBSRGeneratorIJ : public ParBSRGenerator
{
public:
	int size[2];
	int bsize;
	const int* rowptr;
	const global* colind;
	const double* values;
	int pattern_symmetry;

	ParBSRGeneratorIJ(MPI_Comm comm = MPI_COMM_SELF);
	ParBSRGeneratorIJ(MPI_Comm comm, int n, int m, int bsize, const int* rowptr, const global* colind, const double* values, int pattern_symmetry = 0);
	void operator()(ParBSRMatrix& A) const;
};

class ParBSRGenerator7pt : public ParBSRGenerator
{
public:
	int nx;
	int ny;
	int nz;

	int bx;
	int by;
	int bz;

	int Px;
	int Py;
	int Pz;

	double cx;
	double cy;
	double cz;

	ParBSRGenerator7pt(MPI_Comm comm = MPI_COMM_SELF);
	ParBSRGenerator7pt(MPI_Comm comm, int nx, int ny, int nz, int bx, int by, int bz, int Px, int Py, int Pz, double cx, double cy, double cz);
	void operator()(ParBSRMatrix& A) const;
};

class ParBSRGenerator27pt : public ParBSRGenerator
{
public:
	int nx;
	int ny;
	int nz;

	int bx;
	int by;
	int bz;

	int Px;
	int Py;
	int Pz;

	ParBSRGenerator27pt(MPI_Comm comm = MPI_COMM_SELF);
	ParBSRGenerator27pt(MPI_Comm comm, int nx, int ny, int nz, int bx, int by, int bz, int Px, int Py, int Pz);
	void operator()(ParBSRMatrix& A) const;
};

class ParBSRPrecond : public ParOperator
{
public:
	virtual ~ParBSRPrecond();
	virtual void Setup(const ParBSRMatrix& A, int REUSE = 0) = 0;

};

class ParBSRPrecondJacobi : public ParBSRPrecond
{
private:
	int size;
	int bsize;
	double* D_LU;

public:
	ParBSRPrecondJacobi();
	~ParBSRPrecondJacobi();
	void Free();
	void Setup(const ParBSRMatrix& A, int REUSE = 0);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& b, const ParVector& x) const;
};

class ParBSRPrecondSOR : public ParBSRPrecond
{
private:
	double* D_local_LU;
	BSRMatrix L_local;
	BSRMatrix L_exter;
	BSRMatrix U_local;
	BSRMatrix U_exter;
	Vector recvx;

	bool* interior;

	int label;
	int nnb;
	int* nblab;
	int* nbrank;

	int* recvptr;
	int* sendptr;
	int* sendind;
	double* sendbuf;

public:
	int RelaxationType;
	double RelaxationFactor;

	ParBSRPrecondSOR(int relaxation_type = 2, double relaxation_factor = 1.0);
	~ParBSRPrecondSOR();

	void Free();
	void Setup(const ParBSRMatrix& A, int REUSE = 0);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& b, const ParVector& x) const;
};

class ParBSRPrecondILU0 : public ParBSRPrecond
{
private:
	double* D_local_LU;
	BSRMatrix L_local;
	BSRMatrix L_exter;
	BSRMatrix U_local;
	BSRMatrix U_exter;
	Vector recvx;

	bool* interior;

	int label;
	int nnb;
	int* nblab;
	int* nbrank;

	int* recvptr;
	int* sendptr;
	int* sendind;
	double* sendbuf;

public:
	ParBSRPrecondILU0();
	~ParBSRPrecondILU0();
	
	void Free();
	void Setup(const ParBSRMatrix& A, int REUSE = 0);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& b, const ParVector& x) const;
};

class ParBSRPrecondBJSOR : public ParBSRPrecond
{
private:
	int nthd;
	double* D_LU;
	BSRMatrix L;
	BSRMatrix U;

public:
	int RelaxationType;
	double RelaxationFactor;

	ParBSRPrecondBJSOR(int relaxation_type = 2, double relaxation_factor = 1.0);
	~ParBSRPrecondBJSOR();

	void Free();
	void Setup(const ParBSRMatrix& A, int REUSE = 0);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& b, const ParVector& x) const;
};

class ParBSRPrecondBJILU : public ParBSRPrecond
{
private:
	int nthd;
	double* D_LU;
	BSRMatrix* L;
	BSRMatrix* U;

public:
	int MaxFillins;
	double DropTolerance;

	ParBSRPrecondBJILU(int maxfil = 0, double droptol = 0.01);
	~ParBSRPrecondBJILU();

	void Free();
	int InSize() const;
	int OutSize() const;
	void Setup(const ParBSRMatrix& A, int REUSE = 0);
	void Apply(const ParVector& b, const ParVector& x) const;
};

class ParBSRPrecondAMG : public ParBSRPrecond
{
private:
	int nlev;
	ParBSRMatrix* A;
	ParCSRMatrix* P;
	ParCSRMatrix* R;
	class Smoother;
	Smoother* S;

	ParVector* r;
	ParVector* e;
	ParVector* g;

	void Strength(const ParCSRMatrix& A, bool* locstg, bool* extstg) const;
	void Coarsening(const ParCSRMatrix& A, const bool* locstg, const bool* extstg, int* cfmap, int AGGRESSIVE = 0) const;
	void Interpolation(const ParCSRMatrix& A, const bool* locstg, const bool* extstg, int* cfmap, ParCSRMatrix& P, int AGGRESSIVE = 0) const;
	void RAP(const ParCSRMatrix& R, const ParBSRMatrix& A, const ParCSRMatrix& P, ParBSRMatrix& C) const;
	void Sparsification(const ParBSRMatrix& A, ParBSRMatrix& B) const;
	void SetupSmoother(const ParBSRMatrix& A, Smoother& S, int REUSE = 0) const;
	void PreSmooth(int level, const ParVector& b, const ParVector& x, bool x0zero = 0) const;
	void PostSmooth(int level, const ParVector& b, const ParVector& x) const;
	void CoarseSmooth(int level, const ParVector& b, const ParVector& x, bool x0zero = 0) const;
	void Restriction(int level, const ParVector& b, const ParVector& x, const ParVector& g) const;
	void Prolongation(int level, const ParVector& e, const ParVector& x) const;
	void V_Cycle(int level, const ParVector& b, const ParVector& x, bool x0zero = 0) const;
	void W_Cycle(int level, const ParVector& b, const ParVector& x, bool x0zero = 0) const;
	void F_Cycle(int level, const ParVector& b, const ParVector& x) const;

public:
	int    MaxLevels;
	int    CoarseSize;
	double StrengthThreshold;
	int    AggressiveLevels;
	int    CoarsenType;
	int    InterpType;
	int    InterpMinElements;
	int    InterpMaxElements;
	double TruncationFactor;
	double SparsificationThreshold;
	int    CycleType;
	int    PreSweeps;
	int    PostSweeps;
	int    CoarseSweeps;
	int    SmoothType;
	double JacobiFactor;
	int    SORType;
	double SORFactor;
	int    ILUMaxFillins;
	double ILUDropTolerance;
	double ChebyshevEigenRatio;
	int    PrintStats;

	ParBSRPrecondAMG(
	int    max_levels = 20,
	int    coarse_size = 8,
	double strength_threshold = 0.25,
	int    aggressive_levels = 0,
	int    coarsen_type = 0,
	int    interp_type = 0,
	int    interp_min_elements = 4,
	int    interp_max_elements = 4,	
	double truncation_factor = 0.1,
	double sparsification_threshold = 0.0,
	int    cycle_type = 0,
	int    pre_sweeps = 1,
	int    post_sweeps = 1,
	int    coarse_sweeps = 2,
	int    smooth_type = 1,
	int    print_stats = 0);
	~ParBSRPrecondAMG();
	
	void Free();
	void Setup(const ParBSRMatrix& A, int REUSE = 0);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParVector& b, const ParVector& x) const;
};

struct DenseMatrix : public Operator
{
	int ref;
	int size[2];
	double* values;

	DenseMatrix();
	DenseMatrix(int n, int m, double* values, int ref);
	DenseMatrix(const DenseMatrix& A);
	DenseMatrix(DenseMatrix&& A);
	DenseMatrix(const COOMatrix& A);
	DenseMatrix(const CSRMatrix& A);
	DenseMatrix(const BSRMatrix& A);
	~DenseMatrix();
	DenseMatrix& operator=(const DenseMatrix& A);
	DenseMatrix& operator=(DenseMatrix&& A);
	DenseMatrix& operator=(const COOMatrix& A);
	DenseMatrix& operator=(const CSRMatrix& A);
	DenseMatrix& operator=(const BSRMatrix& A);
	double& operator()(int i, int j) const;

	void Free();
	void Resize(int n, int m);
	void Refer(const DenseMatrix& A);
	int InSize() const;
	int OutSize() const;
	void Apply(const Vector& x, const Vector& y) const;
	void Apply(const MultiVector& X, const MultiVector& Y) const;
};

class DenseInverse : public Operator
{
private:
	int size;
	int* piv;
	double* LU;

public:
	DenseInverse(const DenseMatrix& A);
	~DenseInverse();
	void Free();
	int InSize() const;
	int OutSize() const;
	void Apply(const Vector& x, const Vector& y) const;
	void Apply(const MultiVector& X, const MultiVector& Y) const;
};

struct ComplexVector
{
	int ref;
	int size;
	zomplex* values;

	ComplexVector();
	ComplexVector(int n);
	ComplexVector(int n, zomplex* values, int ref);
	ComplexVector(const Vector& x);
	ComplexVector(const ComplexVector& x);
	ComplexVector(ComplexVector&& x); 
	~ComplexVector();
	ComplexVector& operator=(zomplex a);
	ComplexVector& operator=(const ComplexVector& x);
	ComplexVector& operator=(ComplexVector&& x);
	zomplex& operator[](int i) const;

	void Free();
	void Resize(int n);
	void Fill(zomplex a) const;
	void FillRandom() const;
	void Copy(const ComplexVector& x) const;
	void Scale(zomplex a) const;
	void AddScaled(zomplex a, const ComplexVector& x) const;
	void Add2Scaled(zomplex a, const ComplexVector& x, zomplex b, const ComplexVector& y) const;
	void Refer(const ComplexVector& x);
};

void ComplexVecRead(const char* filename, ComplexVector& x);
void ComplexVecWrite(const char* filename, const ComplexVector& x);
void ComplexVecAXPBY(zomplex alpha, const ComplexVector& x, zomplex beta, const ComplexVector& y);
void ComplexVecAXPBYPCZ(zomplex alpha, const ComplexVector& x, zomplex beta, const ComplexVector& y, zomplex gamma, const ComplexVector& z);
zomplex ComplexVecConjDot(const ComplexVector& x, const ComplexVector& y);
zomplex ComplexVecDot(const ComplexVector& x, const ComplexVector& y);
void ComplexVecElemMul(const ComplexVector& x, const ComplexVector& y);
void ComplexVecElemMulConj(const ComplexVector& x, const ComplexVector& y);
void ComplexVecConj(const ComplexVector& x);
void ComplexVecRecip(const ComplexVector& x);

struct ParComplexVector
{
	MPI_Comm comm;
	ComplexVector local;

	ParComplexVector(MPI_Comm comm = MPI_COMM_SELF);
	ParComplexVector(MPI_Comm comm, int local_size);
	ParComplexVector(MPI_Comm comm, int local_size, zomplex* local_values, int local_ref);
	ParComplexVector(const ComplexVector& x);
	ParComplexVector(const ParVector& x);
	ParComplexVector(const ParComplexVector& x);
	ParComplexVector(ParComplexVector&& x);
	ParComplexVector& operator=(zomplex a);
	ParComplexVector& operator=(const ParComplexVector& x);
	ParComplexVector& operator=(ParComplexVector&& x);
	zomplex& operator[](int i) const;

	void Free();
	void Resize(int n);
	void Refer(const ParComplexVector& x);
	void Fill(zomplex a) const;
	void FillRandom() const;
	void Copy(const ParComplexVector& x) const;
	void Scale(zomplex a) const;
	void AddScaled(zomplex a, const ParComplexVector& x) const;
	void Add2Scaled(zomplex a, const ParComplexVector& x, zomplex b, const ParComplexVector& y) const;
};

void ParComplexVecAXPBY(zomplex alpha, const ParComplexVector& x, zomplex beta, const ParComplexVector& y);
void ParComplexVecAXPBYPCZ(zomplex alpha, const ParComplexVector& x, zomplex beta, const ParComplexVector& y, zomplex gamma, const ParComplexVector& z);
zomplex ParComplexVecConjDot(const ParComplexVector& x, const ParComplexVector& y);
zomplex ParComplexVecDot(const ParComplexVector& x, const ParComplexVector& y);
void ParComplexVecElemMul(const ParComplexVector& x, const ParComplexVector& y);
void ParComplexVecElemMulConj(const ParComplexVector& x, const ParComplexVector& y);
void ParComplexVecConj(const ParComplexVector& x);
void ParComplexVecRecip(const ParComplexVector& x);

class ComplexOperator
{
public:
	virtual ~ComplexOperator();
	virtual int InSize() const = 0;
	virtual int OutSize() const = 0;
	virtual void Apply(const ComplexVector& x, const ComplexVector& y) const = 0;
};

class ParComplexOperator
{
public:
	MPI_Comm comm;

	ParComplexOperator(MPI_Comm comm = MPI_COMM_SELF);
	virtual ~ParComplexOperator();
	virtual int InSize() const = 0;
	virtual int OutSize() const = 0;
	virtual void Apply(const ParComplexVector& x, const ParComplexVector& y) const = 0;
};

class ParComplexOperatorIdent : public ParComplexOperator
{
private:
	int size[2];

public:
	ParComplexOperatorIdent(MPI_Comm comm, int n);
	ParComplexOperatorIdent(MPI_Comm comm, int n, int m);	
	int InSize() const;
	int OutSize() const;
	void Apply(const ParComplexVector& x, const ParComplexVector& y) const;
};

class ParComplexSolver
{
public:
	virtual ~ParComplexSolver();
	virtual void operator()(const ParComplexOperator& A, const ParComplexVector& b, const ParComplexVector& x, int& iter, double& relres) const;
	virtual void operator()(const ParComplexOperator& A, const ParComplexOperator& P, const ParComplexVector& b, const ParComplexVector& x, int& iter, double& relres) const = 0;
};

class ParComplexSolverCG : public ParComplexSolver
{
public:
	int    MaxIters;
	double Tolerance;
	int    PrintStats;
	
	using ParComplexSolver::operator();
	ParComplexSolverCG(int max_iters = 500, double tolerance = 1.0e-08, int print_stats = 0);
	void operator()(const ParComplexOperator& A, const ParComplexOperator& P, const ParComplexVector& b, const ParComplexVector& x, int& iter, double& relres) const;
};

class ParComplexSolverGMRES : public ParComplexSolver
{
public:
	int    Restart;
	int    MaxIters;
	double Tolerance;
	int    PrintStats;

	using ParComplexSolver::operator();
	ParComplexSolverGMRES(int restart = 10, int max_iters = 500, double tolerance = 1.0e-08, int print_stats = 0);
	void operator()(const ParComplexOperator& A, const ParComplexOperator& P, const ParComplexVector& b, const ParComplexVector& x, int& iter, double& relres) const;
};

class ParComplexSolverBCGS : public ParComplexSolver
{
public:
	int    MaxIters;
	double Tolerance;
	int    PrintStats;

	using ParComplexSolver::operator();
	ParComplexSolverBCGS(int max_iters = 500, double tolerance = 1.0e-08, int print_stats = 0);
	void operator()(const ParComplexOperator& A, const ParComplexOperator& P, const ParComplexVector& b, const ParComplexVector& x, int& iter, double& relres) const;
};

struct ComplexCOOMatrix : public ComplexOperator
{
	int ref;
	int size[2];
	int nnz;
	
	int* rowind;
	int* colind;
	zomplex* values;

	ComplexCOOMatrix();
	ComplexCOOMatrix(int n, int m, int nnz, int* rowind, int* colind, zomplex* values, int ref);
	ComplexCOOMatrix(const COOMatrix& A);
	ComplexCOOMatrix(const ComplexCOOMatrix& A);
	ComplexCOOMatrix(ComplexCOOMatrix&& A);
	~ComplexCOOMatrix();
	ComplexCOOMatrix& operator=(const ComplexCOOMatrix& A);
	ComplexCOOMatrix& operator=(ComplexCOOMatrix&& A);
	
	void Free();
	void Refer(const ComplexCOOMatrix& A);
	int InSize() const;
	int OutSize() const;
	void Apply(const ComplexVector& x, const ComplexVector& y) const;
};

void ComplexCOORead(const char* filename, ComplexCOOMatrix& A);
void ComplexCOOWrite(const char* filename, const ComplexCOOMatrix& A);

struct ComplexCSRMatrix : public ComplexOperator
{
	int ref;
	int size[2];

	int* rowptr;
	int* colind;
	zomplex* values;

	ComplexCSRMatrix();
	ComplexCSRMatrix(int n, int m, int* rowptr, int* colind, zomplex* values, int ref);
	ComplexCSRMatrix(const CSRMatrix& A);
	ComplexCSRMatrix(const ComplexCSRMatrix& A);
	ComplexCSRMatrix(const ComplexCOOMatrix& A);
	ComplexCSRMatrix(ComplexCSRMatrix&& A);
	~ComplexCSRMatrix();
	ComplexCSRMatrix& operator=(const ComplexCOOMatrix& A);
	ComplexCSRMatrix& operator=(const ComplexCSRMatrix& A);
	ComplexCSRMatrix& operator=(ComplexCSRMatrix&& A);

	void Free();
	void Refer(const ComplexCSRMatrix& A);
	int InSize() const;
	int OutSize() const;
	void Apply(const ComplexVector& x, const ComplexVector& y) const;
};

void ComplexCSRRead(const char* filename, ComplexCSRMatrix& A);
void ComplexCSRWrite(const char* filename, const ComplexCSRMatrix& A);
void ComplexCSRDiag(const ComplexCSRMatrix& A, const ComplexVector& D);
void ComplexCSREliminZeros(const ComplexCSRMatrix& A);
void ComplexCSRConj(const ComplexCSRMatrix& A);
void ComplexCSRScale(zomplex alpha, const ComplexCSRMatrix& A);
void ComplexCSRScaleRows(const ComplexVector& x, const ComplexCSRMatrix& A);
void ComplexCSRScaleCols(const ComplexVector& x, const ComplexCSRMatrix& A);
void ComplexCSRMatAdd(const ComplexCSRMatrix& A, const ComplexCSRMatrix& B, ComplexCSRMatrix& C);
void ComplexCSRMatMul(const ComplexCSRMatrix& A, const ComplexCSRMatrix& B, ComplexCSRMatrix& C);
void ComplexCSRMatVec(zomplex alpha, const ComplexCSRMatrix& A, const ComplexVector& x, zomplex beta, const ComplexVector& y);
void ComplexCSRTrans(const ComplexCSRMatrix& A, ComplexCSRMatrix& B);
void ComplexCSRConjTrans(const ComplexCSRMatrix& A, ComplexCSRMatrix& B);

class ComplexCSRCreator
{
private:
	int size[2];

public:
	ComplexCSRCreator(int n, int m);
	void operator()(ComplexCSRMatrix& A) const;
};

class ComplexCSRMutator
{
private:
	ComplexCSRMatrix A;
	bool* marks;
	void* data;

public:
	ComplexCSRMutator(const ComplexCSRMatrix& A);
	~ComplexCSRMutator();
	zomplex* Find(int row, int col) const;
	zomplex* Insert(int row, int col) const;
	void Erase(int row, int col) const;
	void operator()(ComplexCSRMatrix& A);
};

class ComplexCSRAccessor
{
private:
	ComplexCSRMatrix A;

public:
	ComplexCSRAccessor(const ComplexCSRMatrix& A);
	zomplex* Find(int row, int col) const;
};

struct ParComplexCSRMatrix : public ParComplexOperator
{
	ComplexCSRMatrix local;
	ComplexCSRMatrix exter;

	int nnb;
	int* nbrank;

	int* recvptr;
	int* recvind;

	int* sendptr;
	int* sendind;
	zomplex* sendbuf;
	
	ComplexVector recvx;

	ParComplexCSRMatrix(MPI_Comm comm = MPI_COMM_SELF);
	ParComplexCSRMatrix(MPI_Comm comm, int local_rows, int local_cols, int exter_cols, 
	int* local_rowptr, int* local_colind, zomplex* local_values, int local_ref, 
	int* exter_rowptr, int* exter_colind, zomplex* exter_values, int exter_ref,
	int nnb, int* nbrank, int* recvptr, int* recvind);
	ParComplexCSRMatrix(const ComplexCSRMatrix& A);
	ParComplexCSRMatrix(const ParCSRMatrix& A);
	ParComplexCSRMatrix(const ParComplexCSRMatrix& A);
	ParComplexCSRMatrix(ParComplexCSRMatrix&& A);
	~ParComplexCSRMatrix();
	ParComplexCSRMatrix& operator=(const ParComplexCSRMatrix& A);
	ParComplexCSRMatrix& operator=(ParComplexCSRMatrix&& A);
	
	void Free();
	void Refer(const ParComplexCSRMatrix& A);
	void ExchangeHalo(const ParComplexVector& x) const;
	void SetupHalo();
	int InSize() const;
	int OutSize() const;
	void Apply(const ParComplexVector& x, const ParComplexVector& y) const;
};

void ParComplexCSRDiag(const ParComplexCSRMatrix& A, const ParComplexVector& D);
void ParComplexCSREliminZeros(const ParComplexCSRMatrix& A);
void ParComplexCSRConj(const ParComplexCSRMatrix& A);
void ParComplexCSRScale(zomplex alpha, const ParComplexCSRMatrix& A);
void ParComplexCSRScaleRows(const ParComplexVector& x, const ParComplexCSRMatrix& A);
void ParComplexCSRScaleCols(const ParComplexVector& x, const ParComplexCSRMatrix& A);
void ParComplexCSRMatAdd(const ParComplexCSRMatrix& A, const ParComplexCSRMatrix& B, ParComplexCSRMatrix& C);
void ParComplexCSRMatMul(const ParComplexCSRMatrix& A, const ParComplexCSRMatrix& B, ParComplexCSRMatrix& C);
void ParComplexCSRMatVec(zomplex alpha, const ParComplexCSRMatrix& A, const ParComplexVector& x, zomplex beta, const ParComplexVector& y);
void ParComplexCSRTrans(const ParComplexCSRMatrix& A, ParComplexCSRMatrix& B);
void ParComplexCSRConjTrans(const ParComplexCSRMatrix& A, ParComplexCSRMatrix& B);

class ParComplexCSRCreator
{
private:
	MPI_Comm comm;
	int size[2];

public:
	ParComplexCSRCreator(MPI_Comm comm, int n, int m);
	void operator()(ParComplexCSRMatrix& A) const;
};

class ParComplexCSRMutator
{
private:
	MPI_Comm comm;
	ComplexCSRMatrix locref;
	ComplexCSRMatrix extref;
	int* recvind;
	int* recvfrom;
	bool* marks;
	void* data;

public:
	ParComplexCSRMutator(const ParComplexCSRMatrix& A);
	~ParComplexCSRMutator();
	zomplex* Find(int row, global col) const;
	zomplex* Insert(int row, global col) const;
	void Erase(int row, global col) const;
	void operator()(ParComplexCSRMatrix& A);
};

class ParComplexCSRAccessor
{
private:
	MPI_Comm comm;
	ComplexCSRMatrix locref;
	ComplexCSRMatrix extref;
	int* recvind;
	int* recvfrom;

public:
	ParComplexCSRAccessor(const ParComplexCSRMatrix& A);
	~ParComplexCSRAccessor();
	zomplex* Find(int row, global col) const;
};

class ParComplexCSRGenerator
{
public:
	MPI_Comm comm;

	ParComplexCSRGenerator(MPI_Comm comm = MPI_COMM_SELF);
	virtual ~ParComplexCSRGenerator();
	virtual void operator()(ParComplexCSRMatrix& A) const = 0;
};

class ParComplexCSRGeneratorIJ : public ParComplexCSRGenerator
{
public:
	int size[2];
	const int* rowptr;
	const global* colind;
	const zomplex* values;
	int pattern_symmetry;

	ParComplexCSRGeneratorIJ(MPI_Comm comm = MPI_COMM_SELF);
	ParComplexCSRGeneratorIJ(MPI_Comm comm, int n, int m, const int* rowptr, const global* colind, const zomplex* values, int pattern_symmetry = 0);
	void operator()(ParComplexCSRMatrix& A) const;
};

class ParComplexCSRPrecond : public ParComplexOperator
{
public:
	virtual ~ParComplexCSRPrecond();
	virtual void Setup(const ParComplexCSRMatrix& A, int REUSE = 0) = 0;
};

class ParComplexCSRPrecondJacobi : public ParComplexCSRPrecond
{
private:
	ComplexVector D_recip;

public:
	void Free();
	void Setup(const ParComplexCSRMatrix& A, int REUSE = 0);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParComplexVector& b, const ParComplexVector& x) const;
};

class ParComplexCSRPrecondSOR : public ParComplexCSRPrecond
{
private:
	ComplexVector D_local_recip;
	ComplexCSRMatrix L_local;
	ComplexCSRMatrix L_exter;
	ComplexCSRMatrix U_local;
	ComplexCSRMatrix U_exter;
	ComplexVector recvx;

	bool* interior;

	int label;
	int nnb;
	int* nblab;
	int* nbrank;

	int* recvptr;
	int* sendptr;
	int* sendind;
	zomplex* sendbuf;

public:
	int RelaxationType;
	double RelaxationFactor;

	ParComplexCSRPrecondSOR(int relaxation_type = 2, double relaxation_factor = 1.0);
	~ParComplexCSRPrecondSOR();

	void Free();
	void Setup(const ParComplexCSRMatrix& A, int REUSE = 0);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParComplexVector& b, const ParComplexVector& x) const;
};

class ParComplexCSRPrecondILU0 : public ParComplexCSRPrecond
{
private:
	ComplexVector D_local_recip;
	ComplexCSRMatrix L_local;
	ComplexCSRMatrix L_exter;
	ComplexCSRMatrix U_local;
	ComplexCSRMatrix U_exter;
	ComplexVector recvx;

	bool* interior;

	int label;
	int nnb;
	int* nblab;
	int* nbrank;

	int* recvptr;
	int* sendptr;
	int* sendind;
	zomplex* sendbuf;

public:
	ParComplexCSRPrecondILU0();
	~ParComplexCSRPrecondILU0();

	void Free();
	void Setup(const ParComplexCSRMatrix& A, int REUSE = 0);
	int InSize() const;
	int OutSize() const;
	void Apply(const ParComplexVector& b, const ParComplexVector& x) const;
};

}

#endif