/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_EQUATIONS_
#define _WGP_STD_EQUATIONS_

#pragma warning(push)
#pragma warning(disable:26495)

#include "interval.h"

namespace wgp {

	class WGP_API StandardIntervalVector {
	public:
		StandardIntervalVector();
		StandardIntervalVector(int degree);
		StandardIntervalVector(const StandardIntervalVector& vt);
		virtual ~StandardIntervalVector();
		int GetDegree() const;
		StandardIntervalVector& operator=(const StandardIntervalVector& vt);
		Interval* Get(int i);
		const Interval* Get(int i) const;
	public:
		StandardIntervalVector Center() const;
		void Center(StandardIntervalVector& vt) const;
		StandardIntervalVector Min() const;
		void Min(StandardIntervalVector& vt) const;
		StandardIntervalVector Max() const;
		void Max(StandardIntervalVector& vt) const;
	private:
		int m_degree;
		Interval* m_data;
	};

	class WGP_API StandardIntervalMatrix {
	public:
		StandardIntervalMatrix();
		StandardIntervalMatrix(int row_count, int col_count);
		StandardIntervalMatrix(const StandardIntervalMatrix& matrix);
		virtual ~StandardIntervalMatrix();
		int GetRowCount() const;
		int GetColCount() const;
		StandardIntervalMatrix& operator=(const StandardIntervalMatrix& matrix);
		Interval* Get(int i, int j);
		const Interval* Get(int i, int j) const;
	public:
		StandardIntervalVector Row(int i) const;
		void Row(int i, StandardIntervalVector& vt) const;
	private:
		int m_row_count;
		int m_col_count;
		Interval* m_data;
	};

	class WGP_API StandardMatrix {
	public:
		StandardMatrix();
		StandardMatrix(int row_count, int col_count);
		StandardMatrix(const StandardMatrix& matrix);
		virtual ~StandardMatrix();
		int GetRowCount() const;
		int GetColCount() const;
		StandardMatrix& operator=(const StandardMatrix& matrix);
		double* Get(int i, int j);
		const double* Get(int i, int j) const;
	public:
		void GaussianElimination(StandardMatrix& adjoint);
	private:
		int m_row_count;
		int m_col_count;
		double* m_data;
	};

	class WGP_API StandardEquationSystem {
	public:
		virtual int GetEquationCount() = 0;
		virtual int GetVariableCount() = 0;
		virtual double GetVariableEpsilon(int i) = 0;
		virtual double GetValueEpsilon(int i) = 0;
		virtual void CalculateValue(const StandardIntervalVector& variable, StandardIntervalVector& value) = 0;
		virtual void CalculatePartialDerivative(const StandardIntervalVector& variable, StandardIntervalMatrix& value) = 0;
	};

}

#pragma warning(pop)

#endif