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
		virtual ~StandardEquationSystem() {}
	public:
		virtual int GetEquationCount() = 0;
		virtual int GetVariableCount() = 0;
		virtual double GetVariableEpsilon(int i) = 0;
		virtual double GetValueEpsilon(int i) = 0;
		virtual void CalculateValue(const StandardIntervalVector& variable, StandardIntervalVector& value) = 0;
		virtual void CalculatePartialDerivative(const StandardIntervalVector& variable, StandardIntervalMatrix& value) = 0;
		virtual bool Transform(const StandardIntervalVector& variable, StandardIntervalMatrix& partial_derivative);
		virtual void Restore();
	};

	template<int degree>
	class WGP_API IntervalVector {
	public:
		IntervalVector() {}
		IntervalVector(int degree) {}

		virtual ~IntervalVector() {}

		int GetDegree() const { return degree; }

		Interval* Get(int i) {
			return m_data + i;
		}

		const Interval* Get(int i) const {
			return m_data + i;
		}
	public:
		IntervalVector Center() const {
			IntervalVector vt;
			for (int i = 0; i < degree; ++i) {
				vt.m_data[i] = m_data[i].Center();
			}
			return vt;
		}

		void Center(IntervalVector& vt) const {
			for (int i = 0; i < degree; ++i) {
				vt.m_data[i] = m_data[i].Center();
			}
		}

		IntervalVector Min() const {
			IntervalVector vt;
			for (int i = 0; i < degree; ++i) {
				vt.m_data[i] = m_data[i].Min;
			}
			return vt;
		}

		void Min(IntervalVector& vt) const {
			for (int i = 0; i < degree; ++i) {
				vt.m_data[i] = m_data[i].Min;
			}
		}

		IntervalVector Max() const {
			IntervalVector vt;
			for (int i = 0; i < degree; ++i) {
				vt.m_data[i] = m_data[i].Max;
			}
			return vt;
		}

		void Max(IntervalVector& vt) const {
			for (int i = 0; i < degree; ++i) {
				vt.m_data[i] = m_data[i].Max;
			}
		}
	private:
		Interval m_data[degree];
	};

	template<int row_count, int col_count>
	class WGP_API IntervalMatrix {
	public:
		IntervalMatrix() {}

		IntervalMatrix(int row_count, int col_count) {}

		virtual ~IntervalMatrix() {}

		int GetRowCount() const { return row_count; }

		int GetColCount() const { return col_count; }

		Interval* Get(int i, int j) {
			return &m_data[i][j];
		}

		const Interval* Get(int i, int j) const {
			return &m_data[i][j];
		}

	public:
		IntervalVector<col_count> Row(int i) const {
			IntervalVector<col_count> vt;
			for (int j = 0; j < col_count; ++j) {
				*vt.Get(j) = m_data[i][j];
			}
			return vt;
		}

		void Row(int i, IntervalVector<col_count>& vt) const {
			for (int j = 0; j < col_count; ++j) {
				*vt.Get(j) = m_data[i][j];
			}
		}
	private:
		Interval m_data[row_count][col_count];
	};

	template<int row_count, int col_count>
	class WGP_API Matrix {
	public:
		Matrix() {}

		Matrix(int row_count, int col_count) {}

		virtual ~Matrix() {}

		int GetRowCount() const { return row_count; }

		int GetColCount() const { return col_count; }

		double* Get(int i, int j) {
			return &m_data[i][j];
		}

		const double* Get(int i, int j) const {
			return &m_data[i][j];
		}

	public:
		template<int adjoint_col_count>
		void GaussianElimination(Matrix<col_count, adjoint_col_count>& adjoint) {
			for (int i = 0; i < row_count - 1; ++i) {
				int k = i;
				double d = abs(*Get(k, i));
				for (int j = i + 1; j < row_count; ++j) {
					double d2 = abs(*Get(j, i));
					if (d2 > d) {
						d = d2;
						k = j;
					}
				}
				if (d > g_double_epsilon) {
					for (int j = i; j < row_count; ++j) {
						double t = *Get(i, j);
						*Get(i, j) = *Get(k, j);
						*Get(k, j) = t;
					}
					for (int j = 0; j < adjoint_col_count; ++j) {
						double t = *adjoint.Get(i, j);
						*adjoint.Get(i, j) = *adjoint.Get(k, j);
						*adjoint.Get(k, j) = t;
					}
					for (int k = i + 1; k < row_count; ++k) {
						double c = *Get(k, i) / *Get(i, i);
						*Get(k, i) = 0;
						for (int j = i + 1; j < row_count; ++j) {
							*Get(k, j) -= *Get(i, j) * c;
						}
						for (int j = 0; j < adjoint_col_count; ++j) {
							*adjoint.Get(k, j) -= *adjoint.Get(i, j) * c;
						}
					}
				}
			}
			for (int i = row_count - 1; i >= 0; --i) {
				if (abs(*Get(i, i)) > g_double_epsilon) {
					for (int k = 0; k < i; ++k) {
						double c = *Get(k, i) / *Get(i, i);
						*Get(k, i) = 0;
						for (int j = 0; j < adjoint_col_count; ++j) {
							*adjoint.Get(k, j) -= *adjoint.Get(i, j) * c;
						}
					}
					double c = 1 / *Get(i, i);
					*Get(i, i) = 1;
					for (int j = 0; j < adjoint_col_count; ++j) {
						*adjoint.Get(i, j) *= c;
					}
				}
			}
		}
	private:
		double m_data[row_count][col_count];
	};

	template<int equation_count, int variable_count>
	class WGP_API EquationSystem {
	public:
		virtual ~EquationSystem() {}
	public:
		int GetEquationCount() { return equation_count; }
		int GetVariableCount() { return variable_count; }
		virtual double GetVariableEpsilon(int i) = 0;
		virtual double GetValueEpsilon(int i) = 0;
		virtual void CalculateValue(const IntervalVector<variable_count>& variable, 
			IntervalVector<equation_count>& value) = 0;
		virtual void CalculatePartialDerivative(const IntervalVector<variable_count>& variable, 
			IntervalMatrix<equation_count, variable_count>& value) = 0;
		virtual bool Transform(const IntervalVector<variable_count>& variable, 
			IntervalMatrix<equation_count, variable_count>& partial_derivative) { return true; }
		virtual void Restore() {}
	};

}

#pragma warning(pop)

#endif