/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_POLYNOMIAL_
#define _WGP_STD_POLYNOMIAL_

#include "interval.h"
#include "solver.h"
#include "equations.h"

#include <assert.h>

namespace wgp {
	
	inline void add_univariate_polynomial(int degree, double* polynomial1, double* polynomial2, double* result) {
		for (int i = 0; i <= degree; ++i) {
			result[i] = polynomial1[i] + polynomial2[i];
		}
	}

	inline void add_univariate_polynomial(int degree, double* polynomial1, double* polynomial2) {
		for (int i = 0; i <= degree; ++i) {
			polynomial1[i] += polynomial2[i];
		}
	}

	inline void sub_univariate_polynomial(int degree, double* polynomial1, double* polynomial2, double* result) {
		for (int i = 0; i <= degree; ++i) {
			result[i] = polynomial1[i] - polynomial2[i];
		}
	}

	inline void sub_univariate_polynomial(int degree, double* polynomial1, double* polynomial2) {
		for (int i = 0; i <= degree; ++i) {
			polynomial1[i] -= polynomial2[i];
		}
	}

	inline void mul_univariate_polynomial(int degree, double* polynomial1, double* polynomial2, double* result_low, double* result_high) {
		for (int i = 0; i < degree; ++i) {
			result_low[i] = 0;
			result_high[i] = 0;
		}
		result_low[degree] = 0;
		for (int i = 0; i <= degree; ++i) {
			for (int j = 0; j <= degree - i; ++j) {
				result_low[i + j] += polynomial1[i] * polynomial2[j];
			}
			for (int j = degree - i + 1; j <= degree; ++j) {
				result_high[i - j - degree] += polynomial1[i] * polynomial2[j];
			}
		}
	}

	inline void mul_univariate_polynomial(int degree1, double* polynomial1, int degree2, double* polynomial2, double* result) {
		for (int i = 0; i <= degree1 + degree2; ++i) {
			result[i] = 0;
		}
		for (int i = 0; i <= degree1; ++i) {
			for (int j = 0; j <= degree2; ++j) {
				result[i + j] += polynomial1[i] * polynomial2[j];
			}
		}
	}

	inline void add_mul_univariate_polynomial(double* polynomial1, int degree2, double* polynomial2, int degree3, double* polynomial3) {
		for (int i = 0; i <= degree2; ++i) {
			for (int j = 0; j <= degree3; ++j) {
				polynomial1[i + j] += polynomial2[i] * polynomial3[j];
			}
		}
	}

	inline void sub_mul_univariate_polynomial(double* polynomial1, int degree2, double* polynomial2, int degree3, double* polynomial3) {
		for (int i = 0; i <= degree2; ++i) {
			for (int j = 0; j <= degree3; ++j) {
				polynomial1[i + j] -= polynomial2[i] * polynomial3[j];
			}
		}
	}

	inline void add_mul_univariate_polynomial(double* polynomial1, int degree2, double* polynomial2, int degree3, double* polynomial3, double d) {
		for (int i = 0; i <= degree2; ++i) {
			for (int j = 0; j <= degree3; ++j) {
				polynomial1[i + j] += polynomial2[i] * polynomial3[j] * d;
			}
		}
	}

	inline void neg_univariate_polynomial(int degree, double* polynomial, double* result) {
		for (int i = 0; i <= degree; ++i) {
			result[i] = -polynomial[i];
		}
	}

	inline void neg_univariate_polynomial(int degree, double* polynomial) {
		for (int i = 0; i <= degree; ++i) {
			polynomial[i] = -polynomial[i];
		}
	}

	inline void mul_univariate_polynomial(int degree, double* polynomial, double d, double* result) {
		for (int i = 0; i <= degree; ++i) {
			result[i] = polynomial[i] * d;
		}
	}

	inline void mul_univariate_polynomial(int degree, double* polynomial, double d) {
		for (int i = 0; i <= degree; ++i) {
			polynomial[i] *= d;
		}
	}

	inline void add_mul_univariate_polynomial(double* polynomial, int degree2, double* polynomial2, double d) {
		for (int i = 0; i <= degree2; ++i) {
			polynomial[i] += polynomial2[i] * d;
		}
	}

	inline void div_univariate_polynomial(int degree, double* polynomial, double d, double* result) {
		for (int i = 0; i <= degree; ++i) {
			result[i] = polynomial[i] / d;
		}
	}

	inline void div_univariate_polynomial(int degree, double* polynomial, double d) {
		for (int i = 0; i <= degree; ++i) {
			polynomial[i] /= d;
		}
	}

	inline void univariate_polynomial_dx(int degree, double* polynomial, double* dx_polynomial) {
		for (int i = 1; i <= degree; ++i) {
			dx_polynomial[i - 1] = i * polynomial[i];
		}
	}

	inline void univariate_polynomial_dx(int degree, double* polynomial) {
		for (int i = 1; i <= degree; ++i) {
			polynomial[i - 1] = i * polynomial[i];
		}
	}

	inline int degree_numerator_univariate_rational_polynomial_dx(int numerator_degree, int denominator_degree) {
		return numerator_degree + denominator_degree - 1;
	}

	inline void numerator_univariate_rational_polynomial_dx(int numerator_degree, double* numerator_polynomial,
		int denominator_degree, double* denominator_polynomial, double* result) {
		int degree_result = degree_numerator_univariate_rational_polynomial_dx(numerator_degree, denominator_degree);
		for (int i = 0; i <= degree_result; ++i) {
			result[i] = 0;
		}
		for (int i = 1; i <= numerator_degree; ++i) {
			for (int j = 0; j <= denominator_degree; ++j) {
				result[i + j - 1] += i * numerator_polynomial[i] * denominator_polynomial[j];
			}
		}
		for (int i = 0; i <= numerator_degree; ++i) {
			for (int j = 1; j <= denominator_degree; ++j) {
				result[i + j - 1] -= j * numerator_polynomial[i] * denominator_polynomial[j];
			}
		}
	}

	inline int degree_denominator_univariate_rational_polynomial_dx(int denominator_degree) {
		return denominator_degree * 2;
	}

	inline void denominator_univariate_rational_polynomial_dx(int denominator_degree, double* denominator_polynomial, double* result) {
		int degree_result = degree_denominator_univariate_rational_polynomial_dx(denominator_degree);
		for (int i = 0; i <= degree_result; ++i) {
			result[i] = 0;
		}
		for (int i = 0; i <= denominator_degree; ++i) {
			for (int j = 0; j <= denominator_degree; ++j) {
				result[i + j] += denominator_polynomial[i] * denominator_polynomial[j];
			}
		}
	}

	inline double calculate_univariate_polynomial_value(int degree, double* polynomial, double x) {
		double result = polynomial[0];
		double d = 1;
		for (int i = 1; i <= degree; ++i) {
			d *= x;
			result += polynomial[i] * d;
		}
		return result;
	}

	inline Interval estimate_univariate_polynomial_interval(int degree, double* polynomial, const Interval& x) {
		if (x.Max == x.Min) {
			return calculate_univariate_polynomial_value(degree, polynomial, x.Min);
		}
		else {
			Interval result = polynomial[0];
			double d0 = 1;
			double d1 = 1;
			if (x.Max * x.Min < 0) {
				if (abs(x.Max) < abs(x.Min)) {
					for (int i = 1; i <= degree; ++i) {
						d0 *= x.Min;
						d1 *= x.Max;
						if (i & 1) {
							result = result + polynomial[i] * Interval(d0, d1);
						}
						else {
							result = result + polynomial[i] * Interval(0, d0);
						}
					}
				}
				else {
					for (int i = 1; i <= degree; ++i) {
						d0 *= x.Min;
						d1 *= x.Max;
						if (i & 1) {
							result = result + polynomial[i] * Interval(d0, d1);
						}
						else {
							result = result + polynomial[i] * Interval(0, d1);
						}
					}
				}
			}
			else {
				if (x.Max > 0) {
					for (int i = 1; i <= degree; ++i) {
						d0 *= x.Min;
						d1 *= x.Max;
						result = result + polynomial[i] * Interval(d0, d1);
					}
				}
				else {
					for (int i = 1; i <= degree; ++i) {
						d0 *= x.Min;
						d1 *= x.Max;
						if (i & 1) {
							result = result + polynomial[i] * Interval(d0, d1);
						}
						else {
							result = result + polynomial[i] * Interval(d1, d0);
						}
					}
				}
			}
			return result;
		}
	}

	inline double calculate_univariate_rational_polynomial_value(int numerator_degree, double* numerator_polynomial,
		int denominator_degree, double* denominator_polynomial, double x) {
		return calculate_univariate_polynomial_value(numerator_degree, numerator_polynomial, x) /
			calculate_univariate_polynomial_value(denominator_degree, denominator_polynomial, x);
	}

	inline Interval estimate_univariate_rational_polynomial_interval(int numerator_degree, double* numerator_polynomial,
		int denominator_degree, double* denominator_polynomial, const Interval& x) {
		return estimate_univariate_polynomial_interval(numerator_degree, numerator_polynomial, x) /
			estimate_univariate_polynomial_interval(denominator_degree, denominator_polynomial, x);
	}

	class UnivariablePolynomialVariable {
	public:
		UnivariablePolynomialVariable() : m_d0_dirty(true), m_dt_dirty(true) {
		}

		UnivariablePolynomialVariable(int degree) : UnivariablePolynomialVariable() {
		}

		Interval Get(int index) const {
			return m_t;
		}

		void Set(int index, const Interval& t) {
			m_t = t;
		}

		void Split(int index, UnivariablePolynomialVariable& vt1, UnivariablePolynomialVariable& vt2) const {
			double m = m_t.Center();
			vt1.m_t = Interval(m_t.Min, m);
			vt1.m_d0_dirty = true;
			vt1.m_dt_dirty = true;
			vt2.m_t = Interval(m, m_t.Max);
			vt2.m_d0_dirty = true;
			vt2.m_dt_dirty = true;
		}
	private:
		Interval m_t;
		mutable bool m_d0_dirty;
		mutable Interval m_d0;
		mutable bool m_dt_dirty;
		mutable Interval m_dt;
	private:
		friend class UnivariablePolynomialEquation;
	};

	class UnivariablePolynomialEquation {
	public:
		UnivariablePolynomialEquation(int degree, double* polynomial, double* d_polynomial);
	public:
		int GetEquationCount();
		int GetVariableCount();
		double GetVariableEpsilon(int i);
		double GetValueEpsilon(int i, bool is_checking);
		void CalculateValue(const UnivariablePolynomialVariable& variable, IntervalVector<1>& value);
		void CalculatePartialDerivative(const UnivariablePolynomialVariable& variable, IntervalMatrix<1, 1>& value);
		int GetSplitIndex(const UnivariablePolynomialVariable& variable, int prev_split_index, double size);
		int CompareIteratePriority(const UnivariablePolynomialVariable& variable1, double size1, const UnivariablePolynomialVariable& variable2, double size2);
		bool PreIterate(UnivariablePolynomialVariable* variable, SolverIteratedResult& result, double& size);
		bool CheckFinished(const Array<SolverHeapItem<UnivariablePolynomialVariable, double>>& heap);
	private:
		Interval CalculateD0(const UnivariablePolynomialVariable& variable);
		Interval CalculateDt(const UnivariablePolynomialVariable& variable);
	private:
		int m_degree;
		double* m_polynomial;
		double* m_d_polynomial;
	};

	typedef Solver<UnivariablePolynomialEquation, UnivariablePolynomialVariable,
		IntervalVector<1>, IntervalVector<1>, IntervalMatrix<1, 1>> UnivariablePolynomialEquationSolver;


	inline UnivariablePolynomialEquation::UnivariablePolynomialEquation(int degree, double* polynomial, double* d_polynomial) :
		m_degree(degree),
		m_polynomial(polynomial),
		m_d_polynomial(d_polynomial) {
	}

	inline int UnivariablePolynomialEquation::GetEquationCount() {
		return 1;
	}

	inline int UnivariablePolynomialEquation::GetVariableCount() {
		return 1;
	}

	inline double UnivariablePolynomialEquation::GetVariableEpsilon(int i) {
		return 1E-6;
	}

	inline double UnivariablePolynomialEquation::GetValueEpsilon(int i, bool is_checking) {
		return is_checking ? 1E-6 : 1E-12;
	}

	inline void UnivariablePolynomialEquation::CalculateValue(const UnivariablePolynomialVariable& variable, IntervalVector<1>& value) {		
		value.Set(0, CalculateD0(variable));
	}

	inline void UnivariablePolynomialEquation::CalculatePartialDerivative(const UnivariablePolynomialVariable& variable, IntervalMatrix<1, 1>& value) {		
		*value.Get(0, 0) = CalculateDt(variable);
	}

	inline int UnivariablePolynomialEquation::GetSplitIndex(const UnivariablePolynomialVariable& variable, int prev_split_index, double size) {
		return 0;
	}

	inline int UnivariablePolynomialEquation::CompareIteratePriority(const UnivariablePolynomialVariable& variable1, double size1,
		const UnivariablePolynomialVariable& variable2, double size2) {
		if (size1 < size2) {
			return -1;
		}
		if (size1 > size2) {
			return 1;
		}
		return 0;
	}

	inline bool UnivariablePolynomialEquation::PreIterate(UnivariablePolynomialVariable* variable, SolverIteratedResult& result, double& size) {
		return false;
	}

	inline bool UnivariablePolynomialEquation::CheckFinished(const Array<SolverHeapItem<UnivariablePolynomialVariable>>& heap) {
		return heap.GetCount() >= m_degree;
	}

	inline Interval UnivariablePolynomialEquation::CalculateD0(const UnivariablePolynomialVariable& variable) {
		if (variable.m_d0_dirty) {
			variable.m_d0 = estimate_univariate_polynomial_interval(m_degree, m_polynomial, variable.m_t);
		}
		return variable.m_d0;
	}

	inline Interval UnivariablePolynomialEquation::CalculateDt(const UnivariablePolynomialVariable& variable) {
		if (variable.m_dt_dirty) {
			variable.m_dt = estimate_univariate_polynomial_interval(m_degree - 1, m_d_polynomial, variable.m_t);
		}
		return variable.m_dt;
	}


}

#endif

