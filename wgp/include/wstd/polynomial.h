/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_POLYNOMIAL_
#define _WGP_STD_POLYNOMIAL_

#include "interval.h"
#include "solver.h"
#include "equations.h"

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
		Interval result = polynomial[0];
		for (int i = 1; i <= degree; ++i) {
			result = result + polynomial[i] * pow(x, i);
		}
		return result;
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

	class UnivariablePolynomialEquation {
	public:
		UnivariablePolynomialEquation(int degree, double* polynomial, double* d_polynomial);
	public:
		int GetEquationCount();
		int GetVariableCount();
		double GetVariableEpsilon(int i);
		double GetValueEpsilon(int i, bool is_checking);
		void CalculateValue(const IntervalVector<1>& variable, IntervalVector<1>& value);
		void CalculatePartialDerivative(const IntervalVector<1>& variable, IntervalMatrix<1, 1>& value);
		int GetSplitIndex(const IntervalVector<1>& variable, int prev_split_index, double size);
		int CompareIteratePriority(const IntervalVector<1>& variable1, double size1, const IntervalVector<1>& variable2, double size2);
		bool PreIterate(IntervalVector<1>* variable, SolverIteratedResult& result, double& size);
		bool CheckFinished(const Array<SolverHeapItem<IntervalVector<1>, double>>& heap);
	private:
		int m_degree;
		double* m_polynomial;
		double* m_d_polynomial;
	};

	typedef Solver<UnivariablePolynomialEquation, IntervalVector<1>, 
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

	inline void UnivariablePolynomialEquation::CalculateValue(const IntervalVector<1>& variable, IntervalVector<1>& value) {
		value.Set(0, estimate_univariate_polynomial_interval(m_degree, m_polynomial, variable.Get(0)));
	}

	inline void UnivariablePolynomialEquation::CalculatePartialDerivative(const IntervalVector<1>& variable, IntervalMatrix<1, 1>& value) {
		*value.Get(0, 0) = estimate_univariate_polynomial_interval(m_degree - 1, m_d_polynomial, variable.Get(0));
	}

	inline int UnivariablePolynomialEquation::GetSplitIndex(const IntervalVector<1>& variable, int prev_split_index, double size) {
		return 0;
	}

	inline int UnivariablePolynomialEquation::CompareIteratePriority(const IntervalVector<1>& variable1, double size1,
		const IntervalVector<1>& variable2, double size2) {
		if (size1 < size2) {
			return -1;
		}
		if (size1 > size2) {
			return 1;
		}
		return 0;
	}

	inline bool UnivariablePolynomialEquation::PreIterate(IntervalVector<1>* variable, SolverIteratedResult& result, double& size) {
		return false;
	}

	inline bool UnivariablePolynomialEquation::CheckFinished(const Array<SolverHeapItem<IntervalVector<1>>>& heap) {
		return heap.GetCount() >= m_degree * 4;
	}


}

#endif

