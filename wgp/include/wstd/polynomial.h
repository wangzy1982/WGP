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
	
	inline void add_univariate_polynomial(int highest_degree, double* polynomial1, double* polynomial2, double* result) {
		for (int i = 0; i <= highest_degree; ++i) {
			result[i] = polynomial1[i] + polynomial2[i];
		}
	}

	inline void add_univariate_polynomial(int highest_degree, double* polynomial1, double* polynomial2) {
		for (int i = 0; i <= highest_degree; ++i) {
			polynomial1[i] += polynomial2[i];
		}
	}

	inline void sub_univariate_polynomial(int highest_degree, double* polynomial1, double* polynomial2, double* result) {
		for (int i = 0; i <= highest_degree; ++i) {
			result[i] = polynomial1[i] - polynomial2[i];
		}
	}

	inline void sub_univariate_polynomial(int highest_degree, double* polynomial1, double* polynomial2) {
		for (int i = 0; i <= highest_degree; ++i) {
			polynomial1[i] -= polynomial2[i];
		}
	}

	inline void mul_univariate_polynomial(int highest_degree, double* polynomial1, double* polynomial2, double* result_low, double* result_high) {
		for (int i = 0; i < highest_degree; ++i) {
			result_low[i] = 0;
			result_high[i] = 0;
		}
		result_low[highest_degree] = 0;
		for (int i = 0; i <= highest_degree; ++i) {
			for (int j = 0; j <= highest_degree - i; ++j) {
				result_low[i + j] += polynomial1[i] * polynomial2[j];
			}
			for (int j = highest_degree - i + 1; j <= highest_degree; ++j) {
				result_high[i - j - highest_degree] += polynomial1[i] * polynomial2[j];
			}
		}
	}

	inline void mul_univariate_polynomial(int highest_degree1, double* polynomial1, int highest_degree2, double* polynomial2, double* result) {
		for (int i = 0; i <= highest_degree1 + highest_degree2; ++i) {
			result[i] = 0;
		}
		for (int i = 0; i <= highest_degree1; ++i) {
			for (int j = 0; j <= highest_degree2; ++j) {
				result[i + j] += polynomial1[i] * polynomial2[j];
			}
		}
	}

	inline void add_mul_univariate_polynomial(double* polynomial1, int highest_degree2, double* polynomial2, int highest_degree3, double* polynomial3) {
		for (int i = 0; i <= highest_degree2; ++i) {
			for (int j = 0; j <= highest_degree3; ++j) {
				polynomial1[i + j] += polynomial2[i] * polynomial3[j];
			}
		}
	}

	inline void neg_univariate_polynomial(int highest_degree, double* polynomial, double* result) {
		for (int i = 0; i <= highest_degree; ++i) {
			result[i] = -polynomial[i];
		}
	}

	inline void neg_univariate_polynomial(int highest_degree, double* polynomial) {
		for (int i = 0; i <= highest_degree; ++i) {
			polynomial[i] = -polynomial[i];
		}
	}

	inline void mul_univariate_polynomial(int highest_degree, double* polynomial, double d, double* result) {
		for (int i = 0; i <= highest_degree; ++i) {
			result[i] = polynomial[i] * d;
		}
	}

	inline void mul_univariate_polynomial(int highest_degree, double* polynomial, double d) {
		for (int i = 0; i <= highest_degree; ++i) {
			polynomial[i] *= d;
		}
	}

	inline void add_mul_univariate_polynomial(double* polynomial, int highest_degree2, double* polynomial2, double d) {
		for (int i = 0; i <= highest_degree2; ++i) {
			polynomial[i] += polynomial2[i] * d;
		}
	}

	inline void div_univariate_polynomial(int highest_degree, double* polynomial, double d, double* result) {
		for (int i = 0; i <= highest_degree; ++i) {
			result[i] = polynomial[i] / d;
		}
	}

	inline void div_univariate_polynomial(int highest_degree, double* polynomial, double d) {
		for (int i = 0; i <= highest_degree; ++i) {
			polynomial[i] /= d;
		}
	}

	inline void univariate_polynomial_dx(int highest_degree, double* polynomial, double* dx_polynomial) {
		for (int i = 1; i <= highest_degree; ++i) {
			dx_polynomial[i - 1] = i * polynomial[i];
		}
	}

	inline void univariate_polynomial_dx(int highest_degree, double* polynomial) {
		for (int i = 1; i <= highest_degree; ++i) {
			polynomial[i - 1] = i * polynomial[i];
		}
	}

	inline Interval calculate_univariate_polynomial_value(int highest_degree, double* polynomial, double* temp_polynomial, double* temp_powers, const Interval& x) {
		double xc = x.Center();
		temp_powers[0] = 1;
		for (int i = 1; i <= highest_degree; ++i) {
			temp_powers[i] = temp_powers[i - 1] * xc;
		}
		double d = 0;
		for (int i = 0; i <= highest_degree; ++i) {
			d += polynomial[i] * temp_powers[i];
		}
		temp_polynomial[0] = d;
		for (int i = 1; i <= highest_degree; ++i) {
			temp_polynomial[i] = polynomial[i] * i;
		}
		double c = 1;
		for (int k = 1; k <= highest_degree; ++k) {
			d = 0;
			for (int i = k; i <= highest_degree; ++i) {
				d += temp_polynomial[i] * temp_powers[i - k];
			}
			temp_polynomial[k] = d / c;
			for (int i = k + 1; i <= highest_degree; ++i) {
				temp_polynomial[i] = temp_polynomial[i] * (i - k);
			}
			c *= (k + 1);
		}
		Interval x2 = x - xc;
		Interval y = temp_polynomial[highest_degree];
		for (int i = highest_degree - 1; i >= 0; --i) {
			y = y * x2 + temp_polynomial[i];
		}
		return y;
	}

	inline double calculate_univariate_polynomial_value(int highest_degree, double* polynomial, double x) {
		double p = 1;
		double d = 0;
		for (int i = 0; i <= highest_degree; ++i) {
			d += polynomial[i] * p;
			p *= x;
		}
		return d;
	}

	inline Interval calculate_univariate_polynomial_value(int highest_degree, double* polynomial, const Interval& x) {
		if (highest_degree < 16) {
			double temp_polynomial[16];
			double temp_powers[16];
			return calculate_univariate_polynomial_value(highest_degree, polynomial, temp_polynomial, temp_powers, x);
		}
		else {
			double* temp_memory = new double[(highest_degree + 1) * 2];
			Interval y = calculate_univariate_polynomial_value(highest_degree, polynomial, temp_memory, temp_memory + (highest_degree + 1), x);
			delete[] temp_memory;
			return y;
		}
	}

	class UnivariablePolynomialEquationSystem {
	public:
		UnivariablePolynomialEquationSystem(int highest_degree, double* polynomial, double* d_polynomial) :
			m_highest_degree(highest_degree),
			m_polynomial(polynomial),
			m_d_polynomial(d_polynomial) {
		}
	public:
		int GetEquationCount() { return 1; }
		int GetVariableCount() { return 1; }
		double GetVariableEpsilon(int i) { return 1E-6; }
		double GetValueEpsilon(int i) { return 1E-6; }
		void CalculateValue(const IntervalVector<1>& variable, IntervalVector<1>& value) {
			value.Set(0, calculate_univariate_polynomial_value(m_highest_degree, m_polynomial, variable.Get(0)));
		}
		void CalculatePartialDerivative(const IntervalVector<1>& variable, IntervalMatrix<1, 1>& value) {
			*value.Get(0, 0) = calculate_univariate_polynomial_value(m_highest_degree - 1, m_d_polynomial, variable.Get(0));
		}
		void Transform(const IntervalVector<1>& variable, IntervalVector<1>& value, IntervalMatrix<1, 1>& partial_derivative,
			bool& recheck_value, bool& use_default_transform) {
			recheck_value = false;
			use_default_transform = false;
		}
		void Restore() {}
		int GetSplitIndex(const IntervalVector<1>& variable, int prev_split_index, double size) { return 0; }
		int CompareIteratePriority(const IntervalVector<1>& variable1, double size1, const IntervalVector<1>& variable2, double size2) {
			if (size1 < size2) {
				return -1;
			}
			if (size1 > size2) {
				return 1;
			}
			return 0;
		}
		bool SpeciallySolve(IntervalVector<1>* variable, SolverIteratedResult& result, double& size) {
			return false;
		}
	private:
		int m_highest_degree;
		double* m_polynomial;
		double* m_d_polynomial;
	};

	inline Interval calculate_univariate_polynomial_value_ex(int highest_degree, double* polynomial, double* d_polynomial, const Interval& x) {
		Solver<UnivariablePolynomialEquationSystem, IntervalVector<1>, IntervalVector<1>, 
			IntervalVector<1>, IntervalMatrix<1, 1>, Matrix<1, 1>> solver;
		UnivariablePolynomialEquationSystem equations(highest_degree, polynomial, d_polynomial);
		solver.SetEquationSystem(&equations);
		solver.SetMaxFuzzyRootCount(10);
		IntervalVector<1> initial_variable;
		initial_variable.Set(0, x);
		solver.SetInitialVariable(initial_variable);
		const Array<IntervalVector<1>>& clear_roots = solver.GetClearRoots();
		const Array<IntervalVector<1>>& fuzzy_roots = solver.GetFuzzyRoots();
		Interval result = calculate_univariate_polynomial_value(highest_degree, polynomial, x.Min);
		result.Merge(calculate_univariate_polynomial_value(highest_degree, polynomial, x.Max));
		for (int i = 0; i < clear_roots.GetCount(); ++i) {
			result.Merge(calculate_univariate_polynomial_value(highest_degree, polynomial, clear_roots.Get(i).Get(0)));
		}
		for (int i = 0; i < fuzzy_roots.GetCount(); ++i) {
			result.Merge(calculate_univariate_polynomial_value(highest_degree, polynomial, fuzzy_roots.Get(i).Get(0)));
		}
		return result;
	}
}

#endif

