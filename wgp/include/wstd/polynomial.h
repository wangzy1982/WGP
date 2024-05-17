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

	inline void univariate_polynomial_dt(int degree, double* polynomial, double* dt_polynomial) {
		for (int i = 1; i <= degree; ++i) {
			dt_polynomial[i - 1] = i * polynomial[i];
		}
	}

	inline void univariate_polynomial_dt(int degree, double* polynomial) {
		for (int i = 1; i <= degree; ++i) {
			polynomial[i - 1] = i * polynomial[i];
		}
	}

	inline double calculate_univariate_polynomial_value(int degree, double* polynomial, double t) {
		double result = polynomial[0];
		double d = 1;
		for (int i = 1; i <= degree; ++i) {
			d *= t;
			result += polynomial[i] * d;
		}
		return result;
	}

	inline void mul_two_univariate_polynomial(int u_degree, double* u_polynomial, int v_degree, double* v_polynomial, double* polynomial) {
		for (int i = 0; i <= u_degree; ++i) {
			double* p = polynomial + i * v_degree;
			for (int j = 0; j <= v_degree; ++j) {
				p[j] = u_polynomial[i] * v_polynomial[j];
			}
		}
	}

	inline void mul_two_univariate_polynomial(int u_degree, double* u_polynomial, int v_degree, double* v_polynomial, double d, double* polynomial) {
		for (int i = 0; i <= u_degree; ++i) {
			double* p = polynomial + i * v_degree;
			for (int j = 0; j <= v_degree; ++j) {
				p[j] = u_polynomial[i] * v_polynomial[j] * d;
			}
		}
	}

	inline void add_mul_two_univariate_polynomial(double* polynomial, int u_degree, double* u_polynomial, int v_degree, double* v_polynomial) {
		for (int i = 0; i <= u_degree; ++i) {
			double* p = polynomial + i * v_degree;
			for (int j = 0; j <= v_degree; ++j) {
				p[j] += u_polynomial[i] * v_polynomial[j];
			}
		}
	}

	inline void add_mul_two_univariate_polynomial(double* polynomial, int u_degree, double* u_polynomial, int v_degree, double* v_polynomial, double d) {
		for (int i = 0; i <= u_degree; ++i) {
			double* p = polynomial + i * v_degree;
			for (int j = 0; j <= v_degree; ++j) {
				p[j] += u_polynomial[i] * v_polynomial[j] * d;
			}
		}
	}

	inline void two_polynomial_du(int u_degree, int v_degree, double* polynomial, double* du_polynomial) {
		for (int i = 0; i < u_degree; ++i) {
			double* p1 = polynomial + (i + 1) * v_degree;
			double* p2 = du_polynomial + i * v_degree;
			for (int j = 0; j <= v_degree; ++j) {
				p2[j] = (i + 1) * p1[j];
			}
		}
	}

	inline void two_polynomial_dv(int u_degree, int v_degree, double* polynomial, double* dv_polynomial) {
		for (int i = 0; i <= u_degree; ++i) {
			double* p1 = polynomial + i * v_degree;
			double* p2 = dv_polynomial + i * v_degree;
			for (int j = 1; j <= v_degree; ++j) {
				p2[j - 1] = j * p1[j];
			}
			p2[v_degree] = 0;
		}
	}

	inline double calculate_two_polynomial_value(int u_degree, int v_degree, double* polynomial, double u, double v) {
		double result = 0;
		double du = 1;
		for (int i = 0; i <= u_degree; ++i) {
			double* p = polynomial + i * v_degree;
			double dv = 1;
			for (int j = 0; j <= v_degree; ++j) {
				result += p[j] * du * dv;
				dv *= v;
			}
			du *= u;
		}
		return result;
	}
	
	inline Interval estimate_univariate_polynomial_interval(int degree, double* polynomial, const Interval& t) {
		if (t.Max == t.Min) {
			return calculate_univariate_polynomial_value(degree, polynomial, t.Min);
		}
		else {
			Interval result = polynomial[0];
			double d0 = 1;
			double d1 = 1;
			if (t.Max * t.Min < 0) {
				if (abs(t.Max) < abs(t.Min)) {
					for (int i = 1; i <= degree; ++i) {
						d0 *= t.Min;
						d1 *= t.Max;
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
						d0 *= t.Min;
						d1 *= t.Max;
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
				if (t.Max > 0) {
					for (int i = 1; i <= degree; ++i) {
						d0 *= t.Min;
						d1 *= t.Max;
						result = result + polynomial[i] * Interval(d0, d1);
					}
				}
				else {
					for (int i = 1; i <= degree; ++i) {
						d0 *= t.Min;
						d1 *= t.Max;
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

	inline Interval estimate_two_polynomial_interval(int u_degree, int v_degree, double* polynomial, const Interval& u, const Interval& v) {
		if (u.Min == u.Max) {
			if (v.Min == v.Max) {
				double result = 0;
				double du = 1;
				for (int i = 0; i <= u_degree; ++i) {
					double* p = polynomial + i * v_degree;
					double dv = 1;
					for (int j = 0; j <= v_degree; ++j) {
						result += p[j] * du * dv;
						dv *= v.Min;
					}
					du *= u.Min;
				}
				return result;
			}
			else {
				bool bv = v.Min * v.Max < 0;
				Interval result = 0;
				double du = 1;
				for (int i = 0; i <= u_degree; ++i) {
					double* p = polynomial + i * v_degree;
					double dv0 = 1;
					double dv1 = 1;
					for (int j = 0; j <= v_degree; ++j) {
						if (bv) {
							if (j & 1) {
								result = result + p[j] * du * Interval(dv0, dv1);
							}
							else if (dv0 < dv1) {
								result = result + p[j] * du * Interval(0, dv1);
							}
							else {
								result = result + p[j] * du * Interval(0, dv0);
							}
						} 
						else {
							if (dv0 < dv1) {
								result = result + p[j] * du * Interval(dv0, dv1);
							}
							else {
								result = result + p[j] * du * Interval(dv1, dv0);
							}
						}
						dv0 *= v.Min;
						dv1 *= v.Max;
					}
					du *= u.Min;
				}
				return result;
			}
		}
		else {
			if (v.Min == v.Max) {
				bool bu = u.Min * u.Max < 0;
				Interval result = 0;
				double du0 = 1;
				double du1 = 1;
				Interval du;
				for (int i = 0; i <= u_degree; ++i) {
					if (bu) {
						if (i & 1) {
							du.Min = du0;
							du.Max = du1;
						}
						else if (du0 > du1) {
							du.Min = 0;
							du.Max = du0;
						}
						else {
							du.Min = 0;
							du.Max = du1;
						}
					}
					else {
						if (du0 > du1) {
							du.Min = du1;
							du.Max = du0;
						}
						else {
							du.Min = du0;
							du.Max = du1;
						}
					}
					double* p = polynomial + i * v_degree;
					double dv = 1;
					for (int j = 0; j <= v_degree; ++j) {
						result = result + p[j] * dv * du;
						dv *= v.Min;
					}
					du0 *= u.Min;
					du1 *= u.Max;
				}
				return result;
			}
			else {
				bool bu = u.Min * u.Max < 0;
				bool bv = v.Min * v.Max < 0;
				Interval result = 0;
				double du0 = 1;
				double du1 = 1;
				Interval du;
				for (int i = 0; i <= u_degree; ++i) {
					if (bu) {
						if (i & 1) {
							du.Min = du0;
							du.Max = du1;
						}
						else if (du0 > du1) {
							du.Min = 0;
							du.Max = du0;
						}
						else {
							du.Min = 0;
							du.Max = du1;
						}
					}
					else {
						if (du0 > du1) {
							du.Min = du1;
							du.Max = du0;
						}
						else {
							du.Min = du0;
							du.Max = du1;
						}
					}
					double* p = polynomial + i * v_degree;
					double dv0 = 1;
					double dv1 = 1;
					for (int j = 0; j <= v_degree; ++j) {
						if (bv) {
							if (j & 1) {
								result = result + p[j] * du * Interval(dv0, dv1);
							}
							else if (dv0 < dv1) {
								result = result + p[j] * du * Interval(0, dv1);
							}
							else {
								result = result + p[j] * du * Interval(0, dv0);
							}
						}
						else {
							if (dv0 < dv1) {
								result = result + p[j] * du * Interval(dv0, dv1);
							}
							else {
								result = result + p[j] * du * Interval(dv1, dv0);
							}
						}
						dv0 *= v.Min;
						dv1 *= v.Max;
					}
					du0 *= u.Min;
					du1 *= u.Max;
				}
				return result;
			}
		}
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
		return heap.GetCount() >= m_degree;
	}
	
	class TwoPolynomialEquation {
	public:
		TwoPolynomialEquation(int u_degree1, int v_degree1, double* polynomial1, double* du_polynomial1, double* dv_polynomial1, 
			int u_degree2, int v_degree2, double* polynomial2, double* du_polynomial2, double* dv_polynomial2);
	public:
		int GetEquationCount();
		int GetVariableCount();
		double GetVariableEpsilon(int i);
		double GetValueEpsilon(int i, bool is_checking);
		void CalculateValue(const IntervalVector<2>& variable, IntervalVector<2>& value);
		void CalculatePartialDerivative(const IntervalVector<2>& variable, IntervalMatrix<2, 2>& value);
		int GetSplitIndex(const IntervalVector<2>& variable, int prev_split_index, double size);
		int CompareIteratePriority(const IntervalVector<2>& variable1, double size1, const IntervalVector<2>& variable2, double size2);
		bool PreIterate(IntervalVector<2>* variable, SolverIteratedResult& result, double& size);
		bool CheckFinished(const Array<SolverHeapItem<IntervalVector<2>, double>>& heap);
	private:
		int m_u_degree0;
		int m_v_degree0;
		double* m_polynomial0;
		double* m_du_polynomial0;
		double* m_dv_polynomial0;
		int m_u_degree1;
		int m_v_degree1;
		double* m_polynomial1;
		double* m_du_polynomial1;
		double* m_dv_polynomial1;
	};

	typedef Solver<TwoPolynomialEquation, IntervalVector<2>,
		IntervalVector<2>, IntervalVector<2>, IntervalMatrix<2, 2>> TwoPolynomialEquationSolver;


	inline TwoPolynomialEquation::TwoPolynomialEquation(
		int u_degree0, int v_degree0, double* polynomial0, double* du_polynomial0, double* dv_polynomial0,
		int u_degree1, int v_degree1, double* polynomial1, double* du_polynomial1, double* dv_polynomial1) :
		m_u_degree0(u_degree0),
		m_v_degree0(v_degree0),
		m_polynomial0(polynomial0),
		m_du_polynomial0(du_polynomial0),
		m_dv_polynomial0(dv_polynomial0),
		m_u_degree1(u_degree1),
		m_v_degree1(v_degree1),
		m_polynomial1(polynomial1),
		m_du_polynomial1(du_polynomial1),
		m_dv_polynomial1(dv_polynomial1) {
	}

	inline int TwoPolynomialEquation::GetEquationCount() {
		return 2;
	}

	inline int TwoPolynomialEquation::GetVariableCount() {
		return 2;
	}

	inline double TwoPolynomialEquation::GetVariableEpsilon(int i) {
		return 1E-6;
	}

	inline double TwoPolynomialEquation::GetValueEpsilon(int i, bool is_checking) {
		return is_checking ? 1E-6 : 1E-12;
	}

	inline void TwoPolynomialEquation::CalculateValue(const IntervalVector<2>& variable, IntervalVector<2>& value) {
		Interval u = variable.Get(0);
		Interval v = variable.Get(1);
		value.Set(0, estimate_two_polynomial_interval(m_u_degree0, m_v_degree0, m_polynomial0, u, v));
		value.Set(1, estimate_two_polynomial_interval(m_u_degree1, m_v_degree1, m_polynomial1, u, v));
	}

	inline void TwoPolynomialEquation::CalculatePartialDerivative(const IntervalVector<2>& variable, IntervalMatrix<2, 2>& value) {
		Interval u = variable.Get(0);
		Interval v = variable.Get(1);
		*value.Get(0, 0) = estimate_two_polynomial_interval(m_u_degree0 - 1, m_v_degree0, m_du_polynomial0, u, v);
		*value.Get(0, 1) = estimate_two_polynomial_interval(m_u_degree0, m_v_degree0, m_dv_polynomial0, u, v);
		*value.Get(1, 0) = estimate_two_polynomial_interval(m_u_degree1 - 1, m_v_degree1, m_du_polynomial1, u, v);
		*value.Get(1, 1) = estimate_two_polynomial_interval(m_u_degree1, m_v_degree1, m_dv_polynomial1, u, v);
	}

	inline int TwoPolynomialEquation::GetSplitIndex(const IntervalVector<2>& variable, int prev_split_index, double size) {
		return prev_split_index == 0 ? 1 : 0;
	}

	inline int TwoPolynomialEquation::CompareIteratePriority(const IntervalVector<2>& variable1, double size1,
		const IntervalVector<2>& variable2, double size2) {
		if (size1 < size2) {
			return -1;
		}
		if (size1 > size2) {
			return 1;
		}
		return 0;
	}

	inline bool TwoPolynomialEquation::PreIterate(IntervalVector<2>* variable, SolverIteratedResult& result, double& size) {
		return false;
	}

	inline bool TwoPolynomialEquation::CheckFinished(const Array<SolverHeapItem<IntervalVector<2>>>& heap) {
		return heap.GetCount() >= (m_u_degree0 > m_u_degree1 ? m_u_degree0 : m_u_degree1) * (m_v_degree0 > m_v_degree1 ? m_v_degree0 : m_v_degree1);
	}

}

#endif

