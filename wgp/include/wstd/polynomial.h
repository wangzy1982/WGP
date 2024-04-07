/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_POLYNOMIAL_
#define _WGP_STD_POLYNOMIAL_

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
}

#endif

