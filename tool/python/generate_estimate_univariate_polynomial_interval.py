import sympy
import numpy


def format_num(d):
    sd = str(d)
    if sd.find('.') != -1:
        sd = sd.rstrip('0').rstrip('.')
    return sd


def generate_matrix_code1(max_degree, base_tab_count):
    s = '\n'
    for degree in range(2, max_degree + 1):
        s += 'const double bezier_standard_fitting_matrix_{}d[{}][{}] = {{\n'.format(degree, degree + 1, degree + 1)
        u = sympy.symbols('u')
        exprs = []
        for i in range(degree + 1):
            exprs.append(sympy.binomial(degree, i) * (u ** i) * ((1 - u) ** (degree - i)))
        matrix = []
        for i in range(degree + 1):
            row = []
            for j in range(degree + 1):
                row.append(exprs[j].subs(u, sympy.Rational(i, degree)))
            matrix.append(row)
        matrix_inv = numpy.linalg.inv(numpy.array(sympy.Matrix(matrix), dtype=numpy.float64))
        s2 = ''
        for i in range(degree + 1):
            s2 += '\t{ '
            s2 += format_num(matrix_inv[i][0])
            for j in range(1, degree + 1):
                s2 += ', '
                s2 += format_num(matrix_inv[i][j])
            s2 += ' },\n'
        s2 = s2.rstrip('\n').rstrip(',')
        s2 += '\n};\n\n'
        s += s2

    base_tab = ''
    for i in range(base_tab_count):
        base_tab += '    '
    return s.replace('\n', '\n' + base_tab)


def generate_normal_func_code1(max_degree, base_tab_count):
    s = '''
typedef Interval (*estimate_univariate_polynomial_interval_nd_func)(double* polynomial, const Interval& t);

Interval estimate_univariate_polynomial_interval_0d(double* polynomial, const Interval& t) {
    return polynomial[0];
}

Interval estimate_univariate_polynomial_interval_1d(double* polynomial, const Interval& t) {
    double s[2] = {
        calculate_univariate_polynomial_value(1, polynomial, t.Min),
        calculate_univariate_polynomial_value(1, polynomial, t.Max)
    };
    Interval result;
    if (s[0] < s[1]) {
        result = Interval(s[0], s[1]);
    }
    else {
        result = Interval(s[1], s[0]);
    }
    return result;
}

'''
    for degree in range(2, max_degree + 1):
        s2 = ''
        s2 += 'Interval estimate_univariate_polynomial_interval_{}d(double* polynomial, const Interval& t) {{\n'.format(
            degree)
        s2 += '    double tlen = t.Length();\n'
        s2 += '    double s[{}] = {{\n'.format(degree + 1)
        s2 += '        calculate_univariate_polynomial_value({}, polynomial, t.Min),\n'.format(degree)
        for i in range(1, degree):
            s2 += '        calculate_univariate_polynomial_value({}, polynomial, {}.0 / {} * tlen + t.Min),\n'.format(
                degree, i, degree)
        s2 += '        calculate_univariate_polynomial_value({}, polynomial, t.Max)\n'.format(degree)
        s2 += '    };\n'
        s2 += '    Interval result;\n'
        s2 += '    if (s[0] < s[{}]) {{\n'.format(degree)
        s2 += '        result = Interval(s[0], s[{}]);\n'.format(degree)
        s2 += '    }\n'
        s2 += '    else {\n'
        s2 += '        result = Interval(s[{}], s[0]);\n'.format(degree)
        s2 += '    }\n'
        s2 += '    for (int i = 1; i < {}; ++i) {{\n'.format(degree)
        s2 += '        double p = 0;\n'
        s2 += '        for (int j = 0; j <= {}; ++j) {{\n'.format(degree)
        s2 += '            p += bezier_standard_fitting_matrix_{}d[i][j] * s[j];\n'.format(degree)
        s2 += '        }\n'
        s2 += '        result.Merge(p);\n'
        s2 += '    }\n'
        s2 += '    return result;\n'
        s2 += '}\n\n'
        s += s2
    s2 = ''
    s2 += 'const estimate_univariate_polynomial_interval_nd_func estimate_univariate_polynomial_interval_nd_funcs[{}] = {{\n'.format(max_degree + 1)
    for degree in range(0, max_degree):
        s2 += '    estimate_univariate_polynomial_interval_{}d,\n'.format(degree)
    s2 += '    estimate_univariate_polynomial_interval_{}d\n'.format(max_degree)
    s2 += '};\n\n'
    s += s2
    s2 = ''
    s2 += '''Interval estimate_univariate_polynomial_interval(int degree, double* polynomial, const Interval& t) {{
    if (degree > {}) {{
        throw "degree is too large";
    }}
    return estimate_univariate_polynomial_interval_nd_funcs[degree](polynomial, t);   
}}

'''.format(max_degree)
    s += s2

    base_tab = ''
    for i in range(base_tab_count):
        base_tab += '    '
    return s.replace('\n', '\n' + base_tab)


def generate_rational_func_code1(max_degree, base_tab_count):
    s = '''
typedef Interval (*estimate_univariate_rational_polynomial_interval_nd_func)(double* n_polynomial, double* d_polynomial, const Interval& t);

Interval estimate_univariate_rational_polynomial_interval_0d(double* n_polynomial, double* d_polynomial, const Interval& t) {
    return Interval(n_polynomial[0]) / Interval(d_polynomial[0]);
}

Interval estimate_univariate_rational_polynomial_interval_1d(double* n_polynomial, double* d_polynomial, const Interval& t) {
    double n[2] = {
        calculate_univariate_polynomial_value(1, n_polynomial, t.Min),
        calculate_univariate_polynomial_value(1, n_polynomial, t.Max)
    };
    double d[2] = {
        calculate_univariate_polynomial_value(1, d_polynomial, t.Min),
        calculate_univariate_polynomial_value(1, d_polynomial, t.Max)
    };
    double a = n[0] / d[0];
    double b = n[1] / d[1];
    Interval result;
    if (a < b) {
        result = Interval(a, b);
    }
    else {
        result = Interval(b, a);
    }
    return result;
}

'''
    for degree in range(2, max_degree + 1):
        s2 = ''
        s2 += 'Interval estimate_univariate_rational_polynomial_interval_{}d(double* n_polynomial, double* d_polynomial, const Interval& t) {{\n'.format(
            degree)
        s2 += '    double tlen = t.Length();\n'
        s2 += '    double n[{}] = {{\n'.format(degree + 1)
        s2 += '        calculate_univariate_polynomial_value({}, n_polynomial, t.Min),\n'.format(degree)
        for i in range(1, degree):
            s2 += '        calculate_univariate_polynomial_value({}, n_polynomial, {}.0 / {} * tlen + t.Min),\n'.format(
                degree, i, degree)
        s2 += '        calculate_univariate_polynomial_value({}, n_polynomial, t.Max)\n'.format(degree)
        s2 += '    };\n'
        s2 += '    double d[{}] = {{\n'.format(degree + 1)
        s2 += '        calculate_univariate_polynomial_value({}, d_polynomial, t.Min),\n'.format(degree)
        for i in range(1, degree):
            s2 += '        calculate_univariate_polynomial_value({}, d_polynomial, {}.0 / {} * tlen + t.Min),\n'.format(
                degree, i, degree)
        s2 += '        calculate_univariate_polynomial_value({}, d_polynomial, t.Max)\n'.format(degree)
        s2 += '    };\n'
        s2 += '    double a = n[0] / d[0];\n'
        s2 += '    double b = n[{}] / d[{}];\n'.format(degree, degree)
        s2 += '    Interval result;\n'
        s2 += '    if (a < b) {\n'
        s2 += '        result = Interval(a, b);\n'
        s2 += '    }\n'
        s2 += '    else {\n'
        s2 += '        result = Interval(b, a);\n'
        s2 += '    }\n'
        s2 += '    for (int i = 1; i < {}; ++i) {{\n'.format(degree)
        s2 += '        double wp = 0;\n'
        s2 += '        double w = 0;\n'
        s2 += '        for (int j = 0; j <= {}; ++j) {{\n'.format(degree)
        s2 += '            wp += bezier_standard_fitting_matrix_{}d[i][j] * n[j];\n'.format(degree)
        s2 += '            w += bezier_standard_fitting_matrix_{}d[i][j] * d[j];\n'.format(degree)
        s2 += '        }\n'
        s2 += '        result.Merge(wp / w);\n'
        s2 += '    }\n'
        s2 += '    return result;\n'
        s2 += '}\n\n'
        s += s2
    s2 = ''
    s2 += 'const estimate_univariate_rational_polynomial_interval_nd_func estimate_univariate_rational_polynomial_interval_nd_funcs[{}] = {{\n'.format(max_degree + 1)
    for degree in range(0, max_degree):
        s2 += '    estimate_univariate_rational_polynomial_interval_{}d,\n'.format(degree)
    s2 += '    estimate_univariate_rational_polynomial_interval_{}d\n'.format(max_degree)
    s2 += '};\n\n'
    s += s2
    s2 = ''
    s2 += '''Interval estimate_univariate_rational_polynomial_interval(int degree, double* n_polynomial, double* d_polynomial, const Interval& t) {{
    if (degree > {}) {{
        throw "degree is too large";
    }}
    return estimate_univariate_rational_polynomial_interval_nd_funcs[degree](n_polynomial, d_polynomial, t);
}}

'''.format(max_degree)
    s += s2

    base_tab = ''
    for i in range(base_tab_count):
        base_tab += '    '
    return s.replace('\n', '\n' + base_tab)


def generate_matrix_code2(max_u_degree, max_v_degree, base_tab_count):
    s = '\n'
    for u_degree in range(1, max_u_degree + 1):
        for v_degree in range(1, max_v_degree + 1):
            s += 'const double bezier_standard_fitting_matrix_{}d_{}d[{}][{}] = {{\n'.format(u_degree, v_degree, (u_degree + 1) * (v_degree + 1), (u_degree + 1) * (v_degree + 1))
            u = sympy.symbols('u')
            v = sympy.symbols('v')
            exprs = []
            for i in range(u_degree + 1):
                for j in range(v_degree + 1):
                    exprs.append((sympy.binomial(u_degree, i) * (u ** i) * ((1 - u) ** (u_degree - i))) * (sympy.binomial(v_degree, j) * (v ** j) * ((1 - v) ** (v_degree - j))))
            matrix = []
            for i in range(u_degree + 1):
                for j in range(v_degree + 1):
                    row = []
                    for k in range((u_degree + 1) * (v_degree + 1)):
                        row.append(exprs[k].subs([(u, sympy.Rational(i, u_degree)), (v, sympy.Rational(j, v_degree))]))
                    matrix.append(row)
            matrix_inv = numpy.linalg.inv(numpy.array(sympy.Matrix(matrix), dtype=numpy.float64))
            s2 = ''
            for i in range((u_degree + 1) * (v_degree + 1)):
                s2 += '\t{ '
                s2 += format_num(matrix_inv[i][0])
                for j in range(1, (u_degree + 1) * (v_degree + 1)):
                    s2 += ', '
                    s2 += format_num(matrix_inv[i][j])
                s2 += ' },\n'
            s2 = s2.rstrip('\n').rstrip(',')
            s2 += '\n};\n\n'
            s += s2

    base_tab = ''
    for i in range(base_tab_count):
        base_tab += '    '
    return s.replace('\n', '\n' + base_tab)


def generate_cs(max_degree, base_tab_count):
    s = ''
    for degree in range(0, max_degree + 1):
        s2 = 'const double g_c{}[{}] = {{\n'.format(degree, degree + 1)
        for i in range(0, degree):
            s2 += '    {},\n'.format(sympy.N(sympy.binomial(degree, i), 19))
        s2 += '    {}\n'.format(sympy.N(sympy.binomial(degree, degree), 19))
        s2 += '};\n\n'
        s += s2
    s2 = 'const double* g_c[{}] = {{\n'.format(max_degree + 1)
    for degree in range(0, max_degree):
        s2 += '    g_c{},\n'.format(degree)
    s2 += '    g_c{}\n'.format(max_degree)
    s2 += '};\n\n'
    s += s2
    base_tab = ''
    for i in range(base_tab_count):
        base_tab += '    '
    return s.replace('\n', '\n' + base_tab)


print('''/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include "wstd/polynomial.h"

namespace wgp {
''')
# print(generate_matrix_code1(30, 1))
# print(generate_normal_func_code1(30, 1))
# print(generate_rational_func_code1(30, 1))

# print(generate_matrix_code2(5, 5, 1))

print(generate_cs(30, 1))

print('''
}
''')
