/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_INTERVAL_
#define _WGP_STD_INTERVAL_

#include "utils.h"

namespace wgp {

    struct WGP_API Interval {
    public:
        double Min;
        double Max;
    public:
        Interval();
        Interval(double d);
        Interval(double min, double max);
        double Center() const;
        double Length() const;
        bool IsIntersected(double d, double epsilon) const;
        bool IsIntersected(const Interval& other, double epsilon) const;
        bool IsInner(const Interval& other) const;
        void Extend(double d);
        void Merge(const Interval& other);
        void Intersect(const Interval& other);
        Interval Secure() const;
    };

    inline Interval operator+(const Interval& x, const Interval& y) {
        return Interval(x.Min + y.Min, x.Max + y.Max).Secure();
    }

    inline Interval operator+(const Interval& x, double y) {
        return Interval(x.Min + y, x.Max + y).Secure();
    }

    inline Interval operator+(double x, const Interval& y) {
        return Interval(y.Min + x, y.Max + x).Secure();
    }

    inline Interval operator-(const Interval& x, const Interval& y) {
        return Interval(x.Min - y.Max, x.Max - y.Min).Secure();
    }

    inline Interval operator-(const Interval& x, double y) {
        return Interval(x.Min - y, x.Max - y).Secure();
    }

    inline Interval operator-(double x, const Interval& y) {
        return Interval(x - y.Max, x - y.Min).Secure();
    }

    inline Interval operator-(const Interval& x) {
        return Interval(-x.Max, -x.Min).Secure();
    }

    inline Interval operator*(const Interval& x, const Interval& y) {
        double min = x.Min * y.Min;
        double max = min;
        double t = x.Max * y.Min;
        if (t < min) {
            min = t;
        } else if (t > max) {
            max = t;
        }
        t = x.Max * y.Max;
        if (t < min) {
            min = t;
        } else if (t > max) {
            max = t;
        }
        t = x.Min * y.Max;
        if (t < min) {
            min = t;
        } else if (t > max) {
            max = t;
        }
        return Interval(min, max).Secure();
    }

    inline Interval operator*(const Interval& x, double y) {
        if (y >= 0) {
            return Interval(x.Min * y, x.Max * y).Secure();
        }
        return Interval(x.Max * y, x.Min * y).Secure();
    }

    inline Interval operator*(double x, const Interval& y) {
        if (x >= 0) {
            return Interval(y.Min * x, y.Max * x).Secure();
        }
        return Interval(y.Max * x, y.Min * x).Secure();
    }

    inline Interval operator/(const Interval& x, const Interval& y) {
        if (y.Max < 0 || y.Min > 0) {
            double min = x.Min / y.Min;
            double max = min;
            double t = x.Max / y.Min;
            if (t < min) {
                min = t;
            } else if (t > max) {
                max = t;
            }
            t = x.Max / y.Max;
            if (t < min) {
                min = t;
            } else if (t > max) {
                max = t;
            }
            t = x.Min / y.Max;
            if (t < min) {
                min = t;
            } else if (t > max) {
                max = t;
            }
            return Interval(min, max).Secure();
        }
        throw "div zero error";
    }

    inline Interval operator/(const Interval& x, double y) {
        if (y > 0) {
            return Interval(x.Min / y, x.Max / y).Secure();
        } 
        if (y < 0) {
            return Interval(x.Max / y, x.Min / y).Secure();
        }
        throw "div zero error";
    }

    inline Interval operator/(double x, const Interval& y) {
        if (y.Max < 0 || y.Min > 0) {
            if (x >= 0) {
                return Interval(x / y.Max, x / y.Min).Secure();
            }
            return Interval(x / y.Min, x / y.Max).Secure();
        }
        throw "div zero error";
    }

    inline Interval::Interval() {
        Min = 0;
        Max = 0;
    }

    inline Interval::Interval(double d) {
        Min = d;
        Max = d;
    }

    inline Interval::Interval(double min, double max) {
        Min = min;
        Max = max;
    }

    inline double Interval::Center() const {
        return (Max + Min) * 0.5;
    }

    inline double Interval::Length() const {
        return Max - Min;
    }

    inline bool Interval::IsIntersected(double d, double epsilon) const {
        return d <= Max + epsilon && d >= Min - epsilon;
    }

    inline bool Interval::IsIntersected(const Interval& other, double epsilon) const {
        return other.Min <= Max + epsilon && other.Max >= Min - epsilon;
    }

    inline bool Interval::IsInner(const Interval& other) const {
        return other.Min <= Min && other.Max >= Max;
    }

    inline void Interval::Extend(double d) {
        Min -= d;
        Max += d;
    }

    inline void Interval::Merge(const Interval& other) {
        if (other.Min < Min) {
            Min = other.Min;
        }
        if (other.Max > Max) {
            Max = other.Max;
        }
    }

    inline void Interval::Intersect(const Interval& other) {
        if (other.Min > Min) {
            Min = other.Min;
        }
        if (other.Max < Max) {
            Max = other.Max;
        }
    }

    inline Interval Interval::Secure() const {
        if (Min * Max <= 0) {
            return Interval(Min + Min * 1.0099995520574283e-15, Max + Max * 1.0099995520574283e-15);
        }
        if (Min < 0) {
            return Interval(Min + Min * 1.0099995520574283e-15, Max);
        }
        return Interval(Min, Max + Max * 1.0099995520574283e-15);
    }

}

inline wgp::Interval sqr(const wgp::Interval& x) {
    double d0 = x.Min * x.Min;
    double d1 = x.Max * x.Max;
    if (x.Max < 0) {
        return wgp::Interval(d1, d0).Secure();
    }
    if (x.Min > 0) {
        return wgp::Interval(d0, d1).Secure();
    }
    if (d0 < d1) {
        return wgp::Interval(0, d1).Secure();
    }
    return wgp::Interval(0, d0).Secure();
}

inline wgp::Interval sqrt(const wgp::Interval& x) {
    if (x.Min >= 0) {
        return wgp::Interval(sqrt(x.Min), sqrt(x.Max)).Secure();
    }
    throw "sqrt error";
}

inline wgp::Interval cos(const wgp::Interval& x) {
    double d = x.Length();
    if (d >= g_pi * 2) {
        return wgp::Interval(-1, 1);
    }
    double min = x.Min - int(x.Min / (g_pi * 2)) * (g_pi * 2);
    if (min < 0) {
        min += g_pi * 2;
    }
    double max = min + d;
    if (min <= g_pi) {
        if (max <= g_pi) {
            return wgp::Interval(cos(max), cos(min)).Secure();
        }
        if (max <= g_pi * 2) {
            double c0 = cos(min);
            double c1 = cos(max);
            return wgp::Interval(-1, c0 > c1 ? c0 : c1).Secure();
        }
        return wgp::Interval(-1, 1);
    } else {
        if (max < g_pi * 2) {
            return wgp::Interval(cos(min), cos(max)).Secure();
        }
        if (max < g_pi * 3) {
            double c0 = cos(min);
            double c1 = cos(max);
            return wgp::Interval(c0 < c1 ? c0 : c1, 1).Secure();
        }
        return wgp::Interval(-1, 1);
    }
}

inline wgp::Interval sin(const wgp::Interval& x) {
    return cos(x - g_pi * 0.5);
}

inline void sincos(const wgp::Interval& x, wgp::Interval* sinx, wgp::Interval* cosx) {
    //todo 优化
    *sinx = sin(x);
    *cosx = cos(x);
}

inline wgp::Interval pow(const wgp::Interval& x, int y) {
    if (y & 1) {
        return wgp::Interval(pow(x.Min, y), pow(x.Max, y)).Secure();
    }
    double a = pow(x.Min, y);
    double b = pow(x.Max, y);
    if (x.Min * x.Max < 0) {
        return a < b ? wgp::Interval(0, b).Secure() : wgp::Interval(0, a).Secure();
    }
    return a < b ? wgp::Interval(a, b).Secure() : wgp::Interval(b, a).Secure();
}

inline wgp::Interval pow(const wgp::Interval& x, const wgp::Interval& y) {
    if (x.Min < 0) {
        if (y.Max != y.Min) {
            throw "domain error";
        }
        if (y.Min != (int)y.Min) {
            throw "domain error";
        }
        double a = pow(x.Min, y.Min);
        double b = pow(x.Max, y.Min);
        if ((int)y.Min & 1) {
            return wgp::Interval(a, b).Secure();
        }
        if (x.Min * x.Max < 0) {
            return a < b ? wgp::Interval(0, b).Secure() : wgp::Interval(0, a).Secure();
        }
        return a < b ? wgp::Interval(a, b).Secure() : wgp::Interval(b, a).Secure();
    }
    else {
        double a = pow(x.Min, y.Min);
        double b = a;
        double c = pow(x.Min, y.Max);
        if (c < a) {
            a = c;
        }
        if (c > b) {
            b = c;
        }
        c = pow(x.Max, y.Min);
        if (c < a) {
            a = c;
        }
        if (c > b) {
            b = c;
        }
        c = pow(x.Max, y.Max);
        if (c < a) {
            a = c;
        }
        if (c > b) {
            b = c;
        }
        if (x.Min * x.Max <= 0) {
            if (0 < a) {
                a = 0;
            }
            if (0 > b) {
                b = 0;
            }
        }
        if (y.Min * y.Max <= 0) {
            if (1 < a) {
                a = 1;
            }
            if (1 > b) {
                b = 1;
            }
        }
        return wgp::Interval(a, b).Secure();
    }
}

inline wgp::Interval log(const wgp::Interval& x) {
    if (x.Min < 0) {
        throw "domain error";
    }
    return wgp::Interval(log(x.Min), log(x.Max)).Secure();
}
    
#endif