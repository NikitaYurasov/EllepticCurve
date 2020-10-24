#include "curve.h"

#include <iostream>
#include <limits>
#include <utility>

Param::Param()
    : p(p_str)
    , a(a_str)
    , x_base(x_base_str)
    , y_base(y_base_str)
    , q(q_str)
    , theta(theta_str)
{}

JacobiCurve::JacobiCurve(const Param& param) {
    this->p = param.p;

    //e = -(3 * theta^2 + 4 * a) / 16
    this->e = (mpz_class(3) * param.theta * param.theta) % this->p;                 // e = 3 * theta^2
    this->e += mpz_class(4) * param.a;                                              // e = 3 * theta^2 + 4 * a
    mpz_neg(this->e.get_mpz_t(), this->e.get_mpz_t());                              // e = -(3 * theta^2 + 4 * a)
    mpz_class invert_16 = 16;
    mpz_invert(invert_16.get_mpz_t(), invert_16.get_mpz_t(), this->p.get_mpz_t());  // invert_16 = 1 / 16
    this->e *= invert_16;                                                           // e = -(3 * theta^2 + 4 * a) / 16
    this->e %= this->p;
    if (mpz_sgn(this->e.get_mpz_t()) == -1)
        this->e += this->p;

    // d = 3 * theta / 4
    this->d = mpz_class(3) * param.theta;                                           // d = 3 * theta
    mpz_class invert_4 = 4;
    mpz_invert(invert_4.get_mpz_t(), invert_4.get_mpz_t(), this->p.get_mpz_t());    // invert_4 = 1 / 4
    this->d *= invert_4;                                                            // d = 3 * theta / 4
    this->d %= this->p;

    // X = 2 * (x_base - theta)
    mpz_class x_base_minux_theta = param.x_base - param.theta;
    this->X = (mpz_class(2) * x_base_minux_theta) % this->p;
    if (mpz_sgn(this->X.get_mpz_t()) == -1)
        this->X += this->p;

    // Y = (2 * x_base + theta) * (x_base - theta)^2 - y_base^2
    this->Y = (2 * param.x_base + param.theta) % this->p;                           // Y = 2 * x_base + theta
    this->Y *= x_base_minux_theta;                                                  // Y = (2 * x_base + theta) * (x_base - theta)
    this->Y %= this->p;
    this->Y *= x_base_minux_theta;                                                  // Y = (2 * x_base + theta) * (x_base - theta)^2
    this->Y %= this->p;
    mpz_class y_base_sqr = (param.y_base * param.y_base) % this->p;                 // y_base_sqr = y_base^2
    this->Y -= y_base_sqr;                                                          // (2 * x_base + theta) * (x_base - theta)^2 - y_base^2
    this->Y %= this->p;
    if (mpz_sgn(this->Y.get_mpz_t()) == -1)
        this->Y += this->p;

    // Z = y_base
    this->Z = param.y_base;
}

Point::Point(const std::string& x, const std::string& y, const std::string& z)
    : X(x)
    , Y(y)
    , Z(z)
{}

Point::Point(int x, int y, int z)
    : X(x)
    , Y(y)
    , Z(z)
{}

Point::Point(mpz_class x, mpz_class y, mpz_class z)
    : X(std::move(x))
    , Y(std::move(y))
    , Z(std::move(z))
{}

void AddPoints(const Point& P1, const Point& P2, Point& P3, const JacobiCurve& curve) {
    mpz_class T1 = P1.X;
    mpz_class T2 = P1.Y;
    mpz_class T3 = P1.Z;
    mpz_class T4 = P2.X;
    mpz_class T5 = P2.Y;
    mpz_class T6 = P2.Z;
    mpz_class T7;
    mpz_class T8;

    // алгоритм сложения

    T7 = (T1 * T3) % curve.p;       // T7 = X1 * Z1
    T7 += T2;                       // T7 = X1 * Z1 + Y1
    T8 = (T4 * T6) % curve.p;       // T8 = X2 * Z2
    T8 += T5;                       // T8 = X2 * Z2 + Y2
    T2 = (T2 * T5) % curve.p;       // T2 = Y1 * Y2
    T7 = (T7 * T8) % curve.p;       // T7 = X3 + Y1 * Y2 + X1 * X2 * Z1 * Z2
    T7 -= T2;                       // T7 = X3 + X1 * X2 * Z1 * Z2
    T5 = (T1 * T4) % curve.p;       // T5 = X1 * X2
    T1 += T3;                       // T1 = X1 + Z1
    T8 = (T3 * T6) % curve.p;       // T8 = Z1 * Z2
    T4 += T6;                       // T4 = X2 + Z2
    T6 = (T5 * T8) % curve.p;       // T6 = X1 * X2 * Z1 * Z2
    T7 -= T6;                       // T7 = X3
    T1 = (T1 * T4) % curve.p;       // T1 = X1 * Z2 + X2 * Z1 + X1 * X2 + Z1 * Z2
    T1 -= T5;                       // T1 = X1 * Z2 + X2 * Z1 + Z1 * Z2
    T1 -= T8;                       // T1 = X1 * Z2 + X2 * Z1
    T3 = (T1 * T1) % curve.p;       // T3 = X1^2 * Z2^2+ X2^2 * Z1^2 + 2 * X1 * X2 * Z1 * Z2
    T6 += T6;                       // T6 = 2 * X1 * X2 * Z1 * Z2
    T3 -= T6;                       // T3 = X1^2 * Z2^2+ X2^2 * Z1^2
    T4 = (curve.e * T6) % curve.p;  // T4 = 2 * e * X1 * X2 * Z1 * Z2
    T3 = (T3 * T4) % curve.p;       // T3 = 2 * e * X1 * X2 * Z1 * Z2 * (X1^2 * Z2^2+ X2^2 * Z1^2)
    T4 = (curve.d * T6) % curve.p;  // T4 = 2 * d * X1 * X2 * Z1 * Z2
    T2 -= T4;                       // T2 = Y1 * Y2 - 2 * d * X1 * X2 * Z1 * Z2
    T4 = (T8 * T8) % curve.p;       // T4 = Z1^2 * Z2^2
    T8 = (T5 * T5) % curve.p;       // T8 = X1^2 * X2^2
    T8 = (curve.e * T8) % curve.p;  // T8 = e * X1^2 * X2^2
    T5 = T4 + T8;                   // T5 = Z1^2 * Z2^2 + e * X1^2 * X2^2
    T2 = (T2 * T5) % curve.p;       // T2 = (Z1^2 * Z2^2 + e * X1^2 * X2^2) * (Y1 * Y2 - 2 * d * X1 * X2 * Z1 * Z2)
    T2 += T3;                       // T2 = Y3
    T5 = T4 - T8;                   // T5 = Z3

    T7 %= curve.p;
    T2 %= curve.p;
    T5 %= curve.p;
    if (mpz_sgn(T7.get_mpz_t()) == -1)
        T7 += curve.p;
    if (mpz_sgn(T2.get_mpz_t()) == -1)
        T2 += curve.p;
    if (mpz_sgn(T5.get_mpz_t()) == -1)
        T5 += curve.p;

    P3.X = T7;                      // X3 = T7
    P3.Y = T2;                      // Y3 = T2
    P3.Z = T5;                      // Z3 = T5
}

void CalculateDegree(Point& kP, const Point& P, const JacobiCurve& curve, const mpz_class& degree) {
    mp_bitcnt_t bitnums = 0;
    for (mp_bitcnt_t i = 0; i != std::numeric_limits<mp_bitcnt_t>::max(); i = mpz_scan1(degree.get_mpz_t(), i + 1)) {
        bitnums = i;
    }
    ++bitnums;

    Point R = P;
    Point Q(0, 1, 1);

    for (int i = bitnums; i > 0; --i) {
        if (mpz_tstbit(degree.get_mpz_t(), i - 1)) {
            AddPoints(Q, R, Q, curve); // Q = R + Q
            AddPoints(R, R, R, curve); // R = R + R
        } else {
            AddPoints(R, Q, R, curve); // R = R + Q
            AddPoints(Q, Q, Q, curve); // Q = Q + Q
        }
    }

    kP = Q;
}

void CastPointToAffine(Point& affpoint, const Point& P, const JacobiCurve& curve) {
    mpz_class x;
    mpz_class y;
    mpz_class z;

    mpz_class z_inverted;

    // x = X / Z
    mpz_invert(z_inverted.get_mpz_t(), P.Z.get_mpz_t(), curve.p.get_mpz_t());   // x = 1 / Z
    x = (z_inverted * P.X) % curve.p;                                           // x = X / Z
    if (mpz_sgn(x.get_mpz_t()) == -1)
        x += curve.p;

    // y = Y / Z^2
    y = (z_inverted * z_inverted) % curve.p;                                    // y = 1 / Z^2
    y = (y * P.Y) % curve.p;                                                    // y = Y / Z^2
    if (mpz_sgn(y.get_mpz_t()) == -1)
        y += curve.p;

    // z = 0
    z = 0;

    affpoint = {std::move(x), std::move(y), std::move(z)};
}

void PrintInAffin(const Point& point, const JacobiCurve& curve) {
    Point affinepoint(0, 1, 1);

    CastPointToAffine(affinepoint, point, curve);
    std::cout << "x = " << affinepoint.X << "\ny = " << affinepoint.Y << "\n\n";
}

void PrintProjective(const Point& P) {
    std::cout << "X = " << P.X << "\nY = " << P.Y << "\nZ = " << P.Z << "\n\n";
}

int CheckPointIsOnCurve(const Point& P, const JacobiCurve& curve) {
    mpz_class left;
    mpz_class right;
    mpz_class buf1;
    mpz_class buf2 = 4;

    /* Y^2 = e * X^4 - 2 * d * X^2 * Z^2 + Z^4 - подставить точки,
     * вычислить значения левой правой стороны
     * если левая и правая сторона равны точка лежит на прямой
     */

    left = (P.Y * P.Y) % curve.p;                                                           // left = Y^2
    mpz_powm(right.get_mpz_t(), P.X.get_mpz_t(), buf2.get_mpz_t(), curve.p.get_mpz_t());    // right = X^4
    right = (right * curve.e) % curve.p;                                                    // right = e * X^4
    mpz_powm(buf2.get_mpz_t(), P.Z.get_mpz_t(), buf2.get_mpz_t(), curve.p.get_mpz_t());     // buf2 = Z^4
    right += buf2;                                                                          // right = e * X^4 + Z^4
    buf2 = (P.X * P.Z) % curve.p;                                                           // buf2 = X * Z
    buf2 = (buf2 * buf2) % curve.p;                                                         // buf2 = X^2 * Z^2
    buf2 = (curve.d * buf2) % curve.p;                                                      // buf2 = d * X^2 * Z^2
    buf2 += buf2;                                                                           // buf2 = 2 * d * X^2 * Z^2
    right -= buf2;                                                                          // right = e * X^4 - 2 * d * X^2 * Z^2 + Z^4
    buf1 = (left - right) % curve.p;                                                        // buf1 = Y^2 - (e * X^4 - 2 * d * X^2 * Z^2 + Z^4)

    int ans = (buf1 == 0);

    return ans;
}
int IsEqual(const Point& P1, const Point& P2, const JacobiCurve& curve) {
    Point affineP1(0, 1, 1);
    Point affineP2(0, 1, 1);
    CastPointToAffine(affineP1, P1, curve);
    CastPointToAffine(affineP2, P2, curve);

    int ans = -100;
    if (affineP1.X == affineP2.X && affineP1.Y == affineP2.Y)
        ans = 0; // точки равны
    else
        ans = -1; // не равны

    return ans;
}
void GetNegativePoint(Point& res, const Point& point)
{
    res = point;
    mpz_neg(res.X.get_mpz_t(), res.X.get_mpz_t());
}

