#include "curve.h"

#include <iostream>
#include <limits>
#include <utility>

Param::Param()
        : p(p_str),
        a(a_str),
        x_base(x_base_str),
        y_base(y_base_str),
        q(q_str),
        theta(theta_str)
        {}

JacobiCurve::JacobiCurve(const Param &param) {
    this->p = param.p;

    this->e = (mpz_class(3) * param.theta * param.theta) % this->p;
    this->e += mpz_class(4) * param.a;
    mpz_neg(this->e.get_mpz_t(), this->e.get_mpz_t());
    mpz_class invert_16 = 16;
    mpz_invert(invert_16.get_mpz_t(), invert_16.get_mpz_t(), this->p.get_mpz_t());
    this->e *= invert_16;
    this->e %= this->p;
    if (mpz_sgn(this->e.get_mpz_t()) == -1)
        this->e += this->p;

    this->d = mpz_class(3) * param.theta;
    mpz_class invert_4 = 4;
    mpz_invert(invert_4.get_mpz_t(), invert_4.get_mpz_t(), this->p.get_mpz_t());
    this->d *= invert_4;
    this->d %= this->p;

    mpz_class x_base_minux_theta = param.x_base - param.theta;
    this->X = (mpz_class(2) * x_base_minux_theta) % this->p;
    if (mpz_sgn(this->X.get_mpz_t()) == -1)
        this->X += this->p;

    this->Y = (2 * param.x_base + param.theta) % this->p;
    this->Y *= x_base_minux_theta;
    this->Y %= this->p;
    this->Y *= x_base_minux_theta;
    this->Y %= this->p;
    mpz_class y_base_sqr = (param.y_base * param.y_base) % this->p;
    this->Y -= y_base_sqr;
    this->Y %= this->p;
    if (mpz_sgn(this->Y.get_mpz_t()) == -1)
        this->Y += this->p;

    this->Z = param.y_base;
}

JacobiPoint::JacobiPoint(const std::string &x, const std::string &y, const std::string &z)
        : X(x), Y(y), Z(z) {}

JacobiPoint::JacobiPoint(int x, int y, int z)
        : X(x), Y(y), Z(z) {}

JacobiPoint::JacobiPoint(mpz_class x, mpz_class y, mpz_class z)
        : X(std::move(x)), Y(std::move(y)), Z(std::move(z)) {}

void AddPoints(const JacobiPoint &P1, const JacobiPoint &P2, JacobiPoint &P_res, const JacobiCurve &curve) {
    mpz_class T1 = P1.X;
    mpz_class T2 = P1.Y;
    mpz_class T3 = P1.Z;
    mpz_class T4 = P2.X;
    mpz_class T5 = P2.Y;
    mpz_class T6 = P2.Z;
    mpz_class T7;
    mpz_class T8;

    T7 = (T1 * T3) % curve.p;
    T7 = (T7 + T2) % curve.p;
    T8 = (T4 * T6) % curve.p;
    T8 = (T8 + T5) % curve.p;
    T2 = (T2 * T5) % curve.p;
    T7 = (T7 * T8) % curve.p;
    T7 = (T7 - T2) % curve.p;
    T5 = (T1 * T4) % curve.p;
    T1 = (T1 + T3) % curve.p;
    T8 = (T3 * T6) % curve.p;
    T4 = (T4 + T6) % curve.p;
    T6 = (T5 * T8) % curve.p;
    T7 = (T7 - T6) % curve.p;
    T1 = (T1 * T4) % curve.p;
    T1 = (T1 - T5) % curve.p;
    T1 = (T1 - T8) % curve.p;
    T3 = (T1 * T1) % curve.p;
    T6 = (T6 + T6) % curve.p;
    T3 = (T3 - T6) % curve.p;
    T4 = (curve.e * T6) % curve.p;
    T3 = (T3 * T4) % curve.p;
    T4 = (curve.d * T6) % curve.p;
    T2 = (T2 - T4) % curve.p;
    T4 = (T8 * T8) % curve.p;
    T8 = (T5 * T5) % curve.p;
    T8 = (curve.e * T8) % curve.p;
    T5 = (T4 + T8) % curve.p;
    T2 = (T2 * T5) % curve.p;
    T2 = (T2 + T3) % curve.p;
    T5 = (T4 - T8) % curve.p;

    T7 %= curve.p;
    T2 %= curve.p;
    T5 %= curve.p;
    if (mpz_sgn(T7.get_mpz_t()) == -1)
        T7 += curve.p;
    if (mpz_sgn(T2.get_mpz_t()) == -1)
        T2 += curve.p;
    if (mpz_sgn(T5.get_mpz_t()) == -1)
        T5 += curve.p;

    P_res.X = T7;
    P_res.Y = T2;
    P_res.Z = T5;
}

void kPowPoint(JacobiPoint &kP, const JacobiPoint &P, const JacobiCurve &curve, const mpz_class &degree) {
    mp_bitcnt_t bit_count;
    for (mp_bitcnt_t i = 0; i != std::numeric_limits<mp_bitcnt_t>::max();
            i = mpz_scan1(degree.get_mpz_t(), i + 1)) {
        bit_count = i;
    }
    ++bit_count;

    JacobiPoint R = P;
    JacobiPoint Q(0, 1, 1);

    for (int i = bit_count; i > 0; --i) {
        if (mpz_tstbit(degree.get_mpz_t(), i - 1)) {
            AddPoints(Q, R, Q, curve);
            AddPoints(R, R, R, curve);
        } else {
            AddPoints(R, Q, R, curve);
            AddPoints(Q, Q, Q, curve);
        }
    }

    kP = Q;
}

void AffineCast(JacobiPoint &affine_repr, const JacobiPoint &P, const JacobiCurve &curve) {
    mpz_class x;
    mpz_class y;
    mpz_class z;

    mpz_class z_inverted;

    mpz_invert(z_inverted.get_mpz_t(), P.Z.get_mpz_t(), curve.p.get_mpz_t());
    x = (z_inverted * P.X) % curve.p;
    if (mpz_sgn(x.get_mpz_t()) == -1)
        x += curve.p;

    y = (z_inverted * z_inverted) % curve.p;
    y = (y * P.Y) % curve.p;
    if (mpz_sgn(y.get_mpz_t()) == -1)
        y += curve.p;

    z = 0;

    affine_repr = {std::move(x), std::move(y), std::move(z)};
}

void AffineRepr(const JacobiPoint &point, const JacobiCurve &curve) {
    JacobiPoint affine_point(0, 1, 1);

    AffineCast(affine_point, point, curve);
    std::cout << "x = " << affine_point.X << "\ny = " << affine_point.Y << "\n\n";
}

void ProjectiveRepr(const JacobiPoint &P) {
    std::cout << "X = " << P.X << "\nY = " << P.Y << "\nZ = " << P.Z << "\n\n";
}

int CheckPoint(const JacobiPoint &P, const JacobiCurve &curve) {
    mpz_class left;
    mpz_class right;
    mpz_class buf1;
    mpz_class buf2 = 4;

    left = (P.Y * P.Y) % curve.p;
    mpz_powm(right.get_mpz_t(), P.X.get_mpz_t(), buf2.get_mpz_t(), curve.p.get_mpz_t());
    right = (right * curve.e) % curve.p;
    mpz_powm(buf2.get_mpz_t(), P.Z.get_mpz_t(), buf2.get_mpz_t(), curve.p.get_mpz_t());
    right += buf2;
    buf2 = (P.X * P.Z) % curve.p;
    buf2 = (buf2 * buf2) % curve.p;
    buf2 = (curve.d * buf2) % curve.p;
    buf2 += buf2;
    right -= buf2;
    buf1 = (left - right) % curve.p;

    int ans = (buf1 == 0);

    return ans;
}

int CheckEqualPoints(const JacobiPoint &P1, const JacobiPoint &P2, const JacobiCurve &curve) {
    JacobiPoint affineP1(0, 1, 1);
    JacobiPoint affineP2(0, 1, 1);
    AffineCast(affineP1, P1, curve);
    AffineCast(affineP2, P2, curve);

    int ans;
    if (affineP1.X == affineP2.X && affineP1.Y == affineP2.Y)
        ans = 1;
    else
        ans = 0;

    return ans;
}

void GetNegativePoint(JacobiPoint &res, const JacobiPoint &point) {
    res = point;
    mpz_neg(res.X.get_mpz_t(), res.X.get_mpz_t());
}

