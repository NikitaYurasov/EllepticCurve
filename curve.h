#ifndef NULL
#define NULL (void*)0
#endif

#ifndef CURVE_H
#define CURVE_H

#include <gmpxx.h>
#include <string>

// Набор параметров id-tc26-gost-3410-2012-256-ParamSetA:
#define p_str               "115792089237316195423570985008687907853269984665640564039457584007913129639319"

#define a_str               "87789765485885808793369751294406841171614589925193456909855962166505018127157"
#define x_base_str          "65987350182584560790308640619586834712105545126269759365406768962453298326056"
#define y_base_str          "22855189202984962870421402504110399293152235382908105741749987405721320435292"
#define q_str               "28948022309329048855892746252171976963338560298092253442512153408785530358887"

// Значение, вычисленное на основе вышеперечисленных параметров с помощью Wolfram Mathematica
#define theta_str           "454069018412434321972378083527459607666454479745512801572100703902391945898"

// структура для хранения параметров в виде больших чисел
struct Param {
    Param();

    mpz_class p;
    mpz_class a;
    mpz_class x_base;
    mpz_class y_base;
    mpz_class q;
    mpz_class theta;
};

// структура для хранения параметров квадрики и порождающих элементов.
struct JacobiCurve {
    JacobiCurve(const Param &param);

    mpz_class Y = 0;
    mpz_class X = 0;
    mpz_class e = 0;
    mpz_class d = 0;
    mpz_class Z = 0;
    mpz_class p = 0;
};

// структура для хранения точек
struct Point {
    Point(const std::string &x, const std::string &y, const std::string &z);

    Point(int x, int y, int z);

    Point(mpz_class x, mpz_class y, mpz_class z);

    Point() = default;

    mpz_class X;
    mpz_class Y;
    mpz_class Z;
};

// реализация сложения двух точек
void AddPoints(const Point &P1, const Point &P2, Point &P3, const JacobiCurve &curve);

// вычисление кратной точки с помощью алгоритма "лесенка Монтгомери"
void CalculateDegree(Point &kP, const Point &P, const JacobiCurve &curve, const mpz_class &degree);

// перевод в аффинные координаты
void CastPointToAffine(Point &affpoint, const Point &P, const JacobiCurve &curve);

// вывод в аффинных координатах
void PrintInAffin(const Point &point, const JacobiCurve &curve);

// вывод в проективных
void PrintProjective(const Point &point);

// проверка, что точка лежит на кривой
int CheckPointIsOnCurve(const Point &P, const JacobiCurve &curve);

// проверка, что точки равны
int IsEqual(const Point &P1, const Point &P2, const JacobiCurve &curve);

//получение обратной точки
void GetNegativePoint(Point &res, const Point &point);

#endif // CURVE_H
