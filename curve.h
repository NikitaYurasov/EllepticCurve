#ifndef NULL
#define NULL (void*)0
#endif

#ifndef CURVE_H
#define CURVE_H

#include <gmpxx.h>
#include <string>

// id-tc26-gost-3410-2012-256-ParamSetA:
#define p_str "115792089237316195423570985008687907853269984665640564039457584007913129639319"

#define a_str "87789765485885808793369751294406841171614589925193456909855962166505018127157"
#define x_base_str "65987350182584560790308640619586834712105545126269759365406768962453298326056"
#define y_base_str "22855189202984962870421402504110399293152235382908105741749987405721320435292"
#define q_str "28948022309329048855892746252171976963338560298092253442512153408785530358887"

#define theta_str "454069018412434321972378083527459607666454479745512801572100703902391945898"

/**
 * Структура хранения параметров стандарта
 */
struct Param {
    Param();

    mpz_class p;
    mpz_class a;
    mpz_class x_base;
    mpz_class y_base;
    mpz_class q;
    mpz_class theta;
};

/**
 * Структура хранения параметров эллиптической кривой в форме квадратики Якоби
 */
struct JacobiCurve {
    JacobiCurve(const Param &param);

    mpz_class Y = 0;
    mpz_class X = 0;
    mpz_class e = 0;
    mpz_class d = 0;
    mpz_class Z = 0;
    mpz_class p = 0;
};

/**
 * Структура для хранения параметров (координат) точки
 */
struct JacobiPoint {
    JacobiPoint(const std::string &x, const std::string &y, const std::string &z);

    JacobiPoint(int x, int y, int z);

    JacobiPoint(mpz_class x, mpz_class y, mpz_class z);

    JacobiPoint() = default;

    mpz_class X;
    mpz_class Y;
    mpz_class Z;
};

/**
 * Складывает две точки P1 и P2. Результат заносится в переменную (точку) P_res
 * @param P1: ссылка на точку №1 для сложения
 * @param P2: ссылка на точку №2 для сложения
 * @param P_res: ссылка на точку, в которую будет записан результат
 * @param curve: структура кривой типа JacobiCurve, в которой хранятся параметры текущей кривой
 */
void AddPoints(const JacobiPoint &P1, const JacobiPoint &P2, JacobiPoint &P_res, const JacobiCurve &curve);

/**
 * Возведение точки Р в степень degree. Используется алгоритм <<Лесенка Монтгомери>>
 * @param kP: ссылка на точку, в которую будет записан результат
 * @param P: ссылка на точку, которая будет возводиться в степень
 * @param curve: структура, хранящая текущие параметры кривой
 * @param degree: значение степени
 */
void kPowPoint(JacobiPoint &kP, const JacobiPoint &P, const JacobiCurve &curve, const mpz_class &degree);

/**
 * Переводит точку из проективных координат в аффинные
 * @param affine_repr: ссылка на точку, в которую будет записан результат в аффинных координатах
 * @param P: ссылка на точку в проективных координатах
 * @param curve: структура, хранящая текущие параметры кривой
 */
void AffineCast(JacobiPoint &affine_repr, const JacobiPoint &P, const JacobiCurve &curve);

/**
 * Выводит на экран координаты точки в аффинном представлении
 * @param point: ссылка на точку
 * @param curve: структура, хранящая текущие параметры кривой
 */
void AffineRepr(const JacobiPoint &point, const JacobiCurve &curve);

/**
 * Выводит на экран координаты точки в проективном представлении
 * @param P: ссылка на точку
 */
void ProjectiveRepr(const JacobiPoint &P);

/**
 * Проверяет, лежит ли точка на кривой.
 * Возвращает 1, если точка лежит на кривой, 0 -- в противном случае
 * @param P: ссылка на точку
 * @param curve: структура, хранящая текущие параметры кривой
 * @return int
 */
int CheckPoint(const JacobiPoint &P, const JacobiCurve &curve);

/**
 * Проверяет равны ли точки друг другу.
 * Возвращает 1, если равны; 0 -- в противном случае
 * @param P1 : ссылка на точку №1
 * @param P2 : ссылка на точку №2
 * @param curve : структура, хранящая текущие параметры кривой
 * @return int
 */
int CheckEqualPoints(const JacobiPoint &P1, const JacobiPoint &P2, const JacobiCurve &curve);

/**
 * Записывает в res -point
 * @param res : ссылка на точку, в которую будет записан результат
 * @param point : ссыка на точку, которую необходимо преставить в отрицательном виде
 */
void GetNegativePoint(JacobiPoint &res, const JacobiPoint &point);

#endif // CURVE_H
