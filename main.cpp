#include "curve.h"

#include <iostream>

void Testing() {
    Param prm;
    // вывожу наборы параметров
    std::cout << "----------------------------------------------------------"
                 "--------------------------------------------------------\n";
    std::cout << "Параметры из стандарта id-tc26-gost-3410-2012-256-ParamSetA:\n";
    std::cout << "p = " << prm.p << '\n'
              << "a = " << prm.a << '\n'
              << "x_base = " << prm.x_base << '\n'
              << "y_base = " << prm.y_base << '\n'
              << "q = " << prm.q << "\n\n";
    std::cout << "Предрасчитанный параметр theta = " << prm.theta << "\n";
    std::cout << "-----------------------------------------------------------"
                 "-------------------------------------------------------\n\n";

    // параметры квадрики
    std::cout << "-----------------------------------------------------------"
                 "-------------------------------------------------------\n";
    std::cout << "Параметры Квадрики Якоби:\n";
    JacobiCurve curve(prm);
    std::cout << "e = " << curve.e << '\n'
              << "d = " << curve.d << '\n'
              << "X_base" << curve.X << '\n'
              << "Y_base" << curve.Y << '\n'
              << "Z_base" << curve.Z << "\n";
    std::cout << "-----------------------------------------------------------"
                 "-------------------------------------------------------\n\n";


    // тестирую работу с нейтральным элементом
    std::cout << "-------------------------------------------------------------"
                 "-----------------------------------------------------\n";
    std::cout << "ТЕСТ 1:\n";
    Point E(0, 1, 1);
    std::cout << "Нейтральный элемент Е:\n"
              << " – в проективных координатах:\n";
    PrintProjective(E);
    std::cout << " - в афинных:\n";
    PrintInAffin(E, curve);

    if (CheckPointIsOnCurve(E, curve))
        std::cout << "Точка E находится на кривой\n";
    else
        std::cout << "Точка E не находится на кривой\n";
    std::cout << "--------------------------------------------------------------"
                 "----------------------------------------------------\n\n";

    // тестирую функций с порождающим элементом
    std::cout << "-----------------------------------------------------------"
                 "-------------------------------------------------------\n";
    std::cout << "ТЕСТ 2:\n";
    Point P_base;
    P_base.X = curve.X;
    P_base.Y = curve.Y;
    P_base.Z = curve.Z;

    std::cout << "Порождающий элемент P_base в аффинных координатах:\n";
    PrintInAffin(P_base, curve);
    if (CheckPointIsOnCurve(P_base, curve))
        std::cout << "Точка P_base находится на кривой\n";
    else
        std::cout << "Точка P_base не находится на кривой\n";
    std::cout << "--------------------------------------------------------"
                 "----------------------------------------------------------\n\n";

    // тест E+P=P
    std::cout << "---------------------------------------------------------"
                 "---------------------------------------------------------\n";
    std::cout << "ТECT 3:\n";
    std::cout << "проверим, что E+P_base = P_base\n";
    Point P1(2, 2, 2);
    AddPoints(E, P_base, P1, curve);
    if (!IsEqual(P_base, P1, curve))
        std::cout << "E+P и P равны\n";
    else
        std::cout << "E+P и P не равны\n";
    std::cout << "-------------------------------------------------------"
                 "-----------------------------------------------------------\n\n";

    // тест выбрать любую точку и проверить лежит ли она кривой
    std::cout << "-------------------------------------------------------"
                 "-----------------------------------------------------------\n";
    std::cout << "ТECT 4:\n";
    std::cout << "проверим пренадлежит ли точка P2=(3:6:8) кривой\n";
    std::cout << "P3 в аффинных:\n";
    Point P2(3, 6, 8);
    PrintInAffin(P2, curve);
    if (CheckPointIsOnCurve(P2, curve))
        std::cout << "точка P2 находится на кривой\n";
    else
        std::cout << "точка P2 не находится на кривой\n";
    std::cout << "----------------------------------------------------"
                 "--------------------------------------------------------------\n\n";

    // тест qP = E
    std::cout << "------------------------------------------------------"
                 "------------------------------------------------------------\n";
    std::cout << "ТECT 5:\n";
    std::cout << "проверим, что qP = E\n";
    Point resPoint(1, 1, 1);
    CalculateDegree(resPoint, P_base, curve, prm.q);
    std::cout << "нейтральный элемент в аффинных координатах:\n";
    PrintInAffin(E, curve);
    std::cout << "qP в аффинных координатах:\n";
    PrintInAffin(resPoint, curve);
    std::cout << "-------------------------------------------------------"
                 "-----------------------------------------------------------\n\n";

    // тест [q+1]P = P и [q-1] = -P
    std::cout << "----------------------------------------------------------"
                 "--------------------------------------------------------\n";
    std::cout << "ТECT 6:\n";
    std::cout << "проверим, что [q+1]P = P и [q-1] = -P\n";
    mpz_class degree = prm.q + mpz_class(1);
    CalculateDegree(resPoint, P_base, curve, degree);
    std::cout << "[q+1]P:\n";
    PrintInAffin(resPoint, curve);
    std::cout << "P:\n";
    PrintInAffin(P_base, curve);
    if(!(IsEqual(resPoint, P_base, curve)))
        std::cout << "[q+1]P равно P\n";
    else
        std::cout << "[q+1]P не равно P\n";

    degree = prm.q - mpz_class(1);
    CalculateDegree(resPoint, P_base, curve, degree);
    std::cout << "[q-1]P:\n";
    PrintInAffin(resPoint, curve);
    std::cout << "-P:\n";
    Point negP;
    GetNegativePoint(negP, P_base);
    PrintInAffin(negP, curve);
    if (!IsEqual(resPoint, negP, curve))
        std::cout << "[q-1]P равно -P\n";
    else
        std::cout << "[q+1]P не равно P\n";
    std::cout << "------------------------------------------------------------"
                 "------------------------------------------------------\n\n";

    // придумать любое k , такое что 0 <= k <q  и посчитать [k]P и проверить пренадлежность кривой
    std::cout << "------------------------------------------------------------"
                 "------------------------------------------------------\n";
    std::cout << "ТECT 7:\n";
    std::cout << "вычислим кратную точку [k]P. пусть k = 100, P = P_base\n";
    degree = 100;
    CalculateDegree(resPoint, P_base, curve, degree);
    PrintInAffin(resPoint, curve);
    if (CheckPointIsOnCurve(resPoint, curve))
        std::cout << "точка [k]P находится на кривой\n";
    else
        std::cout << "точка [k]P не находится на кривой\n";
    std::cout << "-----------------------------------------------------------"
                 "-------------------------------------------------------\n\n";

    //сгенерировать случайное k , такое что 0 <=k < q и проверить пренадлежность кривой
    std::cout << "------------------------------------------------------------"
                 "------------------------------------------------------\n";
    std::cout << "Тест 8:\n";
    std::cout << "сгенерирую случайное k в диапозоне 0 <= k < q\n";
    mpz_class k;
    gmp_randstate_t rnd_state;
    gmp_randinit_mt(rnd_state);
    mpz_urandomm(k.get_mpz_t(), rnd_state, prm.q.get_mpz_t());
    std::cout << "k: " << k << '\n';
    std::cout << "посчитаю [k]P:\n";
    CalculateDegree(resPoint, P_base, curve, k);
    PrintInAffin(resPoint, curve);
    if (CheckPointIsOnCurve(resPoint, curve))
        std::cout << "точка [k]P находится на кривой\n";
    else
        std::cout << "точка [k]P не находится на кривой\n";
    std::cout << "---------------------------------------------------------"
                 "---------------------------------------------------------\n\n";

    // проверить [k1]P + [k2]P = [k1 + k2]P
    std::cout << "----------------------------------------------------------"
                 "--------------------------------------------------------\n";
    std::cout << "Тест 9:\n";
    std::cout << "[k1]P + [k2]P = [k1 + k2]P\n";
    mpz_class k1;
    mpz_class k2;
    std::cout << "сгенерирую k1 и k2\n";
    mpz_class maxrand("100000000000000000");
    mpz_urandomm(k1.get_mpz_t(), rnd_state, maxrand.get_mpz_t());
    mpz_urandomm(k2.get_mpz_t(), rnd_state, maxrand.get_mpz_t());
    std::cout << "k1 = " << k1 << '\n'
              << "k2 = " << k2 << '\n';
    k = k1 + k2;
    std::cout << "k = k1 + k2 = " << k << '\n';
    Point res1(0, 1, 1);
    Point res2(0, 1, 1);
    Point res3(0, 1, 1);
    CalculateDegree(res1, P_base, curve, k1);
    CalculateDegree(res2, P_base, curve, k2);
    CalculateDegree(res3, P_base, curve, k);
    AddPoints(res1, res2, resPoint, curve);

    if (IsEqual(resPoint, res3, curve))
        std::cout << "[k1]P + [k2]P равно [k1 + k2]P\n";
    else
        std::cout << "[k1]P + [k2]P не равно [k1 + k2]P\n";

    if (CheckPointIsOnCurve(res3, curve))
        std::cout << "точка [k]P находится на кривой\n";
    else
        std::cout << "точка [k]P не находится на кривой\n";

    gmp_randclear(rnd_state);
}

int main() {
    Testing();
    return 0;
}
