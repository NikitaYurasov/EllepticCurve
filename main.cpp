#include "curve.h"

#include <iostream>

void Testing() {
    Param prm;
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

    std::cout << "-----------------------------------------------------------"
                 "-------------------------------------------------------\n";
    std::cout << "Параметры Квадрики Якоби:\n";
    JacobiCurve curve(prm);
    std::cout << "e = " << curve.e << '\n'
              << "d = " << curve.d << '\n'
              << "X_base = " << curve.X << '\n'
              << "Y_base = " << curve.Y << '\n'
              << "Z_base = " << curve.Z << "\n";
    std::cout << "-----------------------------------------------------------"
                 "-------------------------------------------------------\n\n";


    std::cout << "-------------------------------------------------------------"
                 "-----------------------------------------------------\n";
    std::cout << "ТЕСТ: ПРОВЕРКА ПРИНАДЛЕЖНОСТИ НЕЙТРАЛЬНОГО ЭЛЕМЕНТА\n";
    JacobiPoint E(0, 1, 1);
    std::cout << "Нейтральный элемент Е:\n"
              << " – в проективных координатах:\n";
    ProjectiveRepr(E);
    std::cout << " - в аффинных координатах:\n";
    AffineRepr(E, curve);
    std::cout << "Ответ: ";
    if (CheckPoint(E, curve))
        std::cout << "Точка E находится на кривой\n";
    else
        std::cout << "Точка E не находится на кривой\n";
    std::cout << "--------------------------------------------------------------"
                 "----------------------------------------------------\n\n";

    std::cout << "-----------------------------------------------------------"
                 "-------------------------------------------------------\n";
    std::cout << "ТЕСТ 2:\n";
    JacobiPoint P_base;
    P_base.X = curve.X;
    P_base.Y = curve.Y;
    P_base.Z = curve.Z;

    std::cout << "Порождающий элемент в аффинных координатах:\n";
    AffineRepr(P_base, curve);

    std::cout << "Ответ: ";
    if (CheckPoint(P_base, curve))
        std::cout << "Точка находится на кривой\n";
    else
        std::cout << "Точка не находится на кривой\n";
    std::cout << "--------------------------------------------------------"
                 "----------------------------------------------------------\n\n";

    std::cout << "---------------------------------------------------------"
                 "---------------------------------------------------------\n";
    std::cout << "ТECT 3: ";
    std::cout << "E+P_base = P_base?\n";
    JacobiPoint P1(2, 2, 2);
    AddPoints(E, P_base, P1, curve);

    std::cout << "Ответ: ";
    if (CheckEqualPoints(P_base, P1, curve))
        std::cout << "E+P == P\n";
    else
        std::cout << "E+P != P\n";
    std::cout << "-------------------------------------------------------"
                 "-----------------------------------------------------------\n\n";

    std::cout << "-------------------------------------------------------"
                 "-----------------------------------------------------------\n";
    std::cout << "ТECT 4:\n";
    std::cout << "Принадлежит ли точка P2=(5:1:4) кривой\n";
    std::cout << "P2 в аффинных:\n";
    JacobiPoint P2(5, 1, 4);
    AffineRepr(P2, curve);
    std::cout << "Ответ: ";
    if (CheckPoint(P2, curve))
        std::cout << "точка P2 находится на кривой\n";
    else
        std::cout << "точка P2 не находится на кривой\n";
    std::cout << "----------------------------------------------------"
                 "--------------------------------------------------------------\n\n";

    std::cout << "------------------------------------------------------"
                 "------------------------------------------------------------\n";
    std::cout << "ТECT 5: ";
    std::cout << "qP = E?\n";
    JacobiPoint resPoint(1, 1, 1);
    kPowPoint(resPoint, P_base, curve, prm.q);
    std::cout << "нейтральный элемент в аффинных координатах:\n";
    AffineRepr(E, curve);
    std::cout << "qP в аффинных координатах:\n";
    AffineRepr(resPoint, curve);
    std::cout << "-------------------------------------------------------"
                 "-----------------------------------------------------------\n\n";

    std::cout << "----------------------------------------------------------"
                 "--------------------------------------------------------\n";
    std::cout << "ТECT 6: ";
    std::cout << "[q+1]P = P и [q-1] = -P\n";
    mpz_class degree = prm.q + mpz_class(1);
    kPowPoint(resPoint, P_base, curve, degree);
    std::cout << "[q+1]P:\n";
    AffineRepr(resPoint, curve);
    std::cout << "P:\n";
    AffineRepr(P_base, curve);

    std::cout << "Ответ: ";
    if((CheckEqualPoints(resPoint, P_base, curve)))
        std::cout << "[q+1]P == P\n";
    else
        std::cout << "[q+1]P != P\n";

    degree = prm.q - mpz_class(1);
    kPowPoint(resPoint, P_base, curve, degree);
    std::cout << "[q-1]P:\n";
    AffineRepr(resPoint, curve);
    std::cout << "-P:\n";
    JacobiPoint negP;
    GetNegativePoint(negP, P_base);
    AffineRepr(negP, curve);

    std::cout << "Ответ: ";
    if (CheckEqualPoints(resPoint, negP, curve))
        std::cout << "[q-1]P == -P\n";
    else
        std::cout << "[q-1]P != -P\n";
    std::cout << "------------------------------------------------------------"
                 "------------------------------------------------------\n\n";

    std::cout << "------------------------------------------------------------"
                 "------------------------------------------------------\n";
    std::cout << "ТECT 7: ";
    std::cout << "Вычисление [k]P при k = 100; P = P_base?\n";
    degree = 100;
    kPowPoint(resPoint, P_base, curve, degree);
    AffineRepr(resPoint, curve);

    std::cout << "Ответ: ";
    if (CheckPoint(resPoint, curve))
        std::cout << "точка [k]P находится на кривой\n";
    else
        std::cout << "точка [k]P не находится на кривой\n";
    std::cout << "-----------------------------------------------------------"
                 "-------------------------------------------------------\n\n";

    std::cout << "------------------------------------------------------------"
                 "------------------------------------------------------\n";
    std::cout << "Тест 8:\n";
    std::cout << "Случайное k в диапазоне 0 <= k < q\n";
    mpz_class k;
    gmp_randstate_t rnd_state;
    gmp_randinit_mt(rnd_state);
    mpz_urandomm(k.get_mpz_t(), rnd_state, prm.q.get_mpz_t());
    std::cout << "k: " << k << '\n';
    std::cout << "[k]P в аффинных координатах:\n";
    kPowPoint(resPoint, P_base, curve, k);
    AffineRepr(resPoint, curve);

    std::cout << "Ответ: ";
    if (CheckPoint(resPoint, curve))
        std::cout << "точка [k]P находится на кривой\n";
    else
        std::cout << "точка [k]P не находится на кривой\n";
    std::cout << "---------------------------------------------------------"
                 "---------------------------------------------------------\n\n";

    std::cout << "----------------------------------------------------------"
                 "--------------------------------------------------------\n";
    std::cout << "Тест 9: ";
    std::cout << "[k1]P + [k2]P = [k1 + k2]P?\n";
    mpz_class k1;
    mpz_class k2;
    std::cout << "Случайные k1 и k2\n";
    mpz_class maxrand("100000000000000000");
    mpz_urandomm(k1.get_mpz_t(), rnd_state, maxrand.get_mpz_t());
    mpz_urandomm(k2.get_mpz_t(), rnd_state, maxrand.get_mpz_t());
    std::cout << "k1 = " << k1 << '\n'
              << "k2 = " << k2 << '\n';
    k = k1 + k2;
    std::cout << "k = k1 + k2 = " << k << '\n';
    JacobiPoint res1(0, 1, 1);
    JacobiPoint res2(0, 1, 1);
    JacobiPoint res3(0, 1, 1);
    kPowPoint(res1, P_base, curve, k1);
    kPowPoint(res2, P_base, curve, k2);
    kPowPoint(res3, P_base, curve, k);
    AddPoints(res1, res2, resPoint, curve);

    std::cout << "Ответ: ";
    if (CheckEqualPoints(resPoint, res3, curve))
        std::cout << "[k1]P + [k2]P == [k1 + k2]P\n";
    else
        std::cout << "[k1]P + [k2]P != [k1 + k2]P\n";

    std::cout << "Ответ: ";
    if (CheckPoint(res3, curve))
        std::cout << "точка [k]P находится на кривой\n";
    else
        std::cout << "точка [k]P не находится на кривой\n";

    gmp_randclear(rnd_state);
}

int main() {
    Testing();
    return 0;
}
