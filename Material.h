#include <iostream>
#include <string>


class Material
{
  public:

    Material(double A, double B, double C, double n, double m, double Kth, double a);
    ~Material();

    double JohnsonCookModel(double eps_eq, double d_eps_eq, double d_eps_0);   ///Application de la loi de comportement
    double JohnsonCookModelTemp(double eps_eq, double d_eps_eq, double d_eps_0, double T_moy, double T0, double T_fus);
    double GetEquivalentStress(){return _sigm_eq;};
    double GetThermalPower(){return _Kth;};
    double GetDiffusivity(){return _a;};

  private:

    const double _A;  ///Limite d'élasticité
    const double _B;   ///Coeff lié à l'écrouissage
    const double _C;   ///Coeff de sensibilité à la vitesse de déformation
    const double _n;  ///Coeff lié à l'écrouissage
    const double _m;   ///Coeff de sensibilité à la température
    double _T_moy;   /// Température moyenne
    double _T_fus;    ///Température de fusion
    double _T0;      ///Température initiale du matériau
    double _sigm_eq; /// Contrainte équivalente

    double _Kth, _a;

};
