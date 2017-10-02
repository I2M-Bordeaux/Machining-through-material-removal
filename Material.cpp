//     Copyright (C) 2017 : Renault Chloe, Cahuc Olivier, Darnis Philippe, Raynald Laheurte, Delos Vincent
//
//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU Lesser General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU Lesser General Public License for more details.
//
//     You should have received a copy of the GNU Lesser General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
/// \file Material.cpp
/// \author Renault Chloe (chloe.renault@u-bordeaux.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include "Material.h"
#include <iostream>
#include <cmath>

using namespace std;

Material::Material(double A, double B, double C, double n, double m, double Kth, double a) :\
_A(A), _B(B), _C(C), _n(n), _m(m), _Kth(Kth), _a(a)
{
}

Material::~Material()
{
}


double Material::JohnsonCookModelTemp(double eps_eq, double d_eps_eq, double d_eps_0, double T_moy, double T0, double T_fus)
{
  _sigm_eq = (_A+_B*(pow(eps_eq,_n)))*(1+_C*log((d_eps_eq)/(d_eps_0)))*(1-pow((T_moy-T0)/(T_fus-T0),_m));
}


double Material::JohnsonCookModel(double eps_eq, double d_eps_eq, double d_eps_0)
{
  _sigm_eq = (_A+_B*pow(eps_eq,_n))*(1+_C*log((d_eps_eq)/(d_eps_0)));
}
