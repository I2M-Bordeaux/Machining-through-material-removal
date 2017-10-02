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
/// \file Machining.cpp
/// \author Renault Chloe (chloe.renault@u-bordeaux.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include "Machining.h"
#include <iostream>
#include <cmath>
#include <vector>


#include <boost/math/special_functions/bessel.hpp>

using namespace std;

// Constructeur par défaut
Machining::Machining(Material& material, Tool& tool, double phi) :
_material(material), _tool(tool), _phi(phi)
{
  _pi = atan(1.)*4;

  _V_c = _tool.GetV_c();
  _f = _tool.GetFeed();
  _gamma_0 = _tool.GetGamma();
  _z = _tool.GetTeeth();
  _omega = _tool.GetAngularVelocity();
  _R_0 = _tool.GetToolEdgeRadius();

  _rbec = 0.04*pow(10,-3);

  _strain_rate_0 = 0.001;
}

// Destructeur par défaut
Machining::~Machining()
{
}



double Machining::Minimize()
{
  double out(0), Q(pow(10,10));
  int n(51);

  double R_min, k;
  R_min = (((2*pow(sin(_pi/4- _gamma_0/2),2)-1)*(n+1)/(n-1)+2)*_f/2*cos(_pi/4 - _gamma_0/2)/cos(_gamma_0))/2;
  k = 1./2 * pow(R_min*cos(_gamma_0), 2) / pow(1 + sin(_gamma_0), 3);

  while (out == 0)
  { _phi += 0.01;

    vector<double> s;
    vector<double> x_phi;
    vector<double> x_phi_pos;
    for (int i = 0 ; i < 1280 ; ++i)
    {
      s.push_back (0.0003+0.0001*i);
      x_phi.push_back (s[i]*cos(_phi-_gamma_0) - (k*sin(_phi)/s[i]));
      x_phi_pos.push_back (s[i]*cos(_phi-_gamma_0) - (k*sin(_phi)/s[i]));
      if(x_phi_pos[i] < 0)x_phi_pos[i] *= -1;
    }
    vector<double>::iterator min =  min_element(x_phi_pos.begin(),x_phi_pos.end());
    double j(distance(x_phi_pos.begin(), min));

    _vect_phi.resize(2);
    _vect_phi[0] = x_phi[j];
    _vect_phi[1] = s[j]*sin(_phi-_gamma_0) + _rbec/sin(_phi) + k*sin(_phi)/s[j];

    double g, f1, f2;

    g = 4*k*cos(_phi-_gamma_0)*sin(_phi);

    f1 = cos(_gamma_0)*(_vect_phi[0]+sqrt(pow(_vect_phi[0],2)+g));
    f2 = 4*cos(_phi-_gamma_0)*sin(_phi)*sqrt(pow(_vect_phi[0],2)+g);
    _eps_eq_zp = sqrt(4/3)* f1/f2;

    f1 = 4*_V_c*k*sin(_phi)*cos(_gamma_0)*pow(_vect_phi[0]+sqrt(pow(_vect_phi[0],2)+g),2);
    f2 = sqrt(pow(_vect_phi[0],2)+g) * pow(g + pow(_vect_phi[0]+sqrt(pow(_vect_phi[0],2)+g),2),2);
    _d_eps_eq_zp = sqrt(4/3)*f1/f2;

    _material.JohnsonCookModel(_eps_eq_zp, _d_eps_eq_zp, _strain_rate_0);
    _sigm_eq_zp = _material.GetEquivalentStress();

    _Q_thq_zp = _sigm_eq_zp * _d_eps_eq_zp;

    if (_Q_thq_zp > Q) {
      out = 1;
    }

    Q = _Q_thq_zp;

  }
  return _phi;
}



OrthogonalMachining::OrthogonalMachining(Material& material, Tool& tool, double phi) :
Machining(material, tool, phi)
{
}


void OrthogonalMachining::Initialize()
{
  _V_cop = _V_c*sin(_phi)/cos(_gamma_0-_phi);   ///Merchant
  _l = _f / (_z*sin(_phi));

  _V_s = _V_c*cos(_gamma_0)/cos(_gamma_0-_phi);

}



void OrthogonalMachining::ResultsZCP()
{
  int n(50);
  double R_min, k;
  R_min = ((2*pow(sin(_pi/4- _gamma_0),2)-1/(n-1))*(n+1)+1)*_f/_l*cos(_pi/4- _gamma_0/2)/cos(_gamma_0)/2;
  k = 1./2 * R_min*pow(cos(_gamma_0), 2) / pow(1 + sin(_gamma_0), 3);

  vector<double> s;
  vector<double> x_phi;
  for (int i = 0 ; i < 1279 ; ++i)
  {
    s.push_back(0.0003+0.0001*i);
    x_phi.push_back(s[i]*cos(_phi-_gamma_0) - (k*sin(_phi)/s[i]));
  }
  vector<double>::iterator min =  min_element(x_phi.begin(),x_phi.end());
  double j(distance(x_phi.begin(), min));

  _vect_phi.resize(2);
  _vect_phi[0] = x_phi[j];
  _vect_phi[1] = s[j]*sin(_phi-_gamma_0) + _rbec/sin(_phi) + k*sin(_phi)/s[j];
  double g, f1, f2;
  g = 4*k*cos(_phi-_gamma_0)*sin(_phi);

  f1 = cos(_gamma_0)*(_vect_phi[0]+sqrt(pow(_vect_phi[0],2)+g));
  f2 = 4*cos(_phi-_gamma_0)*sin(_phi)*sqrt(pow(_vect_phi[0],2)+g);
  _eps_eq_zp = sqrt(4/3)* f1/f2;

  f1 = 4*_V_c*k*sin(_phi)*cos(_gamma_0)*pow(_vect_phi[0]+sqrt(pow(_vect_phi[0],2)+g),2);
  f2 = sqrt(pow(_vect_phi[0],2)+g) * pow(g + pow(_vect_phi[0]+sqrt(pow(_vect_phi[0],2)+g),2),2);
  _d_eps_eq_zp = sqrt(4/3)*f1/f2;

  _strain_rate_0 = 0.001;

  _material.JohnsonCookModel(_eps_eq_zp, _d_eps_eq_zp, _strain_rate_0);
  _sigm_eq_zp = _material.GetEquivalentStress();

  _Q_thq_zp = _sigm_eq_zp * _d_eps_eq_zp;


  double K_th = _material.GetThermalPower();
  double a(_material.GetDiffusivity());
  a = a * pow(10, 6);

  double Q(0);

  for (int i = 0 ; i < n ; ++i)
    {
      _vect_s.resize(2);
      _vect_s[0] = -_l*cos(_phi)*(i-1+0.5)/n;
      _vect_s[1] = _l*sin(_phi)*(i-1+0.5)/n;

      _vect_si.resize(2);
      _vect_si[0] = _vect_s[0];
      _vect_si[1] = 2*_f - _vect_s[1];

      _vect_M.resize(2);
      _vect_M[0] = -_l*cos(_phi)/2;
      _vect_M[1] = _l*sin(_phi)/2;

      double R, Ri, X;
      R = sqrt(pow((_vect_M[0]-_vect_s[0]),2)+pow((_vect_M[1]-_vect_s[1]),2));
      Ri = sqrt(pow((_vect_M[0]-_vect_si[0]),2)+pow((_vect_M[1]-_vect_si[1]),2));
      X = _vect_M[0]-_vect_s[0];

      double K0, K0i;
      K0 = boost::math::cyl_bessel_k(0, (R*_V_c/(2*a)));
      K0i = boost::math::cyl_bessel_k(0, (Ri*_V_c/(2*a)));

      Q += exp(-X*_V_c/(2*a))/(K0+K0i);
  }

  _T_M = _Q_thq_zs  / (2*_pi*K_th) * Q*_l/n;


}

void OrthogonalMachining::SaveFileZCP(const std::string file_name)
{

}

void OrthogonalMachining::ResultsZCS()
{
  _delta_t2 = 0.3*(_f/_z)*cos(_phi - _gamma_0)/sin(_phi);
  _delta_lc = 0.9*(cos(_phi - _gamma_0)/sin(_phi)+2)*(sin(_phi-_gamma_0)+(0.3*cos(_phi-_gamma_0))/sin(_phi));
  /// Aparemment le 0.9 est = f/2, f l'avance
}

void OrthogonalMachining::SaveFileZCS(const std::string file_name)
{
}

void OrthogonalMachining::ResultsZD()
{
}

void OrthogonalMachining::SaveFileZD(const std::string file_name)
{
}



Milling::Milling(Material& material, Tool& tool, double phi) :
Machining(material, tool, phi)
{
    _a_p = 2;
    _kappa_r = _pi/4;
    _lambda_s = 6*_pi/180;
    _theta = 3 * _pi /2;
}

/// Calculs dont les résultats sont utiles dans les trois zones
void Milling::Initialize()
{
  _hP2 = (_a_p/sin(_kappa_r)*cos(_lambda_s) - _rbec);
  _V_f = _f*_omega/(2.*_pi);

  _V_P_R1.resize(3);
  _V_P1_R1.resize(3);
  _V_P2inf_R1.resize(3);
  _V_P2sup_R1.resize(3);

  _V_P_R1[0] = _V_f*cos(_theta)-(_R_0*_omega);
  _V_P_R1[1] = _V_f*sin(_theta);
  _V_P_R1[2] = 0.;

  _V_P1_R1[0] = _V_f*cos(_theta)-((_R_0-_rbec*(2-sqrt(2))/2)*_omega);
  _V_P1_R1[1] = _V_f*sin(_theta);
  _V_P1_R1[2] = 0;

  _V_P2sup_R1 = _V_P_R1;
  _V_P2inf_R1 = _V_P_R1;

  _V_P_R4.resize(3);
  _V_P1_R4.resize(3);
  _V_P2inf_R4.resize(3);
  _V_P2sup_R4.resize(3);

  _V_P_R4[0] = (_V_f*cos(_theta)-(_R_0*_omega))*cos(_lambda_s)+_V_f*sin(_lambda_s)*cos(_kappa_r)*sin(_theta);
  _V_P_R4[1] = _V_f*sin(_theta)*sin(_kappa_r);
  _V_P_R4[2] = sin(_lambda_s)*(_V_f*cos(_theta)-_R_0*_omega)-_V_f*sin(_lambda_s)*cos(_kappa_r)*sin(_theta);

  _V_P1_R4 = _V_P_R4;
  _V_P1_R4[0] += _rbec*(2.-sqrt(2.))/2.*_omega*(cos(_gamma_0)*cos(_lambda_s)*sin(_kappa_r)-cos(_kappa_r));
  _V_P1_R4[1] -= _rbec*(2.-sqrt(2.))/2.*_omega* \
  (sin(_gamma_0)*cos(_lambda_s)*sin(_kappa_r)+sin(_kappa_r)*sin(_lambda_s));
  _V_P1_R4[2] += _rbec*(2.-sqrt(2.))/2.*_omega* \
  (cos(_gamma_0)*sin(_lambda_s)*sin(_kappa_r)+sin(_gamma_0)*cos(_kappa_r));

  _V_P2sup_R4 = _V_P_R4;
  _V_P2sup_R4[0] -= (_rbec + _hP2)*cos(_kappa_r)*_omega;
  _V_P2sup_R4[1] -= (_rbec + _hP2)*sin(_kappa_r)*sin(_lambda_s)*_omega;

  _V_P2inf_R1 = _V_P_R1;
  _V_P2inf_R4[0] -= _rbec*cos(_kappa_r)*_omega;
  _V_P2inf_R4[1] -= _rbec*sin(_kappa_r)*sin(_lambda_s)*_omega;

  _V_corth_P1 = sqrt(pow(_V_P1_R4[0],2)+pow(_V_P1_R4[1],2));
  _V_corth_P2sup = sqrt(pow(_V_P2sup_R4[0],2)+pow(_V_P2sup_R4[1],2));



  _a_pp = _a_p / (sin(_kappa_r)*cos(_lambda_s)) - _rbec*sin(_kappa_r);

  _V_N_P2inf = sqrt(pow(_V_P2inf_R4[0],2)+pow(_V_P2inf_R4[1],2))*sin(_phi);
  _V_N_P2sup = sqrt(pow(_V_P2sup_R4[0],2)+pow(_V_P2sup_R4[1],2))*sin(_phi);
  _dV_N = _V_N_P2sup - _V_N_P2inf;

  _V_S_P2inf = sqrt(pow(_V_P2inf_R4[0],2)+pow(_V_P2inf_R4[1],2))*cos(_gamma_0)/sin(_phi-_gamma_0);
  _V_S_P2sup = sqrt(pow(_V_P2sup_R4[0],2)+pow(_V_P2sup_R4[1],2))*cos(_gamma_0)/sin(_phi-_gamma_0);
  _dV_S = _V_S_P2sup - _V_S_P2inf;

  _V_gz = 1.;
  _dV_z = 1.;

  _V_cop_P2inf = sqrt(pow(_V_P2inf_R4[0],2)+pow(_V_P2inf_R4[1],2))*sin(_phi)/cos(_phi-_gamma_0);
  _V_cop_P2sup = sqrt(pow(_V_P2sup_R4[0],2)+pow(_V_P2sup_R4[1],2))*sin(_phi)/cos(_phi-_gamma_0);
  _dV_cop = _V_cop_P2sup - _V_cop_P2inf;

  _V_y_P2inf = sqrt(pow(_V_P2inf_R4[0],2)+pow(_V_P2inf_R4[1],2))*cos(_gamma_0);
  _V_y_P2sup = sqrt(pow(_V_P2sup_R4[0],2)+pow(_V_P2sup_R4[1],2))*cos(_gamma_0);
  _dV_y = _V_y_P2sup - _V_y_P2inf;


}


void Milling::ResultsZCP()
{

  _vect_phi.resize(3);
  _vect_phi[0] = 1.;
  _vect_phi[1] = 1.;
  _vect_phi[2] = 1.;
  double t = _h_moy / _V_N_P2inf;   /// Temps pris en compte dans les calculs (élément de volume parcourt la zone)

  _U_phi.resize(3);
  _U_phi[0] = (_dV_N / (_a_pp*_l) * _vect_phi[1]*_vect_phi[2] + _V_N_P2inf / _l * _vect_phi[1])*t;
  _U_phi[1] = (_dV_S / (_a_pp*_h_moy) * _vect_phi[0]*_vect_phi[2] + _V_S_P2inf / _h_moy * _vect_phi[0])*t;
  _U_phi[2] = 1.;


  double k1, k2, k3;
  k1 = pow((_dV_N / (_a_pp*_l)*_vect_phi[2]+_V_N_P2inf/_l +_dV_S / (_a_pp*_h_moy)*_vect_phi[2]+_V_S_P2inf/_h_moy),2);
  k2 = pow((_dV_N / (_a_pp*_l)*_vect_phi[1]+_dV_z / (_l*_h_moy) * _vect_phi[1] +_V_gz / _h_moy),2);
  k3 = pow((_dV_S / (_a_pp*_h_moy)*_vect_phi[0]+_dV_z / (_l*_h_moy) * _vect_phi[0]),2);
  _eps_eq_zp = t*sqrt(1./3 *(k1 + k2 + k3));

  _d_eps_eq_zp = _eps_eq_zp / t;

  _material.JohnsonCookModel(_eps_eq_zp, _d_eps_eq_zp, _strain_rate_0);
  _sigm_eq_zp = _material.GetEquivalentStress();
  _Q_thq_zp = _sigm_eq_zp * _d_eps_eq_zp * _a_pp * _h_moy * _l;
}

void Milling::SaveFileZCP(const std::string file_name)
{
  _file = file_name + "_ZCP.txt";
  _file_out_ZCP.open(_file);
  for (int i=0 ; i < 3 ; i++)
  {
      _file_out_ZCP << _U_phi[i] << " ";
   }
  _file_out_ZCP << _pi;


  _file_out_ZCP.close();
}


void Milling::ResultsZCS()
{
  _delta_t2 = 0.3*(_f/_z)*cos(_phi - _gamma_0)/sin(_phi);
  _delta_lc = 0.9*(cos(_phi - _gamma_0)/sin(_phi)+2)*(sin(_phi-_gamma_0)+(0.3*cos(_phi-_gamma_0))/sin(_phi));

  _vect_c.resize(3);
  _vect_c[0] = 1.;
  _vect_c[1] = 1.;
  _vect_c[2] = 1.;
  double t = _delta_lc / _V_cop_P2inf;     /// Temps pris en compte dans les calculs (élément de volume parcourt la zone)

  _U_c.resize(3);
  _U_c[0] = (-_dV_cop / (_a_pp*_delta_lc) * _vect_c[0]*_vect_c[2] - _V_cop_P2inf / _delta_lc * _vect_phi[0] \
    +_dV_cop / _a_pp * _vect_c[2]+_V_cop_P2inf)*t;
  _U_c[1] = (_dV_y / _a_pp * (1-_vect_c[0]/_delta_lc) * _vect_c[2])*t;
  _U_c[2] = (_dV_z / _delta_t2 * (_vect_c[1] - _vect_c[0]*_vect_c[1] / _delta_lc))*t ;

    double k1, k2, k3, k4;
    k1 = pow(_dV_y / (_a_pp*_delta_lc)*_vect_c[2],2);
    k2 = 4./3*pow(_dV_cop / (_a_pp*_delta_lc)*_vect_c[2]+_V_cop_P2inf/_delta_lc,2);
    k3 = pow((-_dV_cop / (_a_pp*_delta_lc)*_vect_c[0]+_dV_cop / _a_pp - _dV_z / (_delta_lc*_delta_t2) * _vect_c[1] ),2);
    k4 = pow(_dV_y / _a_pp*_h_moy * (1 - _vect_c[0]/_delta_lc)+_dV_z / (_delta_t2) * (1 - _vect_c[0]/_delta_lc),2);

    _eps_eq_zs = t*sqrt(1./3 * (k1 + k2 +k3 + k4));

   _d_eps_eq_zs = 1.;

  _material.JohnsonCookModel(_eps_eq_zs, _d_eps_eq_zs, _strain_rate_0);
  _sigm_eq_zs = _material.GetEquivalentStress();
  _Q_thq_zs = _sigm_eq_zs * _d_eps_eq_zs * _a_pp * _delta_t2 * _delta_lc;


}

void Milling::SaveFileZCS(const std::string file_name)
{
  _file = file_name + "_ZCS.txt";
  _file_out_ZCS.open(_file);
  for (int i=0 ; i < 3 ; i++)
  {
      _file_out_ZCS << _U_c[i] << " ";
  }
  _file_out_ZCS.close();
}

void Milling::ResultsZD()
{
}

void Milling::SaveFileZD(const std::string file_name)
{
}
