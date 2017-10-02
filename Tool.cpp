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
/// \file Tool.cpp
/// \author Renault Chloe (chloe.renault@u-bordeaux.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include "Tool.h"
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// Constructeur par défaut
Tool::Tool(DataFile* data_file)
{
  _V_c = data_file->GetV_c();
  _V_c *= pow(10, 3);

  _f = data_file->GetFeed();
  _gamma_0 = data_file->GetGamma();
  _z = data_file->GetTeeth();

  _omega = data_file->GetAngularVelocity();
  _R_0 = data_file->GetToolEdgeRadius()*pow(10,-3);
}
