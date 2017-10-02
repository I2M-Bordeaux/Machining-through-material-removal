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
/// \file DataFile.cpp
/// \author Renault Chloe (chloe.renault@u-bordeaux.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include "DataFile.h"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

DataFile::DataFile(std::string file_name): _file_name(file_name)
{}

void DataFile::ReadDataFile()
{
  ifstream data_file(_file_name.data());
  if (!data_file.is_open())
  {
    cout << "Unable to open file" << _file_name << endl;
    abort();
  }

  string file_line;

  while (!data_file.eof())
  {
    getline(data_file, file_line);

    if (file_line.find("Vc") != std::string::npos)
    {
      data_file >> _Vc;
    }

    if (file_line.find("gamma_0") != std::string::npos)
    {
      data_file >> _gamma_0;
    }

    if (file_line.find("R0") != std::string::npos)
    {
      data_file >> _R_0;
    }

    if (file_line.find("rbec") != std::string::npos)
    {
      data_file >> _rbec;
    }

    if (file_line.find("_f") != std::string::npos)
    {
      data_file >> _f;
    }

    if (file_line.find("teeth") != std::string::npos)
    {
      data_file >> _z;
    }

    if (file_line.find("choice") != std::string::npos)
    {
      data_file >> _machining_choice;
      if (_machining_choice != "orthogonal")
      {
        _choice = 1;
      }
      if (_machining_choice != "milling")
      {
        _choice = 2;
      }
    }
  }
}
