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
/// \file main.cc
/// \author Renault Chloe (chloe.renault@u-bordeaux.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include "Machining.h"
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
using namespace std;


int main(int argc, char** argv)
{
  if (argc < 2)
  {
    cout << "Please, enter the name of your data file." << endl;
    abort();
  }
  const string data_file_name = argv[1];

  ///Fichier de données
  DataFile* data_file = new DataFile(data_file_name);
  data_file->ReadDataFile();





  double pi(atan(1.)*4);

  /// Paramètre du matériau
  double A(595), B(768), C(0.0137), n(0.209), m(0.807); /// Constantes de la loi de comportement
  double Kth(48), a(15*pow(10, -6));

  Material* material = new Material(A, B, C, n, m, Kth, a);



  Tool* tool = new Tool(data_file);

  ///Paramètres cinématiques
  double phi(0.16); /// Angle de cisaillement primaire

  //double theta (1.); /// Angle de rotation de la fraise
  //double h(0.1); /// Epaisseur de la ZCP


  string results = data_file->GetMachiningChoice(); // nom du fichier résultat
  int userChoiceMachining = data_file->GetChoice(); // Choix de l’utilisateur

  Machining* mod = new OrthogonalMachining(*material, *tool, phi);

  phi = mod->Minimize();
  cout << phi <<endl;


  Machining* model(0);
  switch(userChoiceMachining)
  {
  case 1:
    model = new OrthogonalMachining(*material, *tool, phi);

    break;

    case 2:
    model = new Milling(*material, *tool, phi);

    break;
  }

  Machining& machining = *model;

  machining.Initialize();
  machining.ResultsZCP();
  machining.SaveFileZCP(results);
  machining.ResultsZCS();
  machining.SaveFileZCS(results);


  cout << "------------------------------------" << endl;
  cout << "Zone de cisaillement primaire" << endl;
  cout << "Résultats enregistrés dans le fichier " << results << "_ZCP.txt" << endl;
  cout << "------------------------------------" << endl;

  cout << "------------------------------------" << endl;
  cout << "Zone de cisaillement secondaire" << endl;
  cout << "Résultats enregistrés dans le fichier " << results << "_ZCS.txt" << endl;
  cout << "------------------------------------" << endl;

   if (userChoiceMachining == 1)
   {
     machining.ResultsZD();
     machining.SaveFileZD(results);
     cout << "------------------------------------" << endl;
     cout << "Zone de dépouille" << endl;
     cout << "Résultats enregistrés dans le fichier " << results << "_ZD.txt" << endl;
     cout << "------------------------------------" << endl;
  }



  delete model;
  delete mod;
  delete tool;
  delete material;


return 0;
}
