#include <string>
#include <vector>
#include <iostream>

class DataFile
{

private:
  std::string _file_name;
  double _Vc;
  double _f;
  double _gamma_0;
  double _z;
  double _rbec;
  double _omega;
  double _R_0;
  std::string _machining_choice;
  int _choice;


public:
  DataFile(std::string file_name);
  void ReadDataFile();
  std::string GetFileName(){return _file_name;}
  double GetV_c(){return _Vc;}
  double GetGamma(){return _gamma_0;}
  double GetTeeth(){return _z;}
  double GetFeed(){return _f;}
  double GetAngularVelocity(){return _omega;}
  double GetToolEdgeRadius(){return _R_0;}
  double GetChoice(){return _choice;}
  std::string GetMachiningChoice(){return _machining_choice;}

};
