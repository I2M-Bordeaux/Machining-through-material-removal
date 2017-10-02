#include <fstream>
#include <vector>
#include "DataFile.h"

class Tool
{
private:
    double _V_c;
    double _f;
    double _omega;
    double _z;
    double _R_0;
    double _r_beta;
    double _gamma_0;
    double _r_eps;
    double _kappa_r;
    double _lambda_s;



  public:
    // Constructeur par défaut
    Tool(DataFile* data_file);
    double GetV_c(){return _V_c;}
    double GetGamma(){return _gamma_0;}
    double GetTeeth(){return _z;}
    double GetFeed(){return _f;}
    double GetAngularVelocity(){return _omega;}
    double GetToolEdgeRadius(){return _R_0;}

};
