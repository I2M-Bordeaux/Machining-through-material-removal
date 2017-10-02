#include "Material.h"
#include "Tool.h"



#include <fstream>
#include <vector>

class Machining
{
  protected:

    double _pi;

    /// Variables cinématiques
    double _gamma_0, _phi; /// Angle de coupe, angle de cisaillement primaire
    double  _f, _omega; /// avance, vitesse de rotation
    double _V_f, _V_c;  ///Vitesse d'avance, vitesse de coupe
    double _l; /// Longueur de la ZCP
    double _delta_t2, _delta_lc; /// Longueur de contact outil matière, épaisseur de la ZCS
    double _strain_rate_0;

    ///Pointeur vers le matériau choisi
    Material& _material;
    ///Pointeur vers l'outil choisi
    Tool& _tool;

    int _z; /// Nombre de dents
    double _R_0, _rbec; /// Rayon de l'outil, Rayon de bec


    /// Écriture du fichier
    std::string _file;
    std::ofstream _file_out_ZCP, _file_out_ZCS;

    std::vector<double> _vect_phi ; ///Base du repère lié à la ZCP
    double _eps_eq_zp; /// Grandeur généralisée de la def dans la ZCP
    double _d_eps_eq_zp;    ///Vitesse de déformation généralisée dans la ZCP
    double _sigm_eq_zp; /// Contrainte équivalent dans la ZCP
    double _Q_thq_zp;   /// Puissance thermique de la ZCP

    double _eps_eq_zs; /// Grandeur généralisée de la def dans la ZCS
    double _d_eps_eq_zs;    ///Vitesse de déformation généralisée dans la ZCS
    double _sigm_eq_zs; /// Contrainte équivalent dans la ZCS
    double _Q_thq_zs;   /// Puissance thermique de la ZCS

    double _test;


  public:
    // Constructeur par défaut
    Machining(Material& material, Tool& tool, double phi);
    // Destructeur par défaut
    virtual ~Machining();



    double Minimize();


    ///Initialisation, calculs préliminaires (champs cinématiques ...)
    virtual void Initialize() = 0;
    // Calculs actions mécaniques (déplacements, vitesse de déplacement, déformation)
    virtual void ResultsZCP() = 0;
    // Sauvegarde les solutions dans un fichier txt
    virtual void SaveFileZCP(const std::string file_name) = 0;
    // Calculs des vitesses, déplacements, vitesse de déplacement, déformation
    virtual void ResultsZCS() = 0;
    // Sauvegarde les solutions dans un fichier txt
    virtual void SaveFileZCS(const std::string file_name) = 0;
    // Calculs des vitesses, déplacements, vitesse de déplacement, déformation
    virtual void ResultsZD() = 0;
    // Sauvegarde les solutions dans un fichier txt
    virtual void SaveFileZD(const std::string file_name) = 0;
};

class OrthogonalMachining : public Machining
{
  private:

    double _V_cop;   ///Vitesse du copeau par Merchant
    double _V_s;

    double _T_M;

    std::vector<double> _vect_s ;
    std::vector<double> _vect_si ;
    std::vector<double> _vect_M ;

  public:
    OrthogonalMachining(Material& material, Tool& tool, double phi);

    void Initialize();
    void ResultsZCP();
    void SaveFileZCP(const std::string file_name);
    void ResultsZCS();
    void SaveFileZCS(const std::string file_name);
    void ResultsZD();
    void SaveFileZD(const std::string file_name);
};


class Milling : public Machining
  {

  private:

    /// Variables cinématiques
    double _kappa_r, _lambda_s ;  /// Angle d'attaque, Vitesse de rotation de fraise, angle d'inclinaison d'arête
    double _theta, _a_p ; /// Angle de rotation de la fraise

    std::vector<double> _V_P_R1 ;   ///Vitesse du point P dans le repère R1
    std::vector<double> _V_P1_R1 ;   ///Vitesse du point P1 dans le repère R1
    std::vector<double> _V_P2inf_R1 ;   ///Vitesse du point P2_inf dans le repère R1
    std::vector<double> _V_P2sup_R1 ;   ///Vitesse du point P2_sup dans le repère R1
    std::vector<double> _V_P_R4 ;   ///Vitesse du point P dans le repère R4
    std::vector<double> _V_P1_R4 ;   ///Vitesse du point P1 dans le repère R4
    std::vector<double> _V_P2inf_R4 ;   ///Vitesse du point P2_inf dans le repère R4
    std::vector<double> _V_P2sup_R4 ;   ///Vitesse du point P2_sup dans le repère R4
    double _a_pp, _hP2; ///profondeur de passe réelle, fonction linéaire
    double _V_corth_P1, _V_corth_P2sup; ///Vitesse de coupe
    double _V_N_P2inf, _V_N_P2sup, _dV_N;
    double _V_S_P2inf, _V_S_P2sup, _dV_S;
    double _V_gz, _dV_z;
    double _V_cop_P2inf, _V_cop_P2sup, _dV_cop;
    double _V_y_P2inf, _V_y_P2sup, _dV_y;



    double _h_moy; /// Epaisseur de la ZCP

    std::vector<double> _U_phi ;        ///Champ de déplacement

    std::vector<double> _vect_c ; ///Base du repère lié à la ZCS
    std::vector<double> _U_c ;        ///Champ de déplacement



  public:
    Milling(Material& material, Tool& tool, double phi);
    void Initialize();
    void ResultsZCP();
    void SaveFileZCP(const std::string file_name);
    void ResultsZCS();
    void SaveFileZCS(const std::string file_name);
    void ResultsZD(); /// Pour l'instant cette fonction ne renvoie rien
    void SaveFileZD(const std::string file_name);
};
