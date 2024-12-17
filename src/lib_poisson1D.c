/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

/**
 * @brief Initialise l'opérateur de Poisson 1D au format bande (GB) en stockage colonne par colonne (col-major).
 * @param AB    Pointeur vers le tableau de la matrice bande au format col-major (dimension (lab x la)).
 * @param lab   Nombre total de bandes significatives (par exemple 3 pour une tridiagonale).
 * @param la    Nombre de points (inconnues) de la grille 1D (taille de la matrice carrée la x la).
 * @param kv    Position de la diagonale principale dans la matrice bande (pour lab = 3, kv = 1).
 */
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  // Initialisation à zéro
  for (int i = 0; i < lab * la; i++) {
    AB[i] = 0.0;
  }

  // Remplissage de la diagonale principale et des bandes sous/sur diagonale
  for (int i = 0; i < la; i++) {
    // Bande sous-diagonale
    if (i > 0) {
      AB[(kv - 1) + i * lab] = -1.0;
    }
    // Diagonale principale
    AB[(kv) + i * lab] = 2.0;
    // Bande sur-diagonale
    if (i < la - 1) {
      AB[(kv + 1) + i * lab] = -1.0;
    }
  }
}
/**
 * @brief Initialise la matrice identité au format bande (GB) en stockage colonne par colonne (col-major).
 * @param AB    Pointeur vers le tableau de la matrice bande au format col-major (dimension (lab x la)).
 * @param lab   Nombre total de bandes significatives (par exemple 3 pour une tridiagonale).
 * @param la    Nombre de points (inconnues) de la grille 1D (taille de la matrice carrée la x la).
 * @param kv    Position de la diagonale principale dans la matrice bande (pour lab = 3, kv = 1).
 */
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  // Initialisation à zéro
  for (int i = 0; i < (*lab) * (*la); i++) {
    AB[i] = 0.0;
  }
  // Remplissage de la diagonale principale avec des 1
  for (int i = 0; i < *la; i++) {
    AB[(*kv) + i * (*lab)] = 1.0;
  }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
}  

void set_grid_points_1D(double* x, int* la){
}

double relative_forward_error(double* x, double* y, int* la){
}

int indexABCol(int i, int j, int *lab){
  return 0;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}
