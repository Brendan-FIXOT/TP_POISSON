/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

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
