/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=-1.0;
    AB[kk+ *kv+1]=2.0;
    AB[kk+ *kv+2]=-1.0;
  }
  AB[0]=0.0;
  if (*kv == 1) {AB[1]=0;}
  
  AB[(*lab)*(*la)-1]=0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}  

void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}
/**
 * @brief Calcule l'erreur relative directe entre deux vecteurs.
 * @param x   Pointeur vers le vecteur de la solution approchée de dimension *la.
 * @param y   Pointeur vers le vecteur de la solution de référence de dimension *la.
 * @param la  Pointeur vers le nombre de points (inconnues) de la grille 1D (taille de la matrice carrée la x la).
 * @return La valeur de l'erreur relative.
 */
double relative_forward_error(double* x, double* y, int* la) {
  double num = 0.0;
  double den = 0.0;
  for (int i = 0; i < *la; i++) {
    num += (x[i] - y[i]) * (x[i] - y[i]);  // Norme 2 du vecteur des erreurs
    den += y[i] * y[i];                    // Norme 2 de la solution de référence (on prend y car dans l'appel de la fonction l'argument y est EX_SOL)
  }
  return sqrt(num) / sqrt(den);
}

int indexABCol(int i, int j, int* lab) {
  return 0;
}

/**
 * @brief Effectue la factorisation LU d'une matrice tridiagonale au format bande (GB).
 * @param la     Nombre de points (inconnues) de la grille 1D (taille de la matrice carrée la x la).
 * @param n      Taille de la matrice A (logiquement la même que la).
 * @param kl     Nombre de sous-diagonales (pour une tridiagonale, kl = 1).
 * @param ku     Nombre de sur-diagonales (pour une tridiagonale, ku = 1).
 * @param AB     Pointeur vers le tableau de la matrice bande au format col-major (dimension (lab x la)).
 * @param lab    Nombre total de bandes significatives (par exemple 3 pour une tridiagonale).
 * @param ipiv   Tableau de pivots de la factorisation (non utilisé ici, mais demandé pour la compatibilité LAPACK).
 * @param info   Indicateur de succès de la factorisation.
 * @return info  Retourne 0 si la factorisation réussit, sinon une erreur.
 */
int dgbtrftridiag(int* la, int* n, int* kl, int* ku, double* AB, int* lab, int* ipiv, int* info) {
  *info = 0;
  for (int k = 0; k < (*la) - 1; k++) {
    // Vérification du pivot
    double pivot = AB[1 + k * (*lab)];
    if (pivot == 0.0) {
      *info = k + 1;
      return *info;  // Retourne immédiatement l'erreur
    }

    // Stockage de la permutation (pas utile ici, mais évite le segmentation fault)
    ipiv[k] = k + 1;

    if ((k + 1) < *la) {
      // 1. Mise à jour des coefficients de L (bande sous-diagonale)
      double l = AB[0 + (k + 1) * (*lab)] / pivot;
      AB[0 + (k + 1) * (*lab)] = l;

      // 2. Mise à jour des coefficients de U (diagonale principale et bande au-dessus)
      AB[1 + (k + 1) * (*lab)] -= l * AB[2 + k * (*lab)];
    }
  }
  return *info;
}
