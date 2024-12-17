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

/**
 * @brief Initialise le second membre (RHS) pour le problème de Poisson 1D avec conditions de Dirichlet.
 * @param RHS   Pointeur vers le tableau (de dimension la) contenant le second membre (RHS) à construire.
 * @param la    Pointeur vers le nombre de points (inconnues) de la grille 1D (taille de la matrice carrée la x la).
 * @param BC0   Pointeur vers la valeur de la condition de Dirichlet au bord gauche (u(0) = BC0).
 * @param BC1   Pointeur vers la valeur de la condition de Dirichlet au bord droit (u(1) = BC1).
 */
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  // Initialisation à zéro
  for (int i = 0; i < *la; i++) {
    RHS[i] = 0.0;
  }
  // Condition au bord de DIRICHLET
  RHS[0] = *BC0; // Gauche
  RHS[*la - 1] = *BC1; // Droit
}

/**
 * @brief Calcule la solution analytique de l'équation de Poisson 1D avec conditions de Dirichlet.
 * @param EX_SOL  Pointeur vers le tableau contenant la solution analytique de dimension la.
 * @param X       Pointeur vers le tableau contenant les positions des points de la grille de dimension la.
 * @param la      Pointeur vers le nombre de points (inconnues) de la grille 1D (taille de la matrice carrée la x la).
 * @param BC0     Pointeur vers la valeur de la condition de Dirichlet au bord gauche (u(0) = BC0).
 * @param BC1     Pointeur vers la valeur de la condition de Dirichlet au bord droit (u(1) = BC1).
 */
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  for (int i = 0; i < *la; i++) {
    *EX_SOL[i] = *BC0 + X[i] * (*BC1 - *BC0); // Solution linéaire de l'énoncé
  }
}

/**
 * @brief Calcule les points de la grille 1D uniformément espacés sur l'intervalle [0, 1].
 * @param X   Pointeur vers le tableau contenant les positions des points de la grille de dimension la.
 * @param la  Pointeur vers le nombre de points (inconnues) de la grille 1D (taille de la matrice carrée la x la).
 */
void set_grid_points_1D(double* x, int* la){
  double h = 1.0 / (*la + 1); // pas de la grille
  for (int i = 0; i < *la; i++) {
    x[i] = (i + 1) * h;
  }
}

/**
 * @brief Calcule l'erreur relative directe entre deux vecteurs.
 * @param x   Pointeur vers le vecteur de la solution approchée de dimension *la.
 * @param y   Pointeur vers le vecteur de la solution de référence de dimension *la.
 * @param la  Pointeur vers le nombre de points (inconnues) de la grille 1D (taille de la matrice carrée la x la).
 * @return La valeur de l'erreur relative.
 */
double relative_forward_error(double* x, double* y, int* la){
  double num = 0.0;
  double den = 0.0;
  for (int i = 0; i < *la; i++) {
    num += (x[i] - y[i]) * (x[i] - y[i]); // Norme 2 du vecteur des erreurs
    den += y[i] * y[i]; // Norme 2 de la solution de référence (on prend y car dans l'appel de la fonction l'argument y est EX_SOL)
  }
  return sqrt(num) / sqrt(den);
}

int indexABCol(int i, int j, int *lab){
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
int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  *info = 0;

  for (int i = 0; i < *la - 1; i++) {    
    // 1. Mise à jour des coefficients de L (bande sous-diagonale)
    double l = AB[0 + (k + 1) * (*lab)] / AB[1 + k * (*lab)];
    AB[0 + (k + 1) * (*lab)] = l;

    // 2. Mise à jour des coefficients de U (diagonale principale et bande au-dessus)
    AB[1 + (k + 1) * (*lab)] -= l * AB[2 + k * (*lab)];
  }
  return *info;
}
