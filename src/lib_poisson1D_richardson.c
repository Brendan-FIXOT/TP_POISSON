/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

/**
 * @brief Calcule les valeurs propres de la matrice associée au problème 1D de Poisson.
 *
 * Les valeurs propres de la matrice de Poisson 1D sont données par :
 * λ_k = 4 * sin²(kπh / 2), où k = 1, 2, ..., la.
 *
 * @param la      Pointeur vers le nombre de points (inconnues) de la grille 1D.
 * @param eigval  Pointeur vers les valeurs propres de la matrice associée au problème 1D de Poisson.
 */
void eig_poisson1D(int* la, double *eigval) {
  double h = 1.0 / (*la + 1);  // pas de discrétisation
  double sin_val;              // sin(kπh / 2), où k = 1, 2, ..., la
  for (int k = 1; k <= (*la); ++k) {
    sin_val = sin(k * M_PI * h / 2);
    eigval[k - 1] = 4 * sin_val * sin_val;  // Calcul de la valeur propre pour le k correspondant
    //printf("Valeur propre %d : %lf\n", k, eigval[k - 1]);
  }
}

/**
 * @brief Calcule la valeur propre maximale de la matrice associée au problème 1D de Poisson.
 *
 * Le tableau eigval est parcouru pour déterminer la valeur propre maximale
 *
 * @param la      Pointeur vers le nombre de points (inconnues) de la grille 1D.
 * @param eigval  Pointeur vers les valeurs propres de la matrice associée au problème 1D de Poisson.
 * @return        Return la valeur propre maximale.
 */
double eigmax_poisson1D(int* la, double *eigval) {
  double eigmax = eigval[0];
  for (int i = 1; i < *la; ++i) {
    if (eigmax < eigval[i]) {
      eigmax = eigval[i];
    }
  }
  return eigmax;
}

/**
 * @brief Calcule la valeur propre minimale de la matrice associée au problème 1D de Poisson.
 *
 * Le tableau eigval est parcouru pour déterminer la valeur propre minimale
 *
 * @param la      Pointeur vers le nombre de points (inconnues) de la grille 1D.
 * @param eigval  Pointeur vers les valeurs propres de la matrice associée au problème 1D de Poisson.
 * @return        Retourne la valeur propre minimale.
 */
double eigmin_poisson1D(int* la, double *eigval) {
  double eigmin = eigval[0];  // Les valeurs propres de la matrice ne sont pas négatives
  for (int i = 1; i < *la; ++i) {
    if (eigmin > eigval[i]) {
      eigmin = eigval[i];
    }
  }
  return eigmin;
}

/**
 * @brief     Calcule la valeur optimale d'alpha pour la fonction richardson_alpha()
 *
 * alpha_opt = 2 / (lambdamin + lambdamax)
 *
 * @param la  Pointeur vers le nombre de points (inconnues) de la grille 1D (taille de la matrice carrée la x la).
 * @return    Returne la valeur optimale d'alpha
 */
double richardson_alpha_opt(int* la, double *eigval) {
  double eigmax = eigmax_poisson1D(la, eigval);
  double eigmin = eigmin_poisson1D(la, eigval);

  if (eigmax + eigmin == 0) {
    fprintf(stderr, "Erreur : Somme des valeurs propres maximale et minimale est nulle.\n");
    exit(EXIT_FAILURE);
  }
  
  return 2.0 / (eigmax + eigmin);
}

/**
 * @brief Résout le système linéaire AX = B en utilisant l'algorithme de Richardson.
 *
 * @param AB          Pointeur vers le tableau de la matrice bande au format col-major (dimension (lab x la)).
 * @param RHS         Pointeur vers le tableau (de dimension la) contenant le second membre (RHS) à construire.
 * @param X           Pointeur vers le tableau contenant les positions des points de la grille de dimension la.
 * @param alpha_rich  Pas de relaxation alpha.
 * @param lab         Pointeur vers le nombre total de bandes significatives (par exemple 3 pour une tridiagonale).
 * @param la          Pointeur vers le nombre de points (inconnues) de la grille 1D (taille de la matrice carrée la x la).
 * @param ku          Nombre de sur-diagonales.
 * @param kl          Nombre de sous-diagonales.
 * @param tol         Tolérance pour le critère de convergence basé sur la norme du résidu.
 * @param maxit       Nombre maximal d'itérations.
 * @param resvec      Tableau pour stocker l'historique de la norme du résidu à chaque itération.
 * @param nbite       Nombre total d'itérations effectuées avant d'atteindre la convergence.
 */
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
  double res_norm;

  // Allocation du résidu temporaire
  double *res = (double *)malloc((*la) * sizeof(double));
  if (res == NULL) {
    perror("Erreur d'allocation mémoire pour le vecteur res");
    exit(EXIT_FAILURE);
  }

  // Norme de RHS pour normaliser le résidu
  double norm_rhs = cblas_dnrm2(*la, RHS, 1);

  for (int iter = 0; iter < *maxit; iter++) {
    // Initialisation de res avec B (sans memcpy)
    /*
    for (int i = 0; i < *la; i++) {
      res[i] = RHS[i];
    }
    */
    // Initialisation de res avec B (avec memcpy)
    memcpy(res, RHS, (*la) * sizeof(double));

    // Calcul du résidu res = B - AX
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, res, 1);

    // Calcul de la norme du résidu (sans cblas)
    /*
    res_norm = 0.0;
    for (int i = 0; i < *la; i++) {
      res_norm += res[i] * res[i];
    }
    res_norm = sqrt(res_norm);
    */
    // Calcul de la norme relative du résidu (avec cblas)
    res_norm = cblas_dnrm2(*la, res, 1) / norm_rhs;
    resvec[iter] = res_norm;

    // Vérification de convergence
    if (res_norm < *tol) {
      *nbite = iter + 1;
      break;
    }

    // Mise à jour de X : X = X + alpha * res
    cblas_daxpy(*la, *alpha_rich, res, 1, X, 1);
  }

  // Si convergence non atteinte, nbite = maxit
  if (res_norm >= *tol) {
    *nbite = *maxit;
  }

  // Libération de la mémoire
  free(res);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) {
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) {
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
}
