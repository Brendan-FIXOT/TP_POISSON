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
void eig_poisson1D(int *la, double *eigval) {
  double h = 1.0 / (*la + 1);  // pas de discrétisation
  double sin_val;              // sin(kπh / 2), où k = 1, 2, ..., la
  for (int k = 1; k <= (*la); ++k) {
    sin_val = sin(k * M_PI * h / 2);
    eigval[k - 1] = 4 * sin_val * sin_val;  // Calcul de la valeur propre pour le k correspondant
    // printf("Valeur propre %d : %lf\n", k, eigval[k - 1]);
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
double eigmax_poisson1D(int *la, double *eigval) {
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
double eigmin_poisson1D(int *la, double *eigval) {
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
double richardson_alpha_opt(int *la, double *eigval) {
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

  save_convergence_history("data/convergence_history_richardson.dat", resvec, *nbite);

  // Libération de la mémoire
  free(res);
}

/**
 * @brief Extrait une matrice tridiagonale au format MB (Matrice Bande) à partir d'une matrice bande générale au format AB.
 *
 * La matrice MB est définie comme (M = D), où :
 * D est la diagonale principale.
 *
 * @param AB  Pointeur vers la matrice bande d'origine au format col-major (dimension (lab x la)).
 * @param MB  Pointeur vers la matrice bande tridiagonale au format col-major (dimension (lab x la)).
 * @param lab Nombre total de bandes significatives dans la matrice AB.
 * @param la  Taille de la matrice (nombre d'inconnues, c'est-à-dire la taille de la grille).
 * @param ku  Nombre de sur-diagonales dans la matrice bande AB.
 * @param kl  Nombre de sous-diagonales dans la matrice bande AB.
 * @param kv  Indice de la diagonale principale dans la matrice bande AB.
 */
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) {
  // Copie uniquement de la diagonale principale de AB dans MB
  for (int i = 0; i < *la; i++) {
    MB[(*lab / 2) + i * (*lab)] = AB[*kv + i * (*lab)];
    // printf("lab = %d et lab/2 = %d\n", *lab, (*lab/2));
  }

  /*
  // Impression de la matrice AB (format GB)
  printf("Matrice AB au format bande (lab=%d, la=%d):\n", *lab, *la);
  for (int i = 0; i < *lab; i++) {
    for (int j = 0; j < *la; j++) {
      printf("%6.2f ", AB[i + j * (*lab)]);
    }
    printf("\n");
  }

  // Impression de la matrice MB
  printf("Matrice MB au format bande (lab=%d, la=%d):\n", *lab, *la);
  for (int i = 0; i < *lab; i++) {
    for (int j = 0; j < *la; j++) {
      printf("%6.2f ", MB[i + j * (*lab)]);
    }
    printf("\n");
  }
  */
}

/**
 * @brief Implémente la méthode de Jacobi pour résoudre un système linéaire en utilisant le format GB.
 *
 * @param AB          Pointeur vers la matrice bande d'origine au format col-major (dimension (lab x la)).
 * @param RHS         Pointeur vers le tableau contenant le second membre (RHS).
 * @param X           Pointeur vers le tableau contenant les positions des points de la grille de dimension la.
 * @param la          Taille de la matrice (nombre d'inconnues, c'est-à-dire la taille de la grille).
 * @param lab         Nombre total de bandes significatives dans MB (3 pour une matrice tridiagonale, mais seules les diagonales principales sont utilisées ici).
 * @param tol         Tolérance pour le critère de convergence.
 * @param maxit       Nombre maximal d'itérations.
 * @param resvec      Tableau pour stocker l'historique de la norme du résidu à chaque itération.
 * @param nbite       Pointeur vers le nombre total d'itérations effectuées avant convergence.
 */
void jacobi_GB(double *AB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
  double res_norm;

  // Allocation pour la nouvelle solution et le résidu temporaire
  double *X_new = (double *)malloc((*la) * sizeof(double));
  double *res = (double *)malloc((*la) * sizeof(double));
  if (X_new == NULL || res == NULL) {
    perror("Erreur d'allocation mémoire");
    exit(EXIT_FAILURE);
  }

  // Norme de RHS pour normaliser le résidu
  double norm_rhs = cblas_dnrm2(*la, RHS, 1);

  // Initialisation du résidu, puis calcul du résidu avant d'entrer dans la boucle pour avoir 1 en première valeur de resvec
  memcpy(res, RHS, (*la) * sizeof(double));
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, res, 1);
  res_norm = cblas_dnrm2(*la, res, 1) / norm_rhs;

  resvec[0] = res_norm;

  // Itérations de Jacobi
  for (int iter = 0; iter < *maxit; iter++) {
    res_norm = 0.0;

    for (int i = 0; i < *la; i++) {
      double diag = AB[(*lab / 2) + i * (*lab)];  // Diagonale principale
      double sum = RHS[i];

      // Contribution de la sous-diagonale
      if (i > 0) {
        sum -= AB[((*lab / 2) - 1) + i * (*lab)] * X[i - 1];
      }

      // Contribution de la sur-diagonale
      if (i < *la - 1) {
        sum -= AB[((*lab / 2) + 1) + i * (*lab)] * X[i + 1];
      }

      // Mise à jour de X_new[i]
      X_new[i] = sum / diag;
    }

    // Calcul du résidu : res = RHS - AB * X_new (avec CBLAS)
    memcpy(res, RHS, (*la) * sizeof(double));  // Initialisation avec RHS
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X_new, 1, 1.0, res, 1);

    // Calcul de la norme du résidu avec CBLAS
    res_norm = cblas_dnrm2(*la, res, 1) / norm_rhs;
    resvec[iter + 1] = res_norm;

    // Vérification du critère de convergence
    if (res_norm < *tol) {
      *nbite = iter + 1;
      memcpy(X, X_new, (*la) * sizeof(double));  // Mise à jour finale de X
      break;
    }

    // Mise à jour pour l'itération suivante
    memcpy(X, X_new, (*la) * sizeof(double));
  }

  if (res_norm >= *tol) {
    *nbite = *maxit;  // Convergence non atteinte
  }

  save_convergence_history("data/convergence_history_jacobi.dat", resvec, *nbite + 1);

  // Libération des ressources allouées
  free(X_new);
  free(res);
}

/**
 * @brief Extrait la matrice MB pour la méthode de Gauss-Seidel à partir de la matrice AB au format tridiagonal.
 *
 * La matrice MB est définie comme (M = D - E), où :
 * - D est la diagonale principale.
 * - E est la sous-diagonale.
 *
 * @param AB  Pointeur vers la matrice bande originale au format GB (dimension lab x la).
 * @param MB  Pointeur vers la matrice bande extraite pour Gauss-Seidel (dimension lab x la).
 * @param lab Nombre total de bandes significatives dans AB.
 * @param la  Taille de la matrice carrée (nombre d'inconnues).
 * @param ku  Nombre de sur-diagonales dans AB.
 * @param kl  Nombre de sous-diagonales dans AB.
 * @param kv  Indice de la diagonale principale dans AB.
 */
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) {
  // Copie de la diagonale principale de AB dans MB
  for (int i = 0; i < *la; i++) {
    MB[(*lab / 2) + i * (*lab)] = AB[*kv + i * (*lab)];
  }

  // Copie de la sous-diagonale (E, négative pour M = D - E)
  if (*kl > 0) {
    for (int i = 0; i < *la; i++) {
      MB[((*lab / 2) + 1) + i * (*lab)] = -AB[(*kv + 1) + i * (*lab)];
    }
  }

  /*
  // Impression de la matrice AB (format GB)
  printf("Matrice AB au format bande (lab=%d, la=%d):\n", *lab, *la);
  for (int i = 0; i < *lab; i++) {
    for (int j = 0; j < *la; j++) {
      printf("%6.2f ", AB[i + j * (*lab)]);
    }
    printf("\n");
  }

  // Impression de la matrice MB
  printf("Matrice MB au format bande (lab=%d, la=%d):\n", *lab, *la);
  for (int i = 0; i < *lab; i++) {
    for (int j = 0; j < *la; j++) {
      printf("%6.2f ", MB[i + j * (*lab)]);
    }
    printf("\n");
  }
  */
}

void gauss_seidel_GB(double *AB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
  double res_norm;

  // Allocation pour le résidu temporaire
  double *res = (double *)malloc((*la) * sizeof(double));
  if (res == NULL) {
    perror("Erreur d'allocation mémoire pour le résidu");
    exit(EXIT_FAILURE);
  }

  // Norme de RHS pour normaliser le résidu
  double norm_rhs = cblas_dnrm2(*la, RHS, 1);

  
  // Initialisation du résidu, puis calcul du résidu avant d'entrer dans la boucle pour avoir 1 en première valeur de resvec
  memcpy(res, RHS, (*la) * sizeof(double));
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, res, 1);
  res_norm = cblas_dnrm2(*la, res, 1) / norm_rhs;

  resvec[0] = res_norm;

  // Boucle d'itérations
  for (int iter = 0; iter < *maxit; iter++) {
    // Mise à jour de X en place
    for (int i = 0; i < *la; i++) {
      double diag = AB[(*lab / 2) + i * (*lab)];
      double sum = RHS[i];

      // Contribution de la sous-diagonale
      if (i > 0) {
        sum -= AB[((*lab / 2) - 1) + i * (*lab)] * X[i - 1];
      }

      // Contribution de la sur-diagonale
      if (i < *la - 1) {
        sum -= AB[((*lab / 2) + 1) + i * (*lab)] * X[i + 1];
      }

      // Mise à jour de X[i]
      X[i] = sum / diag;
    }

    // Calcul du résidu res = RHS - AB * X
    memcpy(res, RHS, (*la) * sizeof(double));
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, res, 1);

    // Calcul de la norme relative du résidu
    res_norm = cblas_dnrm2(*la, res, 1) / norm_rhs;
    resvec[iter + 1] = res_norm;

    // Vérification du critère de convergence
    if (res_norm < *tol) {
      *nbite = iter + 1;
      break;
    }
  }

  // Si convergence non atteinte
  if (res_norm >= *tol) {
    *nbite = *maxit;
  }

  save_convergence_history("data/convergence_history_gauss_seidel.dat", resvec, *nbite + 1);

  // Libération de la mémoire allouée
  free(res);
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
}
