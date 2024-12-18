# **Exercice 3 : Référence et utilisation de BLAS/LAPACK**

## **1. Déclaration et allocation d'une matrice**

- Allouer la matrice au format colonne majeure :
  ```c
  double* A = (double*) malloc(n * m * sizeof(double));
  ```

## **2. Signification de LAPACK\_COL\_MAJOR**

- Indique que la matrice est stockée en **colonne majeure** (colonnes contiguës en mémoire).

## **3. Signification de la dimension principale (ld)**

- **ld** est le nombre d'éléments entre deux entrées successives d'une même colonne.
- Pour une matrice  de taille , **ld = n** en stockage colonne majeure.

## **4. Fonction dgbmv**

- Calcule :  pour une matrice bande **A**.

## **5. Fonction dgbtrf**

- Réalise la **factorisation LU** d'une matrice bande  sans permutation des lignes.

## **6. Fonction dgbtrs**

- Résout le système :  après une factorisation LU de .

## **7. Fonction dgbsv**

- Combine **dgbtrf** et **dgbtrs** pour résoudre directement .

## **8. Calcul de la norme du résidu relatif**

1. Calcul du résidu  avec **dgbmv**.

2. Calcul de la norme relative :

3. Utiliser **dnrm2** pour la norme 2.

---

# **Exercice 4 : Stockage GB et appel à DGBMV**

## **1. Stockage GB pour la matrice de Poisson 1D**

- La matrice de Poisson 1D tridiagonale :
  ```
  [  d  u  0  0  0 ]
  [  l  d  u  0  0 ]
  [  0  l  d  u  0 ]
  [  0  0  l  d  u ]
  [  0  0  0  l  d ]
  ```
- Stockage GB (3 lignes et  colonnes) :
  - Ligne 0 : sur-diagonale (u)
  - Ligne 1 : diagonale principale (d)
  - Ligne 2 : sous-diagonale (l)

## **2. Utilisation de dgbmv**

- Appel de **dgbmv** :
  ```c
  dgbmv("N", &n, &n, &kl, &ku, &alpha, AB, &ldab, x, &incx, &beta, y, &incy);
  ```

## **3. Méthode de validation**

- Comparer **dgbmv** avec le produit classique .
- Calculer la norme de la différence avec **dnrm2**.

---

# **Exercice 5 : DGBTRF, DGBTRS, DGBSV**

## **1. Résolution du système linéaire avec LAPACK**
- Utiliser **dgbsv** directement pour résoudre \(A \cdot X = B\) :
  ```c
  dgbsv(&n, &kl, &ku, &nrhs, AB, &ldab, ipiv, B, &ldb, &info);
  ```

---

## **2. Évaluation des performances et complexité**
- **Complexité de dgbtrf** : \(\mathcal{O}(n \cdot (kl + ku)^2)\)
- **Complexité de dgbtrs** : \(\mathcal{O}(n \cdot (kl + ku))\)
- La complexité dépend fortement de la largeur de la bande \((kl + ku)\).

---

# **Exercice 6 : LU pour les matrices tridiagonales**

## **1. Implémentation de la factorisation LU**
- Soit \(A\) une matrice tridiagonale, la factorisation LU suit les étapes suivantes :
  Pour \(k = 1, \ldots, n-1\) :
  - \(L_{k+1, k} = A_{k+1, k} / A_{k, k}\)
  - \(A_{k+1, k+1} = A_{k+1, k+1} - L_{k+1, k} \cdot A_{k, k+1}\)

---

## **2. Méthode de validation**
- Factorisez \(A\) en \(LU\) avec la méthode manuelle.
- Comparez avec la solution fournie par **dgbtrf**.
- Calculer la norme de la différence entre les deux matrices :
  \(\|A - L \cdot U\|\) à l'aide de **dnrm2**.
- Assurez-vous que la norme de l'erreur est proche de zéro.

