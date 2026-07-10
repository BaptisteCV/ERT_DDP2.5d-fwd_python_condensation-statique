# Documentation technique — options `EPM_CS`, `DDP_CS`, `DDP_EIGEN`

Ce document décrit les trois variantes Eigen ajoutées au code historique
(`EPM`/`DDP` en `ublas`, résolus via UMFPACK) : comment les appeler, ce
qu'elles calculent, et la chaîne exacte de fonctions/fichiers impliqués.

| Option | Principe | Conteneurs |
|---|---|---|
| `EPM_CS` | Milieu poreux équivalent (homogène), avec **condensation statique** (Schur) | Eigen |
| `DDP_CS` | Double porosité (matrice + fractures), avec **condensation statique** | Eigen |
| `DDP_EIGEN` | Double porosité, résolution directe (comme `DDP` classique) mais en conteneurs Eigen natifs | Eigen |

La condensation statique (`_CS`) réduit le système linéaire aux seuls nœuds
de bord ; le bloc intérieur (indépendant de la position de la source) n'est
factorisé **qu'une fois par fréquence** et réutilisé pour toutes les
électrodes — c'est le principal gain de performance par rapport à `DDP`/`EPM`
classiques, qui refont une résolution complète à chaque électrode.

---

## 1. Comment les appeler

Mot-clé `model_option` dans `Simu.txt` du dossier de scénario
(`src/input/<SimuType>/Simu.txt`, ligne 2) : `EPM_CS`, `DDP_CS` ou
`DDP_EIGEN`. Le reste du fichier (`Nx`, `Ny`, `debug_mode`) est identique aux
autres options. Constantes définies dans
`Utilitaries/Constantes.h:79-80` (`DDP_CS`, `DDP_EIGEN`) et l'équivalent pour
`EPM_CS`.

Deux points d'entrée C++, dispatchés sur `model_option` :
- **`MAIN.cpp`** (mode exécutable, binaire `ElectricalResistivity`) : lit
  `src/input/<SimuType>/`, écrit les résultats dans `src/output/<SimuType>/`.
- **`PERFORM.cpp`** (mode librairie, `PERFORM.so` via pybind11, fonction
  `perform()`) : mêmes fichiers d'entrée, mais renvoie directement
  `dict[int, Eigen::SparseMatrix]` à Python (aucun fichier intermédiaire).

Les deux appellent les **mêmes fonctions C++** ; seul un paramètre
`build_result_map` diffère (voir §5).

---

## 2. `EPM_CS`

- **Fichier** : `FourierModel/FourierModelFunctionsCSEigen.cpp:116`
  (namespace `eigen_cs`)
- **Déclaration** : `FourierModel/FourierModelFunctionsCSEigen.h:13-19`
- **Fonction** : `FourierPotentialEPMSourceTermCSEigen(param, domain, singularity, bc_on_border, source_term_positions, build_result_map=false)`
- **Classe** : `EPMSystemCSEigen` (`EPMModel/EPMDefinitionCSEigen.h/.cpp`)

Chaîne d'appels :
1. Construction de la conductivité (`DomainProperties.txt` ou distribution hétérogène)
2. `EPMSystemCSEigen epd_model_base(...)` puis :
   - `InitializeSCIndexing()` — `EPMModel/EPMDefinitionCSEigen.cpp:90`
   - `DetermineBasicLinearSystemCS()` — `EPMDefinitionCSEigen.cpp:115`
3. Boucle fréquences (`WavesValues.txt`) :
   - `UpdateSCBlocksFourier()` — `EPMDefinitionCSEigen.cpp:187`
   - `PrecomputeSCForCurrentFrequency()` — `EPMDefinitionCSEigen.cpp:291` (factorisation du bloc intérieur, 1×/fréquence)
   - Boucle électrodes :
     - `BoundaryConditionsDefinitionCSEigen(...)` — `EPMModel/EPMComputationCSEigen.cpp:141`
     - `UpdateSCBoundaryConditionsAndSources(fourier)` — `EPMDefinitionCSEigen.cpp:212`
     - `SolveCondensedSystem()` — `EPMDefinitionCSEigen.cpp:431`
     - `ReconstructInteriorFromBorder()` — `EPMDefinitionCSEigen.cpp:439`
     - `BuildGlobalSolutionVectorFromSC()` — `EPMDefinitionCSEigen.cpp:456`
4. Si `build_result_map==false` (mode exécutable) : écriture disque `Results<n>.txt` (préfixe attendu par `post_processing/mainEPM.py:61`)
5. Si `build_result_map==true` (mode librairie) : empilement des potentiels (toutes fréquences) en une matrice creuse Eigen par électrode → `res_matrices[elec_nb]`, **sans écriture disque**

### Fichiers lus
`DomainProperties.txt` (+ `porous_cond.txt` si hétérogène), `BoundaryConditions.txt`, `WavesValues.txt`, `ElectrodesPositions.txt`.

---

## 3. `DDP_CS`

- **Fichier** : `FourierModel/FourierModelFunctionsCSEigen.cpp:337`
  (namespace `eigen_cs`)
- **Déclaration** : `FourierModelFunctionsCSEigen.h:21-26`
- **Fonction** : `FourierPotentialDDPSourceTermCSEigen(param, domain, singularity, bc_on_border, source_term_positions, build_result_map=false)`
- **Classe** : `DDPModelCS` (`DDPModel/DDPModelCS.h/.cpp`), hérite de `EPMSystemCSEigen` et ajoute les blocs fracture (`Afi_*`, `Afm_*`, `Amf_*`)

Chaîne d'appels (même principe que `EPM_CS`, + réseau de fractures) :
1. Conductivité + `NetworkMeshes net_mesh(domain, param, simu_option)` (`DFN.txt`)
2. `DDPModelCS ddp_base(...)` puis :
   - `InitializeDFNSCIndexing()` — `DDPModel/DDPModelCS.cpp:155`
   - `DetermineBasicLinearSystemDDPCS()` — `DDPModelCS.cpp:295`
3. Boucle fréquences :
   - `UpdateDDPSCBlocksFourier()` — `DDPModelCS.cpp:497`
   - `PrecomputeDDPSCForCurrentFrequency()` — `DDPModelCS.cpp:853` (Schur + LLT du bloc intérieur, 1×/fréquence)
   - Boucle électrodes (parallèle OpenMP, `DDP_CS_NB_THREADS=14`) :
     - `BoundaryConditionsDefinitionCSEigen(...)` — `EPMModel/EPMComputationCSEigen.cpp:141` (partagée avec `EPM_CS`)
     - `UpdateDDPSCBoundaryConditionsAndSources(fourier)` — `DDPModelCS.cpp:670`
     - `SolveCondensedDDPSystem()` — `DDPModelCS.cpp:1363`
     - `ReconstructDDPInteriorFromBorder()` — `DDPModelCS.cpp:1384`
     - `BuildGlobalSolutionVectorFromSC()` — `EPMModel/EPMDefinitionCSEigen.cpp:456` (héritée)
4. Si `build_result_map==false` (mode exécutable) : écriture disque `Matrix_Potential<n>.txt` (préfixe attendu par `mainDFN.py:46` — renommé depuis `ResultsDDP` pour compatibilité)
5. Si `build_result_map==true` (mode librairie) : même empilement que `EPM_CS` → `res_matrices[elec_nb]`, **sans écriture disque**

### Fichiers lus
`DFN.txt`, `DomainProperties.txt` (+ `porous_cond.txt`), `BoundaryConditions.txt`, `WavesValues.txt`, `ElectrodesPositions.txt`.

---

## 4. `DDP_EIGEN`

- **Fichier** : `FourierModel/FourierModelFunctions.cpp:544`
- **Déclaration** : `FourierModelFunctions.h:24`
- **Fonction** : `FourierPotentialDDPSourceTermEigen(param, domain, singularity, source_term_positions, build_result_map=false)`

Contrairement aux deux options `_CS`, pas de condensation statique : chaque
électrode/fréquence résout le système complet, comme `DDP` classique, mais
avec des conteneurs Eigen natifs de bout en bout (pas de `ublas`).

Chaîne d'appels :
1. `NetworkMeshes net_mesh(domain, param, simu_option)` (`DFN.txt`)
2. `DDPModel ddp_model_init(...)` + `DetermineBasicLinearSystemDDP()` — `DDPModel/DDPModel.cpp:378` (système de base construit en `ublas`, converti ensuite)
3. `DDPEigenBase base = DetermineBasicLinearSystemDDPEigen(ddp_model_init)` — `DDPModel/SystLinearFunctionsEigen.cpp:25` (conversion `ublas→Eigen` du stencil de base, une seule fois)
4. `solver.analyzePattern(...)` — un `DDPEigenSolver` par thread OpenMP (alias `Eigen::SimplicialLLT` ou `Eigen::UmfPackLU` selon macro, `SystLinearFunctionsEigen.h:18/21`)
5. Boucle électrodes (parallèle OpenMP si `DDP_EIGEN_PARALLEL_SOURCES`) × fréquences :
   - `DDPModelFourierBC(...)` — `DDPModel/DDPModel.cpp:620`
   - `UpdateLinearSystemFourierEigen(...)` — `SystLinearFunctionsEigen.cpp:232`
   - `EvaluatePotentialEigen(...)` — `SystLinearFunctionsEigen.cpp:266`
6. Si `build_result_map==false` (mode exécutable) : écriture disque `Matrix_Potential<n>.txt` (réactivée — était commentée)
7. Si `build_result_map==true` (mode librairie) : même empilement que `DDP_CS`/`EPM_CS` → `res_matrices[nb_elec]`, **sans écriture disque**

### Fichiers lus
`DFN.txt`, `DomainProperties.txt`, `BoundaryConditions.txt` (+ `BoundaryConditionsSing.txt`), `WavesValues.txt`, `ElectrodesPositions.txt`.

---

## 5. Le paramètre `build_result_map`

Les trois fonctions ci-dessus (+ `DDP`/`EPM` classiques) sont **symétriques**
sur ce paramètre — il ne contrôle pas seulement l'assemblage de la map de
retour, mais aussi (et surtout) l'écriture disque, pour éviter de payer les
deux à la fois :

| | `build_result_map=false` (défaut, `MAIN.cpp`) | `build_result_map=true` (`PERFORM.cpp`) |
|---|---|---|
| Écriture disque (`Matrix_Potential<n>.txt` / `Results<n>.txt`) | ✅ oui | ❌ non |
| Assemblage de `res_matrices` (retour Python) | ❌ non (map vide) | ✅ oui |

- `MAIN.cpp` ne passe pas ce paramètre (défaut `false`) : les résultats sont
  écrits sur disque, personne ne lit la valeur de retour → inutile de payer
  l'assemblage en mémoire.
- `PERFORM.cpp` passe toujours `true` : Python reçoit le potentiel
  directement via `res_map` (conteneur Eigen, sans passer par le disque),
  donc l'écriture fichier serait un coût I/O pur — d'où sa suppression dans
  ce mode. C'est précisément l'intérêt du mode librairie par rapport au mode
  exécutable : éviter l'aller-retour disque.

Pour `EPM_CS`/`DDP_CS`/`DDP_EIGEN`, la map assemblée est déjà en
`Eigen::SparseMatrix` : `PERFORM.cpp` l'assigne directement à `res_map`, sans
conversion. Seuls `DDP`/`EPM` classiques (fondamentalement en `ublas`)
passent encore par `UblasMapToEigenMap(...)` pour la conversion finale vers
Python — c'est une conversion **une fois, à la fin**.

---

## 6. Solveurs internes

| Option | Solveur | Fichier |
|---|---|---|
| `DDP`/`EPM` classiques | **UMFPACK** (LU direct, SuiteSparse) | `Utilitaries/Math/LinearSystem.cpp:13` (macro `USE_UMFPACK`) |
| `DDP_EIGEN` | `Eigen::SimplicialLLT` (défaut) ou `Eigen::UmfPackLU` | `DDPModel/SystLinearFunctionsEigen.h:18/21` |
| `DDP_CS`/`EPM_CS` | Factorisation LLT du bloc intérieur (Schur), résolution dense pour le système condensé | `DDPModelCS.cpp` / `EPMDefinitionCSEigen.cpp` |

`LinearSystemSolving` (`DDP`/`EPM` classiques) construit le système en
`ublas`, puis extrait un format CSC générique pour appeler UMFPACK
directement — aucune conversion Eigen n'intervient dans la résolution
elle-même (le solveur Eigen (`USE_EIGEN`) reste disponible en alternative
dans ce fichier, mais désactivé par défaut).

---

## 7. Dépendances communes

Types centralisés dans `Utilitaries/Storage/EigenCSStructures.h`
(`eigen_cs::SparseMatrix`, `DenseMatrix`, `Vector`). Parallélisation OpenMP
sur la boucle électrodes. `PERFORM.so` (mode librairie) via pybind11
(`pybind11/eigen.h`) — nécessite un interpréteur Python de la même version
ABI que celle utilisée à la compilation.

Compilation : `make exec` (exécutable) ou `make lib` (librairie) depuis
`src/electrical_resistivity/Release`. Aucune cible dédiée par option — tous
les fichiers sont compilés systématiquement.
