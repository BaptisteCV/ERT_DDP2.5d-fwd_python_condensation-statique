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

### 1.1 Depuis `main.py` (choisir l'option à exécuter)

Le nom du dossier de scénario (`<SimuType>`, ci-dessus) n'est pas passé en
argument direct au C++ : il transite par un fichier intermédiaire écrit par
Python à chaque forward.

- `main.py:32-33` définit deux noms de dossier, un pour le forward EPM et un
  pour le forward DFN/DDP :
  ```python
  expt_name_epm = 'testEPM_CS'   # tabVar[0]
  expt_name_dfn = 'testDFN'      # tabVar[1]
  ```
- `fwdFunctions.fwdEPM()` / `fwdFunctions.fwdDFN()` (`fwdFunctions.py:32-33`
  et `:51-52`) écrivent respectivement `tabVar[0]` / `tabVar[1]` dans
  `src/input/SimuType.txt` juste avant d'appeler le C++ (`perform()` en mode
  librairie, ou l'exécutable). C'est ce fichier que le C++ lit pour savoir
  dans quel sous-dossier de `src/input/` aller chercher `Simu.txt`,
  `DomainProperties.txt`, `DFN.txt`, etc.
- Le `model_option` effectif est donc **entièrement déterminé par le nom de
  dossier choisi**, puisque chaque dossier existant a un `Simu.txt` déjà figé
  sur une option (ligne 2, cf. §1) et son propre jeu complet de fichiers
  d'entrée (`BoundaryConditions.txt`, `DFN.txt`, `DomainProperties.txt`,
  `ElectrodesPositions.txt`, `FileNames.txt`, `WavesValues.txt`) :

  | Dossier (`src/input/...`) | `model_option` (`Simu.txt` ligne 2) |
  |---|---|
  | `testEPM` | `EPM` (classique, ublas/UMFPACK) |
  | `testEPM_CS` | `EPM_CS` |
  | `testDFN` | `DDP` (classique, ublas/UMFPACK) |
  | `testDFN_CS` | `DDP_CS` |
  | `testDFN_EIGEN` | `DDP_EIGEN` |

**Pour changer d'option, il suffit de changer la valeur de
`expt_name_epm`/`expt_name_dfn` dans `main.py`** pour qu'elle pointe vers le
dossier correspondant — aucune autre modification n'est nécessaire côté
Python (pas besoin de toucher `Simu.txt` ni `fwdFunctions.py`). Par exemple,
pour passer le forward DFN de `DDP` à `DDP_CS` :
```python
expt_name_dfn = 'testDFN_CS'   # au lieu de 'testDFN'
```
ou à `DDP_EIGEN` :
```python
expt_name_dfn = 'testDFN_EIGEN'
```
`expt_name_epm` est déjà réglé sur `'testEPM_CS'` dans le fichier actuel ; il
n'existe pas de dossier `testEPM_EIGEN` (l'option `EIGEN` sans condensation
statique n'est proposée que pour DDP, cf. §4).

⚠️ **Piège dans le mode debug.** En cas d'erreur du C++, `main.py` réactive
le mode debug du binaire (dernière ligne de `Simu.txt` mise à `1`) — mais en
rouvrant **en dur** `input/testEPM/Simu.txt` (ligne 55) et
`input/testDFN/Simu.txt` (ligne 82), et non
`input/<expt_name_epm>/Simu.txt` / `input/<expt_name_dfn>/Simu.txt`. Si
`expt_name_epm`/`expt_name_dfn` a été changé vers `testEPM_CS`, `testDFN_CS`
ou `testDFN_EIGEN`, le mode debug est activé sur le mauvais dossier et celui
réellement utilisé n'est jamais basculé en debug après une erreur. À corriger
si ce cas se présente, en remplaçant le nom en dur par `tabVar[0]` /
`tabVar[1]`.

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


Compilation : `make exec` (exécutable) ou `make lib` (librairie) depuis
`src/electrical_resistivity/Release`. Aucune cible dédiée par option — tous
les fichiers sont compilés systématiquement.
