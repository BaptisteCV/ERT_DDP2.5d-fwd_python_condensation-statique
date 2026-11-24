# ert_deps

Dépendances du code `ert_ddp2.5d-fwd_python` (exécutable `ElectricalResistivity`), installées en une commande, sans pagure.

## Installation (Ubuntu)

```bash
sudo apt install build-essential cmake curl xz-utils libblas-dev liblapack-dev

# depuis la racine du depot ert_ddp2.5d-fwd_python
bash ert_deps/install.sh "$PWD"          # ~5 min, ~1 Go dans ert_deps/install

cd src/electrical_resistivity/Release
make clean && make exec -j
```

### Module Python `PERFORM.so` (`make lib`)

Il faut en plus les en-têtes de Python et pybind11. `install.sh` les prend dans le `python3` de la machine, et
prévient s'ils manquent (`make exec` reste possible).

```bash
sudo apt install python3-dev python3-pybind11
bash ert_deps/install.sh "$PWD"          # rapide si déjà installé : rien n'est recompilé
cd src/electrical_resistivity/Release && make clean && make lib -j
```

Attention : `PERFORM.so` ne s'importe qu'avec le `python3` qui a servi à le compiler. Lancé avec un autre Python,
le code passe en mode exécutable sans erreur visible (seulement le message « PERFORM.so not found »).

`install.sh` installe les dépendances dans `ert_deps/install/` (non versionné ; modifiable avec `DEPS_PREFIX=...`), puis écrit
`src/electrical_resistivity/makefile.init` dans le dépôt du code. Ce fichier, propre à la machine, remplace les
chemins par défaut du makefile. Si le dépôt est déplacé, il faut relancer `install.sh` : rien n'est recompilé.

## Contenu

| Dépendance | Version | Provenance | Utilisation |
|---|---|---|---|
| SuiteSparse (UMFPACK, CHOLMOD, AMD, CAMD, COLAMD, CCOLAMD) | 5.10.1 | téléchargée, compilée | solveurs creux, édition de liens |
| METIS | 5.1.0 | fournie avec SuiteSparse | édition de liens |
| Boost | 1.84.0 | téléchargée | en-têtes seulement |
| CGAL | 5.5.2 | téléchargée | en-têtes seulement |
| Eigen | 3.4.0 | `headers/Eigen` | en-têtes seulement |
| boost-numeric-bindings | — | `headers/boost/numeric/bindings` | en-têtes seulement (absent de Boost officiel) |

Ce sont les mêmes versions que sur la machine de référence (Ubuntu 24.04, gcc 13.3). Les résultats ont été vérifiés
identiques, fichier par fichier, à ceux de l'exécutable compilé avec pagure sur le cas `testDFN`.

Licences : Eigen sous MPL 2.0, boost-numeric-bindings sous Boost Software License 1.0.
