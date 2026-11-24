#!/bin/bash
# Installation des dependances du code ert_ddp2.5d-fwd_python (executable ElectricalResistivity), sans pagure.
#
#   bash install.sh [CHEMIN_DU_CODE]
#
# CHEMIN_DU_CODE : racine du depot ert_ddp2.5d-fwd_python. Si donne, le script y ecrit
# src/electrical_resistivity/makefile.init, qui donne au makefile les chemins des dependances installees.
#
# Installe dans DEPS_PREFIX (par defaut : ./install a cote de ce script) :
#   - SuiteSparse 5.10.1 (UMFPACK, CHOLMOD, AMD, CAMD, COLAMD, CCOLAMD) et METIS 5.1.0 : compiles ;
#   - Boost 1.84.0 et CGAL 5.5.2 : en-tetes seulement, rien a compiler.
# Eigen 3.4.0 et boost-numeric-bindings sont fournis dans headers/ (introuvables sous cette forme ailleurs).
# Memes versions que sur la machine de reference (Ubuntu 24.04, gcc 13.3).
#
# Prerequis (Ubuntu) : sudo apt install build-essential cmake curl xz-utils libblas-dev liblapack-dev
# Pour make lib (PERFORM.so) en plus : en-tetes de python3 et pybind11 (sudo apt install python3-dev python3-pybind11)
# Variables reglables : DEPS_PREFIX, JOBS, BLAS_LIBS, LAPACK_LIBS.
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CODE_DIR=${1:-}
DEPS_PREFIX=${DEPS_PREFIX:-$HERE/install}
JOBS=${JOBS:-$(nproc)}
BLAS_LIBS=${BLAS_LIBS:--lblas}
LAPACK_LIBS=${LAPACK_LIBS:--llapack}

SS_VERSION=5.10.1
BOOST_VERSION=1.84.0
CGAL_VERSION=5.5.2
SS_PREFIX=$DEPS_PREFIX/suitesparse-$SS_VERSION
BOOST_PREFIX=$DEPS_PREFIX/boost-$BOOST_VERSION
CGAL_PREFIX=$DEPS_PREFIX/CGAL-$CGAL_VERSION

if [ -n "$CODE_DIR" ] && [ ! -d "$CODE_DIR/src/electrical_resistivity/Release" ]; then
    echo "ERREUR : $CODE_DIR n'est pas la racine du depot ert_ddp2.5d-fwd_python" >&2
    exit 1
fi

# ---------- Verification des prerequis ----------
missing=""
for cmd in g++ make cmake curl tar xz; do
    command -v $cmd > /dev/null || missing="$missing $cmd"
done
if ! echo 'int main(){}' | g++ -x c++ - -o /dev/null $LAPACK_LIBS $BLAS_LIBS 2> /dev/null; then
    missing="$missing BLAS/LAPACK($BLAS_LIBS $LAPACK_LIBS)"
fi
if [ -n "$missing" ]; then
    echo "ERREUR : prerequis manquants :$missing" >&2
    echo "  sudo apt install build-essential cmake curl xz-utils libblas-dev liblapack-dev" >&2
    exit 1
fi

mkdir -p "$DEPS_PREFIX/src"
cd "$DEPS_PREFIX/src"

# ---------- Boost (en-tetes seulement) ----------
if [ ! -f "$BOOST_PREFIX/include/boost/version.hpp" ]; then
    echo "=== Boost $BOOST_VERSION (en-tetes)"
    B=boost_${BOOST_VERSION//./_}
    [ -f $B.tar.gz ] || curl -fL -o $B.tar.gz https://archives.boost.io/release/$BOOST_VERSION/source/$B.tar.gz
    tar xzf $B.tar.gz $B/boost
    mkdir -p "$BOOST_PREFIX/include"
    rm -rf "$BOOST_PREFIX/include/boost"
    mv $B/boost "$BOOST_PREFIX/include/"
    rm -rf $B
fi

# ---------- CGAL (en-tetes seulement) ----------
if [ ! -d "$CGAL_PREFIX/include/CGAL" ]; then
    echo "=== CGAL $CGAL_VERSION (en-tetes)"
    [ -f CGAL-$CGAL_VERSION.tar.xz ] || curl -fL -o CGAL-$CGAL_VERSION.tar.xz \
        https://github.com/CGAL/cgal/releases/download/v$CGAL_VERSION/CGAL-$CGAL_VERSION.tar.xz
    tar xJf CGAL-$CGAL_VERSION.tar.xz -C "$DEPS_PREFIX"
fi

# ---------- SuiteSparse + METIS (compiles) ----------
if [ ! -e "$SS_PREFIX/lib/libumfpack.so" ]; then
    echo "=== SuiteSparse $SS_VERSION + METIS"
    if [ ! -d SuiteSparse-$SS_VERSION ]; then
        [ -f SuiteSparse-$SS_VERSION.tar.gz ] || curl -fL -o SuiteSparse-$SS_VERSION.tar.gz \
            https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v$SS_VERSION.tar.gz
        tar xzf SuiteSparse-$SS_VERSION.tar.gz
    fi
    # METIS 5.1.0 demande cmake_minimum_required 2.8, refuse par CMake >= 4 sans cette variable
    export CMAKE_POLICY_VERSION_MINIMUM=3.5
    cd SuiteSparse-$SS_VERSION
    SS_OPTS=(JOBS=$JOBS CUDA=no "BLAS=$BLAS_LIBS" "LAPACK=$LAPACK_LIBS" INSTALL=$SS_PREFIX)

    # Le Makefile de SuiteSparse ignore les erreurs de METIS, d'ou la verification
    make metis "${SS_OPTS[@]}"
    if [ ! -e lib/libmetis.so ]; then
        echo "ERREUR : METIS n'a pas ete compile (voir les messages de cmake ci-dessus)" >&2
        exit 1
    fi
    # Seuls les paquets utilises par le code (le "make library" complet compile aussi GraphBLAS, tres long)
    PACKAGES="SuiteSparse_config AMD CAMD COLAMD CCOLAMD CHOLMOD UMFPACK"
    for p in $PACKAGES; do
        make -C $p library "${SS_OPTS[@]}"
    done
    make metisinstall "${SS_OPTS[@]}"
    for p in $PACKAGES; do
        make -C $p install "${SS_OPTS[@]}"
    done
    cd ..
fi

# ---------- Python + pybind11 (uniquement pour make lib, PERFORM.so) ----------
# En-tetes pris dans le python3 de la machine : PERFORM.so ne s'importera qu'avec ce meme python3.
PY_INCLUDE=$(python3 -c 'import sysconfig; print(sysconfig.get_paths()["include"])' 2> /dev/null || true)
PYBIND_INCLUDE=$(python3 -c 'import pybind11; print(pybind11.get_include())' 2> /dev/null || true)
LIB_OK=1
if [ -z "$PY_INCLUDE" ] || [ ! -f "$PY_INCLUDE/Python.h" ]; then
    LIB_OK=0
    echo "ATTENTION : Python.h introuvable pour python3 (sudo apt install python3-dev)." >&2
fi
if [ -z "$PYBIND_INCLUDE" ]; then
    LIB_OK=0
    echo "ATTENTION : pybind11 introuvable pour python3 (sudo apt install python3-pybind11)." >&2
fi
if [ $LIB_OK = 0 ]; then
    echo "  -> make exec fonctionnera, pas make lib (PERFORM.so). Relancer install.sh une fois corrige." >&2
    PY_INCLUDE=${PY_INCLUDE:-$BOOST_PREFIX/include}
    PYBIND_INCLUDE=${PYBIND_INCLUDE:-$BOOST_PREFIX/include}
fi

# ---------- makefile.init ----------
INIT=$DEPS_PREFIX/makefile.init
cat > "$INIT" <<EOF
# Genere par $HERE/install.sh : chemins des dependances pour src/electrical_resistivity/Release/makefile.
# Ce fichier est propre a la machine (non versionne) ; le placer dans src/electrical_resistivity/makefile.init.
DEPS_PREFIX := $DEPS_PREFIX
PAGURE_ROOT := \$(DEPS_PREFIX)
SUITESPARSE_ROOT := \$(DEPS_PREFIX)/suitesparse-$SS_VERSION
METIS_ROOT := \$(SUITESPARSE_ROOT)
CGAL_ROOT := \$(DEPS_PREFIX)/CGAL-$CGAL_VERSION
BOOST_ROOT := \$(DEPS_PREFIX)/boost-$BOOST_VERSION
# Eigen 3.4.0 + boost/numeric/bindings
EIGEN_INCLUDE := $HERE/headers
EXTRA_INCLUDE_FLAGS :=
# Utilises seulement par make lib (PERFORM.so) : en-tetes de pybind11 et de Python ($(python3 --version 2>&1))
PYMOD_INCLUDE := $PYBIND_INCLUDE
PYTHON_INCLUDE := $PY_INCLUDE
EXTRA_CXXFLAGS :=
EOF

echo
echo "Dependances installees dans $DEPS_PREFIX"
if [ -n "$CODE_DIR" ]; then
    cp "$INIT" "$CODE_DIR/src/electrical_resistivity/makefile.init"
    echo "makefile.init ecrit dans $CODE_DIR/src/electrical_resistivity/"
    echo "Compilation : cd $CODE_DIR/src/electrical_resistivity/Release && make clean && make exec -j$JOBS"
else
    echo "Copier $INIT dans <depot du code>/src/electrical_resistivity/makefile.init"
fi
