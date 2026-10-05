#!/usr/bin/env bash
# Compila genmonomial, gentmodel y genemodel y corre cada uno con su .in de models/parameter-files/,
# para comprobar que la simulación termina sin NaN.
#
# Uso: code/probar_genmodels.sh [directorio_de_build]   (por defecto ./build_genmodels)
#
# Variables opcionales:
#   TEMPLAT_SRC, KOKKOS_SRC  copias locales de TempLat v1.0.2 y Kokkos, para compilar sin red
#   NTHREADS                 hilos de OpenMP (por defecto 8)
#   EXTRA_ARGS               argumentos extra para todas las corridas (p. ej. "N=16 tMax=5")
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BUILD=${1:-$ROOT/build_genmodels}
mkdir -p "$BUILD"
cd "$BUILD"

CMAKE_ARGS=(-DCMAKE_BUILD_TYPE=Release)
[ -n "${TEMPLAT_SRC:-}" ] && CMAKE_ARGS+=(-DFETCHCONTENT_SOURCE_DIR_TEMPLAT="$TEMPLAT_SRC")
[ -n "${KOKKOS_SRC:-}" ] && CMAKE_ARGS+=(-DFETCHCONTENT_SOURCE_DIR_KOKKOS="$KOKKOS_SRC")

estado=0
for m in genmonomial gentmodel genemodel; do
  cmake "$ROOT" -DMODEL=$m "${CMAKE_ARGS[@]}" > cmake_$m.log 2>&1
  make -j"${NTHREADS:-8}" > make_$m.log 2>&1
  mkdir -p run_$m
  (cd run_$m && OMP_NUM_THREADS=${NTHREADS:-8} ../$m input="$ROOT/models/parameter-files/$m.in" \
     overwriteFiles=true ${EXTRA_ARGS:-} > salida.log 2>&1)
  ultima=$(tail -1 run_$m/average_scale_factor.txt)
  if grep -qi nan run_$m/average_scalar_0.txt; then
    echo "$m: NaN (primera fila: $(grep -i -m1 nan run_$m/average_scalar_0.txt | cut -f1))"
    estado=1
  else
    echo "$m: sin NaN; última fila t a a' H = $ultima"
  fi
done
exit $estado
