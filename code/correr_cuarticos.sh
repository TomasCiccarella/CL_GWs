#!/usr/bin/env bash
# Corre en secuencia los tres .in cuárticos (Notas.tex Secs. 4 y 6) y analiza cada corrida.
#
# Uso:  code/correr_cuarticos.sh [N] [kIR]
#   N     tamaño de la grilla (por defecto 64)
#   kIR   por defecto 0.7*64/N: mantiene k_max de N = 64 y agranda la caja (Notas.tex Sec. 6.5)
#
# Variables opcionales:
#   DT=0.007           paso temporal (por defecto el del .in, 0.01)
#   CONTROL=1          agrega al final el monomial cuártico de control (genmonomial_cuartico, Notas.tex
#                      Sec. 7.2): mismo lambda_eff y q que el T-model, sin plateau. No es un modelo del paper.
#   MODELOS="..."      subconjunto y orden de los casos (por defecto los tres del paper, del más corto al
#                      más largo): gentmodel_cuartico genemodel_cuartico_aK5 genemodel_cuartico
#                      (y genmonomial_cuartico si CONTROL=1)
#   NTHREADS=8         hilos de OpenMP
#   FORZAR=1           volver a correr aunque la corrida ya haya terminado
#   EXTRA_ARGS="..."   argumentos extra para CosmoLattice (p. ej. "tMax=100")
#
# Para que siga corriendo al cerrar la terminal:
#   nohup code/correr_cuarticos.sh 128 > corridas_N128.log 2>&1 &
#
# Las salidas quedan en build_corridas/<caso>_N<N>_kIR<kIR>/ (ignorado por git), con el análisis en
# analisis/. Un resumen de tiempos queda en build_corridas/tiempos.txt.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
N=${1:-64}
KIR=${2:-$(python3 -c "print(f'{0.7 * 64 / $N:.4g}')")}
NTHREADS=${NTHREADS:-8}
MODELOS=${MODELOS:-"gentmodel_cuartico genemodel_cuartico_aK5 genemodel_cuartico"}
[ -n "${CONTROL:-}" ] && [[ " $MODELOS " != *" genmonomial_cuartico "* ]] && MODELOS="$MODELOS genmonomial_cuartico"
SALIDA=$ROOT/build_corridas
mkdir -p "$SALIDA"

ejecutable() {  # caso -> modelo de CosmoLattice
  case $1 in
    gentmodel_*) echo gentmodel ;;
    genemodel_*) echo genemodel ;;
    genmonomial_*) echo genmonomial ;;
    *) echo "caso desconocido: $1" >&2; exit 1 ;;
  esac
}

compilar() {  # compila el modelo si no está el ejecutable
  local m=$1 b=$ROOT/build_$1
  [ -x "$b/$m" ] && return
  echo "Compilando $m en $b ..."
  mkdir -p "$b"
  local args=(-DMODEL="$m" -DCMAKE_BUILD_TYPE=Release)
  local deps=$ROOT/build_analisis/_deps  # reutiliza TempLat y Kokkos ya descargados, si están
  [ -d "$deps/templat-src" ] && args+=(-DFETCHCONTENT_SOURCE_DIR_TEMPLAT="$deps/templat-src")
  [ -d "$deps/kokkos-src" ] && args+=(-DFETCHCONTENT_SOURCE_DIR_KOKKOS="$deps/kokkos-src")
  (cd "$b" && cmake "$ROOT" "${args[@]}" > cmake.log 2>&1 && make -j4 > make.log 2>&1) \
    || { echo "Falló la compilación de $m (ver $b/make.log)"; exit 1; }
}

# --- Chequeo de memoria y tiempo antes de empezar (code/recursos_cosmolattice.py) ---------
python3 -c "
import sys; sys.path.insert(0, '$ROOT/code')
import recursos_cosmolattice as r
N, kIR, dt = $N, $KIR, ${DT:-0.01}
tmax = {'gentmodel_cuartico': 500, 'genemodel_cuartico_aK5': 500, 'genemodel_cuartico': 600,
        'genmonomial_cuartico': 500}
total = 0
for caso in '$MODELOS'.split():
    x = r.recursos(N, kIR, 120, tmax[caso], dt=dt)
    total += x['t_total_h']
    print(f'  {caso:24s} tMax = {tmax[caso]}: ~{x[\"t_total_h\"]:.1f} h')
    if not x['estable']:
        sys.exit(f'dt = {dt} es inestable para N = {N}, kIR = {kIR} (omega_max dt > 2)')
libre = r.memoria_libre()
print(f'N = {N}, kIR = {kIR}, dt = {dt}: memoria ~{x[\"mem_GB\"]:.2f} GB por corrida '
      f'(libre ahora: {libre/1e9:.2f} GB); total estimado ~{total:.1f} h')
if libre and x['mem_GB'] * 1e9 > 0.8 * libre:
    sys.exit('No entra en la memoria libre: bajá N o cerrá programas (Notas.tex Sec. 6.2).')
"

echo "Inicio: $(date)"
for caso in $MODELOS; do
  m=$(ejecutable "$caso")
  compilar "$m"
  dir=$SALIDA/${caso}_N${N}_kIR${KIR}${DT:+_dt$DT}
  if [ -z "${FORZAR:-}" ] && grep -qs "finished" "$dir/$m.infos"; then
    echo "$caso: ya terminada en $dir (FORZAR=1 para repetir)"
    continue
  fi
  mkdir -p "$dir"
  echo "$caso: corriendo en $dir ($(date +%H:%M))"
  t0=$(date +%s)
  if (cd "$dir" && OMP_NUM_THREADS=$NTHREADS "$ROOT/build_$m/$m" \
        input="$ROOT/models/parameter-files/$caso.in" overwriteFiles=true \
        N="$N" kIR="$KIR" ${DT:+dt=$DT} ${EXTRA_ARGS:-} > salida.log 2>&1); then
    t1=$(date +%s)
    echo "$caso N=$N kIR=$KIR ${DT:+dt=$DT} $(( (t1 - t0) / 60 )) min ($(date))" >> "$SALIDA/tiempos.txt"
    echo "$caso: terminó en $(( (t1 - t0) / 60 )) min; analizando..."
    python3 "$ROOT/code/analisis_cosmolattice.py" "$dir" | grep -E "Pico|Friedmann|<w>|k_IR" || true
  else
    echo "$caso: CosmoLattice falló (ver $dir/salida.log); sigo con el siguiente."
  fi
done
echo "Fin: $(date)"
