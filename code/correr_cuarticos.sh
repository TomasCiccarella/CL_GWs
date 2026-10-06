#!/usr/bin/env bash
# Corre en secuencia los .in cuárticos (o cuadráticos, con FAMILIA=cuadratico) de los modelos presentes
# en esta rama y analiza cada corrida (Notas.tex, «Parámetros cosmológicos para modelos cuárticos»,
# «... cuadráticos» y «Corridas cuárticas»). code/correr_cuadraticos.sh es el atajo para n = 2.
#
# En main no hay modelos: cada rama (monomial, t-model, e-model) agrega los suyos. El script busca los
# casos conocidos cuyo .in exista en models/parameter-files/.
#
# Uso:  code/correr_cuarticos.sh [N] [kIR]
#   N     tamaño de la grilla (por defecto 64)
#   kIR   por defecto (kIR del .in)*64/N (0.7*64/N en los cuárticos): mantiene k_max de N = 64 y agranda la caja
#         (Notas.tex, «Corridas cuárticas», Recomendación)
#
# Variables opcionales:
#   FAMILIA=cuadratico  casos n = 2 (*_cuadratico.in) en vez de los cuárticos
#   ANALISIS_ARGS="..." argumentos para analisis_cosmolattice.py (p. ej. "--TRH 1e10" para n = 2)
#   CL_BUILD=/ruta     raíz donde están build_<modelo>/ y donde se escribe build_corridas/ (por
#                      defecto, este repositorio); sirve para correr desde un worktree de otra rama
#                      sin recompilar
#   DT=0.007           paso temporal (por defecto el del .in, 0.01)
#   MODELOS="..."      subconjunto y orden de los casos (por defecto, los presentes de esta lista, del
#                      más corto al más largo): gentmodel_cuartico genemodel_cuartico_aK5 genemodel_cuartico
#   CONTROL=1          agrega el monomial cuártico de control (genmonomial_cuartico, rama monomial;
#                      Notas.tex, «El monomial de control»). No es un modelo del paper. En la rama
#                      monomial es el único caso, así que se corre siempre.
#   NTHREADS=8         hilos de OpenMP
#   FORZAR=1           volver a correr aunque la corrida ya haya terminado
#   EXTRA_ARGS="..."   argumentos extra para CosmoLattice (p. ej. "tMax=100 baseSeed=1234")
#   SEMILLA=1234       fija baseSeed y agrega _s1234 al directorio de salida (para comparar semillas
#                      sin pisar otras corridas)
#
# Para que siga corriendo al cerrar la terminal:
#   nohup code/correr_cuarticos.sh 128 > corridas_N128.log 2>&1 &
#
# Las salidas quedan en build_corridas/<caso>_N<N>_kIR<kIR>/ (ignorado por git), con el análisis en
# analisis/. Un resumen de tiempos queda en build_corridas/tiempos.txt.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
CL_BUILD=${CL_BUILD:-$ROOT}
N=${1:-64}
NTHREADS=${NTHREADS:-8}
PARS=$ROOT/models/parameter-files
FAMILIA=${FAMILIA:-cuartico}
if [ -z "${MODELOS:-}" ]; then
  MODELOS=""
  for c in gentmodel_$FAMILIA genemodel_${FAMILIA}_aK5 genemodel_$FAMILIA; do
    [ -f "$PARS/$c.in" ] && MODELOS="$MODELOS $c"
  done
  if [ -f "$PARS/genmonomial_$FAMILIA.in" ] && { [ -n "${CONTROL:-}" ] || [ -z "$MODELOS" ]; }; then
    MODELOS="$MODELOS genmonomial_$FAMILIA"
  fi
fi
MODELOS=$(echo $MODELOS)
if [ -z "$MODELOS" ]; then
  echo "No hay .in cuárticos en esta rama. Cambiá a la rama de un modelo (git checkout t-model, e-model o monomial)."
  exit 1
fi
for c in $MODELOS; do
  [ -f "$PARS/$c.in" ] || { echo "Falta $PARS/$c.in (¿estás en la rama de ese modelo?)"; exit 1; }
done
leer_in() {  # leer_in caso clave -> valor en el .in
  sed -n "s/^$2[[:space:]]*=[[:space:]]*\([^[:space:]#]*\).*/\1/p" "$PARS/$1.in" | head -1
}
primero=${MODELOS%% *}
KIR=${2:-$(python3 -c "print(f'{$(leer_in "$primero" kIR) * 64 / $N:.4g}')")}
DT_IN=$(leer_in "$primero" dt)
SALIDA=$CL_BUILD/build_corridas
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
  local m=$1 b=$CL_BUILD/build_$1
  [ -x "$b/$m" ] && return
  echo "Compilando $m en $b ..."
  mkdir -p "$b"
  local args=(-DMODEL="$m" -DCMAKE_BUILD_TYPE=Release)
  local deps=$CL_BUILD/build_analisis/_deps  # reutiliza TempLat y Kokkos ya descargados, si están
  [ -d "$deps/templat-src" ] && args+=(-DFETCHCONTENT_SOURCE_DIR_TEMPLAT="$deps/templat-src")
  [ -d "$deps/kokkos-src" ] && args+=(-DFETCHCONTENT_SOURCE_DIR_KOKKOS="$deps/kokkos-src")
  (cd "$b" && cmake "$ROOT" "${args[@]}" > cmake.log 2>&1 && make -j4 > make.log 2>&1) \
    || { echo "Falló la compilación de $m (ver $b/make.log)"; exit 1; }
}

# --- Chequeo de memoria y tiempo antes de empezar (code/recursos_cosmolattice.py) ---------
python3 -c "
import re, sys; sys.path.insert(0, '$ROOT/code')
import recursos_cosmolattice as r
N, kIR, dt = $N, $KIR, ${DT:-$DT_IN}
def leer(caso, clave):  # valor de una clave en el .in
    txt = open('$PARS/' + caso + '.in', encoding='utf-8').read()
    return float(re.search(r'^' + clave + r'\s*=\s*(\S+)', txt, re.M).group(1))
total = 0
for caso in '$MODELOS'.split():
    tmax = leer(caso, 'tMax')
    x = r.recursos(N, kIR, leer(caso, 'q'), tmax, dt=dt)
    total += x['t_total_h']
    print(f'  {caso:24s} tMax = {tmax:g}: ~{x[\"t_total_h\"]:.1f} h')
    if not x['estable']:
        sys.exit(f'dt = {dt} es inestable para N = {N}, kIR = {kIR} (omega_max dt > 2)')
libre = r.memoria_libre()
print(f'N = {N}, kIR = {kIR}, dt = {dt}: memoria ~{x[\"mem_GB\"]:.2f} GB por corrida '
      f'(libre ahora: {libre/1e9:.2f} GB); total estimado ~{total:.1f} h')
if libre and x['mem_GB'] * 1e9 > 0.8 * libre:
    sys.exit('No entra en la memoria libre: bajá N o cerrá programas (Notas.tex, «Corridas cuárticas», Memoria).')
"

echo "Inicio: $(date)"
for caso in $MODELOS; do
  m=$(ejecutable "$caso")
  compilar "$m"
  dir=$SALIDA/${caso}_N${N}_kIR${KIR}${DT:+_dt$DT}${SEMILLA:+_s$SEMILLA}
  if [ -z "${FORZAR:-}" ] && grep -qs "finished" "$dir/$m.infos"; then
    echo "$caso: ya terminada en $dir (FORZAR=1 para repetir)"
    continue
  fi
  mkdir -p "$dir"
  echo "$caso: corriendo en $dir ($(date +%H:%M))"
  t0=$(date +%s)
  if (cd "$dir" && OMP_NUM_THREADS=$NTHREADS "$CL_BUILD/build_$m/$m" \
        input="$PARS/$caso.in" overwriteFiles=true \
        N="$N" kIR="$KIR" ${DT:+dt=$DT} ${SEMILLA:+baseSeed=$SEMILLA} ${EXTRA_ARGS:-} > salida.log 2>&1); then
    t1=$(date +%s)
    echo "$caso N=$N kIR=$KIR ${DT:+dt=$DT} ${SEMILLA:+seed=$SEMILLA} $(( (t1 - t0) / 60 )) min ($(date))" >> "$SALIDA/tiempos.txt"
    echo "$caso: terminó en $(( (t1 - t0) / 60 )) min; analizando..."
    python3 "$ROOT/code/analisis_cosmolattice.py" "$dir" ${ANALISIS_ARGS:-} | grep -E "Pico|Friedmann|<w>|k_IR" || true
  else
    echo "$caso: CosmoLattice falló (ver $dir/salida.log); sigo con el siguiente."
  fi
done
echo "Fin: $(date)"
