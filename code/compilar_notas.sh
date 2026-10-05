#!/usr/bin/env bash
# Compila las notas con el nombre que corresponde a la rama actual:
#   main            -> Notas/Notas.pdf          (notas generales)
#   <rama de modelo> -> Notas/Notas_<rama>.pdf   (notas generales + la sección de ese modelo)
# Los PDFs se versionan; en las ramas de modelo se usa otro nombre para que `git merge main` no
# choque con el Notas.pdf de main.
#
# Uso: code/compilar_notas.sh
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
rama=$(git -C "$ROOT" branch --show-current)
if [ "$rama" = "main" ] || [ ! -d "$ROOT/Notas/modelos" ]; then nombre=Notas; else nombre=Notas_$rama; fi
cd "$ROOT/Notas"
latexmk -g -pdf -interaction=nonstopmode -jobname="$nombre" Notas.tex > /dev/null 2>&1 \
  || { echo "Falló la compilación (ver Notas/$nombre.log)"; exit 1; }
indef=$(grep -c "undefined" "$nombre.log" || true)
echo "Notas/$nombre.pdf: $(pdfinfo "$nombre.pdf" 2>/dev/null | awk '/^Pages/{print $2}') páginas, $indef avisos de referencias indefinidas"
