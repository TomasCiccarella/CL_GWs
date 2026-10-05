#!/usr/bin/env bash
# Corre y analiza los casos cuadráticos (n = 2) de esta rama, incluido el monomial m^2 phi^2 de
# control (Notas.tex, «Parámetros cosmológicos para modelos cuadráticos»). Es correr_cuarticos.sh con
# FAMILIA=cuadratico, CONTROL=1 y el reheating de referencia T_RH = 1e10 GeV para el reescaleo a hoy.
#
# Uso:  code/correr_cuadraticos.sh [N] [kIR]     (mismas variables opcionales que correr_cuarticos.sh;
#                                                  TRH=... cambia la temperatura de reheating)
exec env FAMILIA=cuadratico CONTROL=1 ANALISIS_ARGS="--TRH ${TRH:-1e10} ${ANALISIS_ARGS:-}" \
  "$(dirname "$0")/correr_cuarticos.sh" "$@"
