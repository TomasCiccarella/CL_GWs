# Ondas gravitacionales del preheating en modelos de inflación favorecidos por las observaciones

Este repositorio estudia el **fondo estocástico de ondas gravitacionales (GWs) que se produce al final de la inflación**, durante el *preheating*, mediante simulaciones en la red con [CosmoLattice](https://cosmolattice.net). El objetivo es obtener el espectro de GWs **reescalado a hoy**, es decir la amplitud $h^2\Omega_{GW}(f)$ en función de la frecuencia $f$, para modelos de inflación que siguen siendo compatibles con las observaciones del CMB y para sus generalizaciones.

## Índice

1. [La idea](#la-idea): qué se estudia y por qué.
2. [Cómo leer el repositorio](#cómo-leer-el-repositorio): por dónde empezar según lo que quieras hacer.
3. [Organización en ramas](#organización-en-ramas): qué hay en `main` y qué agrega cada rama de modelo.
4. [Guía de carpetas](#guía-de-carpetas): qué tiene cada carpeta y cada archivo.
5. [Cómo usarlo](#cómo-usarlo): comandos para compilar, correr, analizar y compilar las notas.
6. [Créditos](#créditos)

## La idea

Después de la inflación, el inflatón $\phi$ oscila alrededor del mínimo de su potencial. Si está acoplado a otro campo $\chi$ (acá, con $\tfrac12 g^2\phi^2\chi^2$), las oscilaciones amplifican exponencialmente las fluctuaciones de $\chi$ (**resonancia paramétrica**). Las inhomogeneidades resultantes son una fuente de ondas gravitacionales. Ese fondo de GWs queda como una huella de la física del final de la inflación, a frecuencias altas (entre $10^8$ y $10^9$ Hz).

Estudiamos potenciales que cerca del mínimo se comportan como $|\phi|^n$:

| modelo | potencial | comentario |
|---|---|---|
| monomial generalizado | $V = A\lvert\phi\rvert^n$ | caso de referencia; con $n = 4$ es el `lphi4` de CosmoLattice |
| T-model | $V = A\tanh^n(\phi/M)$ | $\alpha$-attractor con plateau |
| E-model | $V = A(1 - e^{-\phi/M})^n$ | $\alpha$-attractor ($\alpha$-Starobinsky generalizado) |

**Modelos favorecidos por las observaciones.** Los parámetros de referencia salen de Ellis, Garcia, Olive & Verner, *Phys. Rev. D* **113**, 063571 (2026) ([doi:10.1103/d35r-7bn8](https://doi.org/10.1103/d35r-7bn8)). Ese trabajo contrasta los T- y E-models generalizados con Planck, BICEP/Keck, ACT DR6 y SPT-3G. Nos concentramos en el **caso cuártico ($n = 4$)** por dos razones:

- **Es compatible con el CMB:** para $\alpha_K = 1$ da $n_s \simeq 0.965$ y $r \simeq 0.004$.
- **El resultado no depende del reheating:** el inflatón oscilando se comporta como radiación ($w = 1/3$). Por eso tanto el número de e-folds $N_*$ como el espectro de GWs de hoy no dependen de la temperatura de reheating, que es desconocida.

Los casos que se simulan son el T-model con $\alpha_K = 1$ y el E-model con $\alpha_K = 1$ y $5$. Se agrega un monomial $\phi^4$ de control (excluido por $r$, pero con el mismo mínimo que el T-model) para aislar el efecto del plateau.

**Del lattice a hoy.** El espectro que da CosmoLattice se lleva a hoy con

```math
f_0 = \frac{\tilde k}{a_e\,\tilde\rho_e^{1/4}}\left(\frac{\omega_*}{f_*}\right)^{1/2}\epsilon^{1/4}\times 4.59\times10^{10}\,\text{Hz},
\qquad
h^2\Omega_{GW,0} = 1.61\times10^{-5}\;\epsilon\;\Omega_{GW,e},
```

con $\epsilon = (a_e/a_{RD})^{2(4-n)/(n+2)}$, que vale exactamente $1$ para $n = 4$. La derivación está verificada con sympy y comparada con Dufaux *et al.* (2007) y Figueroa & Torrentí (2017).

## Cómo leer el repositorio

- **Para entender la física:** leé [La idea](#la-idea) y después las notas (`Notas/Notas.tex`; el PDF se compila con `latexmk`). Las notas están en este orden:
  1. el marco común para potenciales con mínimo $|\phi|^n$: variables de programa, reescaleo temporal, condiciones iniciales y estabilidad;
  2. los parámetros del paper para $n = 4$;
  3. el análisis y el reescaleo a hoy;
  4. los requisitos de cómputo;
  5. la guía de uso, con cada parámetro de los `.in`;
  6. al final, la sección del modelo de la rama en la que estés.
- **Para correr un modelo:** cambiá a su rama (`git checkout t-model`), leé su `MODELO.md` y seguí [Cómo usarlo](#cómo-usarlo).
- **Para entender de dónde sale un número de las notas:** buscalo en `provenance/numbers.json`, que indica qué código lo produjo y qué elecciones se hicieron. `provenance/claims.yaml` dice qué afirmación respalda y cómo se verificó.
- **Para modificar el análisis:** el punto de entrada es `code/analisis_cosmolattice.py`, que usa `code/reescaleo_gws.py`. Si cambiás una fórmula, actualizá `code/verificar_reescaleo.py` y la sección correspondiente de las notas.

## Organización en ramas

Todo lo general está en `main`, y cada tipo de modelo tiene su rama, que es `main` más lo propio de ese modelo:

| rama | qué agrega a `main` |
|---|---|
| `main` | CosmoLattice, el código de análisis, las notas generales y este README |
| `monomial` | `models/genmonomial.h`, `genmonomial.in` ($n = 6$, exploratorio), `genmonomial_cuartico.in` (control $\phi^4$), `Notas/modelos/monomial.tex` y `MODELO.md` |
| `t-model` | `models/gentmodel.h`, `gentmodel.in` ($n = 2.5$, exploratorio), `gentmodel_cuartico.in` (caso del paper), `Notas/modelos/tmodel.tex` y `MODELO.md` |
| `e-model` | `models/genemodel.h`, `genemodel.in` ($n = 2.5$), `genemodel_cuartico.in` y `genemodel_cuartico_aK5.in` (casos del paper), `Notas/modelos/emodel.tex` y `MODELO.md` |

Además, en cada rama `provenance/` suma las afirmaciones propias de ese modelo.

Para trabajar con un modelo: `git checkout t-model` (o `e-model`, `monomial`). Los cambios generales se hacen en `main` y después se llevan a cada rama con `git merge main`. Si el merge da conflicto en `provenance/` (pasa cuando los dos lados agregaron entradas), se resuelve con `python3 code/resolver_procedencia.py` y después `git add provenance && git commit --no-edit`.

## Guía de carpetas

### Lo nuestro

| carpeta / archivo | qué tiene |
|---|---|
| `code/` | herramientas del estudio (ver la tabla de abajo) |
| `Notas/` | notas en LaTeX: `Notas.tex` (documento principal), `Librerías.tex` (paquetes y estilo) y, en cada rama, `modelos/<modelo>.tex` con la sección de su modelo |
| `provenance/` | registro de procedencia: `numbers.json` (cada número: qué código lo produjo, qué se hizo a mano, qué elecciones) y `claims.yaml` (cada afirmación: evidencia y cómo se verificó) |
| `models/gen*.h`, `models/parameter-files/gen*.in` | modelos generalizados y sus archivos de parámetros (sólo en las ramas de modelo) |
| `MODELO.md` | descripción de la rama y cómo correr su modelo (sólo en las ramas de modelo) |
| `README.md` | este archivo |
| `.claude/`, `.codex/` | configuración para asistentes de código: un hook que, al terminar cada tarea, verifica que `provenance/` esté completo y bien formado |

Contenido de `code/`:

| archivo | qué hace |
|---|---|
| `analisis_cosmolattice.py` | analiza una corrida: fondo, campos, espectros, errores, GWs y energías, con los espectros reescalados a hoy; deja gráficos y un `resumen.txt` en `analisis/` |
| `reescaleo_gws.py` | fórmulas del reescaleo a hoy (frecuencia y amplitud) y variables de programa de cada modelo |
| `verificar_reescaleo.py` | verificación con sympy de las ecuaciones de fondo y del reescaleo, y comparación con la literatura |
| `parametros_cuarticos.py` | traduce los parámetros de Ellis *et al.* a los `.in` y comprueba los `.in` presentes |
| `recursos_cosmolattice.py` | estima memoria y tiempo de cómputo según `N` |
| `correr_cuarticos.sh` | corre en secuencia los casos cuárticos de la rama, verificando antes la memoria, y analiza cada uno |
| `verificar_genmodels.py` | derivadas, condiciones iniciales y cotas de estabilidad de los modelos de la rama |
| `probar_genmodels.sh` | compila y prueba los modelos de la rama |
| `resolver_procedencia.py` | resuelve el conflicto de `provenance/` al hacer `git merge main` en una rama |

### CosmoLattice (sin modificar)

| carpeta / archivo | qué tiene |
|---|---|
| `include/CosmoInterface/` | el código de CosmoLattice: evolución, condiciones iniciales, mediciones, GWs (TempLat se descarga al compilar) |
| `source/cosmolattice.cpp` | programa principal |
| `models/*.h`, `models/parameter-files/*.in` | modelos originales de CosmoLattice (`lphi4`, `tanh2`, etc.); `lphi4` es la referencia para $n = 4$ |
| `CMakeLists.txt`, `cmake/` | compilación (`cmake .. -DMODEL=<modelo>`) |
| `tests/`, `profile/`, `format.sh` | pruebas, perfilado y formato del código de CosmoLattice |
| `README_CosmoLattice.md`, `LICENSE.md` | README y licencia originales de CosmoLattice |

### Generado al trabajar (no se versiona)

| carpeta / archivo | qué tiene |
|---|---|
| `build_<modelo>/` | compilación de cada modelo; el ejecutable es `build_<modelo>/<modelo>` |
| `build_corridas/<caso>_N<N>_kIR<kIR>/` | salidas de las corridas de `correr_cuarticos.sh` (los `average_*.txt` y `spectra_*.txt` de CosmoLattice) y su `analisis/` |
| `build_corridas/tiempos.txt` | tiempo real de cada corrida |
| `corridas_*.log` | registro de las corridas lanzadas con `nohup` |
| `Notas/Notas.pdf` | notas compiladas (cambian según la rama) |

## Cómo usarlo

Requisitos: lo que pide CosmoLattice (CMake ≥ 3.16 y un compilador C++20), Python 3 con `numpy`, `scipy`, `sympy`, `matplotlib` y `pyyaml`, y LaTeX (`latexmk`) para las notas.

```bash
git checkout t-model                                   # elegir el modelo
python3 code/recursos_cosmolattice.py --kmax-fijo --N 64 128 192   # memoria y tiempo según N
code/correr_cuarticos.sh 64                            # compila si hace falta, corre y analiza (~30 min por caso)
nohup code/correr_cuarticos.sh 128 > corridas_N128.log 2>&1 &      # producción, en segundo plano
python3 code/analisis_cosmolattice.py build_corridas/gentmodel_cuartico_N64_kIR0.7   # análisis de una corrida
cd Notas && latexmk -g -pdf Notas.tex                  # compilar las notas (-g: rehacer al cambiar de rama)
```

Con 3.7 GB de RAM entra hasta $N = 192$ con GWs (~1 GB). $N = 256$ necesita ~2.5 GB. Los detalles están en las notas, en las secciones «Corridas cuárticas» y «Guía de uso». Todos los parámetros de los `.in` se explican ahí, y cualquiera se puede cambiar desde la terminal, por ejemplo `N=128 kIR=0.35 baseSeed=1234`.

## Créditos

- **CosmoLattice**: D. G. Figueroa, A. Florio, F. Torrentí y W. Valkenburg, [arXiv:2006.15122](https://arxiv.org/abs/2006.15122) y [arXiv:2512.15627](https://arxiv.org/abs/2512.15627); módulo de GWs en la nota técnica II. Si usás este código, citá CosmoLattice como se indica en [cosmolattice.net](https://cosmolattice.net/CLcitation.html). El README original está en [`README_CosmoLattice.md`](README_CosmoLattice.md).
- **Parámetros de los modelos**: J. Ellis, M. A. G. Garcia, K. A. Olive y S. Verner, *Phys. Rev. D* **113**, 063571 (2026).
- **Reescaleo a hoy**: J.-F. Dufaux *et al.*, [arXiv:0707.0875](https://arxiv.org/abs/0707.0875); D. G. Figueroa y F. Torrentí, [arXiv:1707.04533](https://arxiv.org/abs/1707.04533).
