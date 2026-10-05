# Ondas gravitacionales del preheating en modelos de inflación favorecidos por las observaciones

Este repositorio estudia el **fondo estocástico de ondas gravitacionales (GWs) que se produce al final de la inflación**, durante el *preheating*, mediante simulaciones en la red con [CosmoLattice](https://cosmolattice.net). El objetivo es obtener el espectro de GWs **reescalado a hoy**, es decir la amplitud $h^2\Omega_{GW}(f)$ en función de la frecuencia $f$, para modelos de inflación que siguen siendo compatibles con las observaciones del CMB y para sus generalizaciones.

## La idea

Después de la inflación, el inflatón $\phi$ oscila alrededor del mínimo de su potencial. Si está acoplado a otro campo $\chi$ (acá, con $\tfrac12 g^2\phi^2\chi^2$), las oscilaciones amplifican exponencialmente las fluctuaciones de $\chi$ (**resonancia paramétrica**). Las inhomogeneidades resultantes son una fuente de ondas gravitacionales. Ese fondo de GWs queda como una huella de la física del final de la inflación, a frecuencias altas ($f \sim 10^8$–$10^9$ Hz).

Estudiamos potenciales que cerca del mínimo se comportan como $|\phi|^n$:

| modelo | potencial | comentario |
|---|---|---|
| monomial generalizado | $V = A\,\lvert\phi\rvert^n$ | caso de referencia; con $n = 4$ es el `lphi4` de CosmoLattice |
| T-model | $V = A\tanh^n(\phi/M)$ | $\alpha$-attractor con plateau |
| E-model | $V = A\,(1 - e^{-\phi/M})^n$ | $\alpha$-attractor ($\alpha$-Starobinsky generalizado) |

**Modelos favorecidos por las observaciones.** Los parámetros de referencia salen de Ellis, Garcia, Olive & Verner, *Phys. Rev. D* **113**, 063571 (2026) ([doi:10.1103/d35r-7bn8](https://doi.org/10.1103/d35r-7bn8)). Ese trabajo contrasta los T- y E-models generalizados con Planck, BICEP/Keck, ACT DR6 y SPT-3G. Nos concentramos en el **caso cuártico ($n = 4$)** por dos razones:

- **Es compatible con el CMB:** para $\alpha_K = 1$ da $n_s \simeq 0.965$ y $r \simeq 0.004$.
- **El resultado no depende del reheating:** el inflatón oscilando se comporta como radiación ($w = 1/3$). Por eso tanto el número de e-folds $N_*$ como el espectro de GWs de hoy no dependen de la temperatura de reheating, que es desconocida.

Los casos que se simulan son el T-model con $\alpha_K = 1$ y el E-model con $\alpha_K = 1$ y $5$. Se agrega un monomial $\phi^4$ de control (excluido por $r$, pero con el mismo mínimo que el T-model) para aislar el efecto del plateau.

**Del lattice a hoy.** El espectro que da CosmoLattice se lleva a hoy con

$$
f_0 = \frac{\tilde k}{a_e\,\tilde\rho_e^{1/4}}\left(\frac{\omega_*}{f_*}\right)^{1/2}\epsilon^{1/4}\times 4.59\times10^{10}\,\text{Hz},
\qquad
h^2\Omega_{GW,0} = 1.61\times10^{-5}\;\epsilon\;\Omega_{GW,e},
$$

con $\epsilon = (a_e/a_{RD})^{2(4-n)/(n+2)}$, que vale exactamente $1$ para $n = 4$. La derivación está verificada con sympy y comparada con Dufaux *et al.* (2007) y Figueroa & Torrentí (2017).

## Organización del repositorio

Todo lo general está en `main`, y cada tipo de modelo tiene su rama, que es `main` más lo propio de ese modelo:

| rama | qué agrega a `main` |
|---|---|
| `main` | CosmoLattice, el código de análisis, las notas generales y este README |
| `monomial` | `models/genmonomial.h`, `genmonomial.in` ($n = 6$, exploratorio), `genmonomial_cuartico.in` (control $\phi^4$) y su sección de las notas |
| `t-model` | `models/gentmodel.h`, `gentmodel.in` ($n = 2.5$, exploratorio), `gentmodel_cuartico.in` (caso del paper) y su sección de las notas |
| `e-model` | `models/genemodel.h`, `genemodel.in` ($n = 2.5$), `genemodel_cuartico.in` y `genemodel_cuartico_aK5.in` (casos del paper) y su sección de las notas |

Para trabajar con un modelo: `git checkout t-model` (o `e-model`, `monomial`). Los cambios generales se hacen en `main` y después se llevan a cada rama con `git merge main`.

### Qué hay en `main`

```
code/
  analisis_cosmolattice.py   análisis de una corrida: fondo, campos, espectros, errores, GWs y energías,
                             con los espectros reescalados a hoy
  reescaleo_gws.py           fórmulas del reescaleo a hoy (frecuencia y amplitud)
  verificar_reescaleo.py     verificación con sympy y comparación con la literatura
  parametros_cuarticos.py    traducción de los parámetros de Ellis et al. a los .in (y chequeo de los .in)
  recursos_cosmolattice.py   memoria y tiempo de cómputo según N
  correr_cuarticos.sh        corre en secuencia los casos cuárticos de la rama y los analiza
  verificar_genmodels.py     derivadas, condiciones iniciales y cotas de estabilidad de los modelos de la rama
  probar_genmodels.sh        compila y prueba los modelos de la rama
Notas/                       notas en LaTeX (en cada rama se agrega la sección de su modelo)
provenance/                  de dónde sale cada número y cada afirmación de las notas
models/, include/, source/   CosmoLattice (sin modificar) y sus modelos originales
```

## Cómo usarlo

Requisitos: lo que pide CosmoLattice (CMake ≥ 3.16 y un compilador C++20), Python 3 con `numpy`, `scipy`, `sympy`, `matplotlib` y `pyyaml`, y LaTeX (`latexmk`) para las notas.

```bash
git checkout t-model                                   # elegir el modelo
python3 code/recursos_cosmolattice.py --kmax-fijo --N 64 128 192   # memoria y tiempo según N
code/correr_cuarticos.sh 64                            # compila si hace falta, corre y analiza (~30 min por caso)
nohup code/correr_cuarticos.sh 128 > corridas_N128.log 2>&1 &      # producción, en segundo plano
python3 code/analisis_cosmolattice.py build_corridas/gentmodel_cuartico_N64_kIR0.7   # análisis de una corrida
cd Notas && latexmk -pdf Notas.tex                     # compilar las notas
```

Con 3.7 GB de RAM entra hasta $N = 192$ con GWs (~1 GB). $N = 256$ necesita ~2.5 GB. Los detalles están en las notas, en las secciones «Corridas cuárticas» y «Guía de uso». Todos los parámetros de los `.in` se explican ahí, y cualquiera se puede cambiar desde la terminal, por ejemplo `N=128 kIR=0.35 baseSeed=1234`.

## Créditos

- **CosmoLattice**: D. G. Figueroa, A. Florio, F. Torrentí y W. Valkenburg, [arXiv:2006.15122](https://arxiv.org/abs/2006.15122) y [arXiv:2512.15627](https://arxiv.org/abs/2512.15627); módulo de GWs en la nota técnica II. Si usás este código, citá CosmoLattice como se indica en [cosmolattice.net](https://cosmolattice.net/CLcitation.html). El README original está en [`README_CosmoLattice.md`](README_CosmoLattice.md).
- **Parámetros de los modelos**: J. Ellis, M. A. G. Garcia, K. A. Olive y S. Verner, *Phys. Rev. D* **113**, 063571 (2026).
- **Reescaleo a hoy**: J.-F. Dufaux *et al.*, [arXiv:0707.0875](https://arxiv.org/abs/0707.0875); D. G. Figueroa y F. Torrentí, [arXiv:1707.04533](https://arxiv.org/abs/1707.04533).
