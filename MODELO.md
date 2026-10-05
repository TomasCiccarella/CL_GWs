# Rama `t-model`: $V = A\tanh^n(\phi/M)$

Esta rama es `main` más el T-model generalizado. Lo general (análisis, reescaleo a hoy, recursos, guía) está en `main`; ver el [README](README.md).

## Qué agrega

| archivo | qué es |
|---|---|
| `models/gentmodel.h` | modelo de CosmoLattice: $\tilde V = \lvert\tilde M\tanh(\tilde\phi/\tilde M)\rvert^n/n + \tfrac12 q\tilde\phi^2\tilde\chi^2$, con $\omega_* = \sqrt{nA}\,M^{-n/2}\phi_0^{n/2-1}$ y $\alpha = 3(n-2)/(n+2)$ |
| `models/parameter-files/gentmodel.in` | caso exploratorio $n = 2.5$, $M = 10\,m_p$ (sin GWs) |
| `models/parameter-files/gentmodel_cuartico.in` | **caso del paper**: $n = 4$, $\alpha_K = 1$ ($M = \sqrt6\,m_p$, $\lambda = 1.60\times10^{-10}$; $n_s = 0.964$, $r = 0.004$), con GWs |
| `Notas/modelos/tmodel.tex` | sección del modelo en las notas (derivadas, variables de programa, condiciones iniciales) |

Los parámetros del caso cuártico salen de Ellis, Garcia, Olive & Verner (2026); `code/parametros_cuarticos.py` (en `main`) los recalcula y comprueba el `.in`.

## Cómo correrlo

```bash
code/correr_cuarticos.sh 64            # corre gentmodel_cuartico.in y lo analiza (~30 min)
code/correr_cuarticos.sh 128           # producción: kIR = 0.35, ~4 h
python3 code/verificar_genmodels.py    # derivadas y condiciones iniciales
cd Notas && latexmk -g -pdf Notas.tex  # notas con la sección del T-model
```

Para comparar con el monomial de control (rama `monomial`) usá la misma `N`, el mismo `kIR` y la misma `baseSeed`.
