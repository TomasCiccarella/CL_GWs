# Rama `monomial`: $V = A\,|\phi|^n$

Esta rama es `main` más el modelo monomial generalizado. Lo general (análisis, reescaleo a hoy, recursos, guía) está en `main`; ver el [README](README.md).

## Qué agrega

| archivo | qué es |
|---|---|
| `models/genmonomial.h` | modelo de CosmoLattice: $\tilde V = \lvert\tilde\phi\rvert^n/n + \tfrac12 q\tilde\phi^2\tilde\chi^2$, $\omega_* = \sqrt{nA}\,\phi_0^{n/2-1}$, $\alpha = 3(n-2)/(n+2)$ |
| `models/parameter-files/genmonomial.in` | caso exploratorio $n = 6$ (sin GWs); ilustra la inestabilidad para $n > 4$ |
| `models/parameter-files/genmonomial_cuartico.in` | **control** $\phi^4$ con GWs: el $\lambda_{\rm eff}$ y el $q$ del T-model cuártico, sin plateau |
| `Notas/modelos/monomial.tex` | sección del modelo en las notas (derivadas, condiciones iniciales, ejemplo $n = 6$, control) |

El control **no es un modelo del paper**: $\phi^4$ puro predice $n_s = 0.946$ y $r = 0.29$, que $r < 0.036$ excluye. Sirve para comparar con el T-model de la rama `t-model` y aislar el efecto del plateau. Para esa comparación conviene usar la misma `N`, el mismo `kIR` y la misma `baseSeed`.

## Cómo correrlo

```bash
code/correr_cuarticos.sh 64            # corre genmonomial_cuartico.in y lo analiza (~30 min)
python3 code/verificar_genmodels.py    # derivadas, condiciones iniciales y cotas de estabilidad
cd Notas && latexmk -pdf Notas.tex     # notas con la sección del monomial
```
