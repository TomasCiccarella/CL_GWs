# Rama `e-model`: $V = A\,(1 - e^{-\phi/M})^n$

Esta rama es `main` más el E-model generalizado ($\alpha$-Starobinsky). Lo general (análisis, reescaleo a hoy, recursos, guía) está en `main`; ver el [README](README.md).

## Qué agrega

| archivo | qué es |
|---|---|
| `models/genemodel.h` | modelo de CosmoLattice: $\tilde V = \lvert\tilde M(1 - e^{-\tilde\phi/\tilde M})\rvert^n/n + \tfrac12 q\tilde\phi^2\tilde\chi^2$, con $\omega_* = \sqrt{nA}\,M^{-n/2}\phi_0^{n/2-1}$ y $\alpha = 3(n-2)/(n+2)$ |
| `models/parameter-files/genemodel.in` | caso exploratorio $n = 2.5$, $M = 10\,m_p$ (sin GWs) |
| `models/parameter-files/genemodel_cuartico.in` | **caso del paper**: $n = 4$, $\alpha_K = 1$ ($M = \sqrt{3/2}\,m_p$, $\lambda = 1.53\times10^{-10}$; $n_s = 0.965$, $r = 0.004$), con GWs, `tMax = 600` |
| `models/parameter-files/genemodel_cuartico_aK5.in` | **caso del paper**: $n = 4$, $\alpha_K = 5$ ($\lambda = 7.12\times10^{-10}$; $n_s = 0.965$, $r = 0.015$, el más cercano a ACT), con GWs |
| `Notas/modelos/emodel.tex` | sección del modelo en las notas (derivadas, variables de programa, condiciones iniciales) |

El potencial es asimétrico: tiene plateau para $\phi \to +\infty$ y crece exponencialmente para $\phi \to -\infty$, así que $\phi_0$ tiene que ser positivo. Los parámetros cuárticos salen de Ellis, Garcia, Olive & Verner (2026); `code/parametros_cuarticos.py` (en `main`) los recalcula y comprueba los `.in`.

## Cómo correrlo

```bash
code/correr_cuarticos.sh 64            # corre los dos casos cuárticos y los analiza (~1.1 h)
MODELOS=genemodel_cuartico_aK5 code/correr_cuarticos.sh 128   # uno solo, en producción
python3 code/verificar_genmodels.py    # derivadas y condiciones iniciales
cd Notas && latexmk -g -pdf Notas.tex  # notas con la sección del E-model
```
