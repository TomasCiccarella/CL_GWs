#!/usr/bin/env python3
"""Variación entre semillas y efecto de N = 128 en los casos cuárticos (Notas.tex, «Semillas y resolución»).

Uso:
    python3 code/semillas_resolucion.py

Lee de build_corridas/:
  - el T-model con N = 64 y cinco semillas (la original y baseSeed = 20261005, 1, 2, 3);
  - los cuatro casos con N = 64 (semilla original y 20261005) y con N = 128, kIR = 0.35 (20261005).
Imprime lo que citan las Tablas tab:semillas-T y tab:semillas-N128: media y desviación estándar relativa
entre semillas, el cambio de rho_GW/rho de N = 64 a N = 128 con la misma semilla, y la fracción de
rho_GW en cada rango de k~ (para ver si la diferencia está en el IR, en el pico o en el UV). Devuelve
un dict con todo, para registrar la procedencia.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from comparar_corridas import CASOS_N64, RAIZ, cargar  # noqa: E402

B = os.path.join(RAIZ, "build_corridas")
SEMILLA = 20261005
SEMILLAS_T = ["", f"_s{SEMILLA}", "_s1", "_s2", "_s3"]  # "" = la corrida original
RANGOS = [(0, 1.05), (1.05, 3), (3, 8), (8, 20), (20, np.inf)]


def fracciones(c):
    """Fracción de rho_GW/rho en cada rango de k~ (integral de Omega_GW,e d ln k con bins uniformes)."""
    k, O = c["k"], c["Om_e"]
    peso = O * (k[1] - k[0]) / k
    return np.array([peso[(k >= lo) & (k < hi)].sum() for lo, hi in RANGOS]) / peso.sum()


def main():
    res = {}
    T = [cargar(os.path.join(B, f"gentmodel_cuartico_N64_kIR0.7{s}")) for s in SEMILLAS_T]
    print("T-model, N = 64, cinco semillas")
    for clave in ("rgw", "Om_pico", "t_chi"):
        x = np.array([c[clave] for c in T])
        m, d = x.mean(), x.std(ddof=1)
        res[f"T_{clave}"] = dict(valores=x, media=m, desv=d, desv_rel=d / m)
        print(f"  {clave:8s} " + " ".join(f"{v:.4g}" for v in x) + f"   media {m:.4g}, desv {d:.2g} ({100 * d / m:.1f}%)")
    print("  k~ del pico: " + " ".join(f"{c['k_pico']:.3g}" for c in T)
          + "; banda: " + " ".join(f"{c['k_medio'][0]:.3g}-{c['k_medio'][1]:.3g}" for c in T))

    print("\nMisma semilla: N = 64 -> N = 128")
    for caso, _ in CASOS_N64:
        o = cargar(os.path.join(B, f"{caso}_N64_kIR0.7"))
        a = cargar(os.path.join(B, f"{caso}_N64_kIR0.7_s{SEMILLA}"))
        b = cargar(os.path.join(B, f"{caso}_N128_kIR0.35_s{SEMILLA}"))
        fa, fb = fracciones(a), fracciones(b)
        res[caso] = dict(orig=o, n64=a, n128=b, cambio=b["rgw"] / a["rgw"] - 1, fr64=fa, fr128=fb,
                         dfr_max=np.abs(fb - fa).max(),
                         banda_Hz=(b["k_medio"][0] * b["Hz_por_k"], b["k_medio"][1] * b["Hz_por_k"]))
        print(f"  {caso:24s} rgw {o['rgw']:.4g} / {a['rgw']:.4g} / {b['rgw']:.4g} ({100 * res[caso]['cambio']:+.1f}%)"
              f"  k_pico {a['k_pico']:.3g}->{b['k_pico']:.3g}  ult100 {a['ult100']:.2g}/{b['ult100']:.2g}%")
        print(f"  {'':24s} fracciones por rango de k~: " + " ".join(f"{x:.2f}" for x in fa) + " -> "
              + " ".join(f"{x:.2f}" for x in fb) + f"  (cambio máx {res[caso]['dfr_max']:.3f})")
        print(f"  {'':24s} banda hoy N = 128: {res[caso]['banda_Hz'][0]:.3g}-{res[caso]['banda_Hz'][1]:.3g} Hz,"
              f" pico {b['Om_pico']:.3g}")
    return res


if __name__ == "__main__":
    main()
