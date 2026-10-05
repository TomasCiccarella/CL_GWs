#!/usr/bin/env python3
"""Compara varias corridas de CosmoLattice: espectro de GWs final (en la red y hoy) y evolución.

Uso:
    python3 code/comparar_corridas.py [DIR ...] [--out FIGURA] [--etiquetas E1 E2 ...]

Sin directorios, toma las cuatro corridas cuárticas con N = 64 de build_corridas/ (T, E alpha_K = 1,
E alpha_K = 5 y el monomial de control). Escribe la figura (por defecto
Notas/figuras/corridas_N64.pdf) e imprime una tabla con lo que se cita en las notas:
    pico de h^2 Omega_GW,0 y su frecuencia, ancho del espectro (dónde supera la mitad del pico),
    rho_GW/rho final, t~ del máximo de rms(chi~), t~ en que rho_GW/rho llega al 90% del final y
    cuánto crece en los últimos 100 de t~.
El reescaleo a hoy es el de analisis_cosmolattice.py (g_* = 106.75, eps = 1 para n = 4).
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reescaleo_gws as rg  # noqa: E402
from analisis_cosmolattice import CAT, TINTA2, leer_infos, leer_tabla, leer_espectros  # noqa: E402

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASOS_N64 = [("gentmodel_cuartico", r"T, $\alpha_K = 1$"),
             ("genemodel_cuartico", r"E, $\alpha_K = 1$"),
             ("genemodel_cuartico_aK5", r"E, $\alpha_K = 5$"),
             ("genmonomial_cuartico", r"$\phi^4$ (control)")]


def cargar(d):
    """Lo necesario de una corrida para compararla."""
    modelo, p = leer_infos(d)
    fs, om, alpha, n = rg.variables_programa(modelo, p)
    S = leer_tabla(os.path.join(d, "average_scale_factor.txt"))
    E = leer_tabla(os.path.join(d, "average_energies.txt"))
    G = leer_tabla(os.path.join(d, "average_energies_gws.txt"))
    X = leer_tabla(os.path.join(d, "average_scalar_1.txt"))
    b = leer_espectros(os.path.join(d, "spectra_energy_gws.txt"))[-1]
    a_e, rho_e = S["a"][-1], E["E_tot"][-1]
    k, Om_e = b[:, 0], b[:, 1]
    f0 = rg.frecuencia_hoy(k, a_e, rho_e, om, fs, n)
    Om0 = rg.omega_gw_hoy(Om_e, n)
    t, r = G["t"], G["rhoGW_over_rho"]
    jp = int(np.argmax(Om0))
    medio = k[Om0 >= Om0[jp] / 2]
    return dict(k=k, Om_e=Om_e, f0=f0, Om0=Om0, t=t, r=r,
                f_pico=f0[jp], k_pico=k[jp], Om_pico=Om0[jp], k_medio=(medio.min(), medio.max()),
                Hz_por_k=f0[0] / k[0], rgw=r[-1],
                t_chi=X["t"][int(np.argmax(X["rms(phi)"]))],
                t90=t[np.argmax(r >= 0.9 * r[-1])],
                ult100=100 * (r[-1] / np.interp(t[-1] - 100, t, r) - 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="*")
    ap.add_argument("--etiquetas", nargs="*")
    ap.add_argument("--out", default=os.path.join(RAIZ, "Notas", "figuras", "corridas_N64.pdf"))
    args = ap.parse_args()
    if args.dirs:
        dirs, etiquetas = args.dirs, args.etiquetas or [os.path.basename(d.rstrip("/")) for d in args.dirs]
    else:
        dirs = [os.path.join(RAIZ, "build_corridas", f"{c}_N64_kIR0.7") for c, _ in CASOS_N64]
        etiquetas = [e for _, e in CASOS_N64]
    C = [cargar(d) for d in dirs]

    fig, ax = plt.subplots(1, 3, figsize=(10, 3.3), constrained_layout=True)
    for i, (c, et) in enumerate(zip(C, etiquetas)):
        col = CAT[i]
        ax[0].loglog(c["k"], c["Om_e"], color=col, lw=1.4, label=et)
        ax[1].loglog(c["f0"], c["Om0"], color=col, lw=1.4)
        ax[1].plot(c["f_pico"], c["Om_pico"], "o", ms=7, color=col, mec="white", mew=1.5)
        ax[2].semilogy(c["t"], np.where(c["r"] > 0, c["r"], np.nan), color=col, lw=1.4)
        ax[2].axvline(c["t_chi"], color=col, lw=0.8, ls=":")
    ax[0].set(xlabel=r"$\tilde k$", ylabel=r"$\Omega_{GW,e}(\tilde k)$",
              title="(a) Espectro final en la red", ylim=(1e-10, 1e-4))
    ax[1].set(xlabel=r"$f_0$ [Hz]", ylabel=r"$h^2\Omega_{GW,0}(f)$",
              title="(b) Hoy", ylim=(1e-15, 3e-10))
    ax[2].set(xlabel=r"$\tilde t$", ylabel=r"$\rho_{GW}/\rho$",
              title=r"(c) Producción de GWs", ylim=(1e-12, 1e-4))
    ax[0].legend(loc="lower left", fontsize=8)
    for a in ax:
        a.tick_params(which="minor", color=TINTA2)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fig.savefig(args.out)
    print(f"Figura: {args.out}\n")

    cab = ["caso", "f_pico[Hz]", "k~_pico", "h2Om_pico", "k~(>pico/2)", "Hz/k~", "rhoGW/rho",
           "t~_chi", "t~_90", "ult100[%]"]
    print("  ".join(f"{h:>12}" for h in cab))
    for c, d in zip(C, dirs):
        fila = [os.path.basename(d.rstrip("/")).replace("_N64_kIR0.7", ""), f"{c['f_pico']:.3e}",
                f"{c['k_pico']:.3g}", f"{c['Om_pico']:.3e}", f"{c['k_medio'][0]:.3g}-{c['k_medio'][1]:.3g}",
                f"{c['Hz_por_k']:.4e}", f"{c['rgw']:.4e}", f"{c['t_chi']:.4g}", f"{c['t90']:.4g}",
                f"{c['ult100']:.2g}"]
        print("  ".join(f"{x:>12}" for x in fila))


if __name__ == "__main__":
    main()
