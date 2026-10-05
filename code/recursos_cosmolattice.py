#!/usr/bin/env python3
"""Memoria y tiempo de cómputo de una corrida de CosmoLattice en función de N (ver Notas.tex Sec. 6).

Uso:
    python3 code/recursos_cosmolattice.py                      # tabla para N = 32 ... 256
    python3 code/recursos_cosmolattice.py --N 128 --kIR 0.35   # un caso
    python3 code/recursos_cosmolattice.py --sin-gws            # sin el módulo de GWs

Modelo (constantes medidas en esta máquina, i7-8565U 4 núcleos / 8 hilos, 3.7 GB, con 2 escalares;
ver provenance/numbers.json):
    memoria  RSS ~ M0 + b N^3, con b = 148 B/sitio con GWs (63 B/sitio sin GWs) y M0 ~ 10 MB
    tiempo   t_paso ~ c N^3, con c ~ 140 ns/sitio con GWs y 8 hilos (N >= 128; 122 ns en N = 64),
             84 ns/sitio sin GWs
    pasos    tMax/dt
    dt       el menor entre dt_max (precisión, 0.01 por defecto) y 0.75/omega~_max, con
             omega~_max^2 = (2 sqrt(3)/dx)^2 + q  y  dx = 2 pi/(N kIR)
             (frecuencia máxima de la grilla: laplaciano de 2do orden, más la masa efectiva de chi
              con amplitud conforme a phi~ <= 1; para n = 4, alpha = 1 no crecen con a).
             VV2/LF son estables si omega~_max dt < 2. Con dt = 0.01 se probó hasta omega~_max dt = 0.75
             (la dx de N = 192 con kIR = 0.7): rho_GW/rho cambia < 1% respecto de dt = 0.005. El error
             de Friedmann escala como dt^2 y crece con N en la etapa no lineal (5.6e-4 con N = 32,
             1.3e-3 con N = 64, ambos con dt = 0.01).
"""

import argparse

import numpy as np

B_GW, B_SIN = 148.0, 63.0      # bytes por sitio (RSS medido, 2 escalares)
M0 = 10e6                      # bytes fijos
C_PASO = 140e-9                # s por sitio y por paso, con GWs y 8 hilos
C_PASO_SIN = 84e-9             # s por sitio y por paso, sin GWs y 8 hilos (medido en N = 128)


def memoria_libre():
    """MemAvailable de /proc/meminfo, en bytes (None si no se puede leer)."""
    try:
        for linea in open("/proc/meminfo"):
            if linea.startswith("MemAvailable:"):
                return float(linea.split()[1]) * 1024
    except OSError:
        pass
    return None


def omega_max(N, kIR, q):
    dx = 2 * np.pi / (N * kIR)
    return np.sqrt((2 * np.sqrt(3) / dx) ** 2 + q)


def dt_recomendado(N, kIR, q, dt_max=0.01, courant=0.75):
    return min(dt_max, courant / omega_max(N, kIR, q))


def recursos(N, kIR, q, tMax, gws=True, dt=None, dt_max=0.01, n_espectros=None):
    dt = dt or dt_recomendado(N, kIR, q, dt_max)
    sitios = N ** 3
    mem = M0 + (B_GW if gws else B_SIN) * sitios
    pasos = tMax / dt
    t_paso = (C_PASO if gws else C_PASO_SIN) * sitios
    # cada espectro cuesta ~3.4 pasos (medido en N = 128)
    t_total = t_paso * (pasos + 3.4 * (n_espectros if n_espectros else tMax / 10))
    return dict(N=N, kIR=kIR, k_max=kIR * N * np.sqrt(3) / 2, dt=dt, estable=omega_max(N, kIR, q) * dt < 2,
                mem_GB=mem / 1e9, pasos=pasos, t_paso_s=t_paso, t_total_h=t_total / 3600)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, nargs="*", default=[32, 64, 96, 128, 160, 192, 224, 256])
    ap.add_argument("--kIR", type=float, default=0.7)
    ap.add_argument("--kmax-fijo", action="store_true",
                    help="escalar kIR = 0.7*64/N (mismo k_max que N = 64, caja más grande)")
    ap.add_argument("--q", type=float, default=120)
    ap.add_argument("--tMax", type=float, default=500)
    ap.add_argument("--dt", type=float, default=None)
    ap.add_argument("--sin-gws", action="store_true")
    a = ap.parse_args()
    libre = memoria_libre()
    print(f"Memoria disponible ahora: {libre / 1e9:.2f} GB" if libre else "Memoria disponible: ?")
    print(f"{'N':>5} {'kIR':>6} {'k_max':>6} {'dt':>8} {'memoria':>9} {'pasos':>8} {'s/paso':>7} {'total':>9}  ")
    for N in a.N:
        kIR = 0.7 * 64 / N if a.kmax_fijo else a.kIR
        r = recursos(N, kIR, a.q, a.tMax, not a.sin_gws, a.dt)
        aviso = ""
        if libre and r["mem_GB"] * 1e9 > 0.8 * libre:
            aviso = "  <- no entra en la memoria libre (swap)"
        if not r["estable"]:
            aviso += "  <- dt inestable"
        print(f"{N:5d} {kIR:6.3f} {r['k_max']:6.1f} {r['dt']:8.5f} {r['mem_GB']:7.2f}GB {r['pasos']:8.0f} "
              f"{r['t_paso_s']:7.3f} {r['t_total_h']:7.2f} h{aviso}")


if __name__ == "__main__":
    main()
