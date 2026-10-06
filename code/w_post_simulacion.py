#!/usr/bin/env python3
"""Ecuación de estado después de la simulación para los casos cuadráticos (n = 2), y el factor de
dilución de las GWs hasta el reheating que resulta de ella (Notas.tex, «Corridas cuadráticas»,
«Qué pasa después de la simulación»).

Uso:
    python3 code/w_post_simulacion.py [DIR ...] [--TRH 1e10] [--out FIGURA]

Sin directorios, toma las cuatro corridas cuadráticas con N = 64. Para cada una:

  1. Del último espectro de cada campo (Delta_phi~ y Delta_phi~' por log k) arma la energía por modo
     rho_k = (Delta_pi + omega_k^2 Delta_phi)/2, con omega_k^2 = k~^2/a^2 + m~_eff^2, y el número de
     partículas n_k = rho_k/omega_k, que es un invariante adiabático una vez que los campos son libres.
  2. Después de la simulación, cada modo es una partícula libre: n_k a^3 se conserva, el momento físico
     cae como 1/a y la energía es omega_k(a) = sqrt(k~^2/a^2 + m~_eff^2(a)). Las masas efectivas son
     las de campo medio, m~_phi^2 = 1 + q <chi~^2> y m~_chi^2 = q <phi~^2>, recalculadas a cada a con
     <f~^2> = sum n_k/omega_k (a_s/a)^3 (iteración de punto fijo). El condensado homogéneo de phi
     (lo que queda) se diluye como materia.
  3. Con rho(a) y p(a) = sum n_k (k~^2/a^2)/(3 omega_k) sale w(a). El reheating ocurre cuando
     rho(a_RD) = (pi^2 g_*/30) T_RH^4, y el factor de dilución de las GWs es
         eps = (rho_s/rho_RD) (a_s/a_RD)^4,
     que se reduce a (a_s/a_RD)^(1-3w) si w es constante.

Supuestos (que la figura y la tabla no pueden verificar): los campos dejan de interactuar salvo por
las masas de campo medio (sin dispersión ni decaimiento hasta T_RH), el reheating es instantáneo en
a_RD, y se usa la relación de dispersión del continuo (con N = 64 el UV no está resuelto y la energía
de gradiente del continuo es ~40% mayor que la de la red al final de la corrida). La energía del modelo
se normaliza a la total de la simulación en a_s, para usar sólo su forma.
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
from analisis_cosmolattice import CAT, leer_infos, leer_tabla, leer_espectros  # noqa: E402
from comparar_corridas import CASOS_N64_CUADRATICOS, RAIZ  # noqa: E402


def modos_finales(d):
    """Estado al último tiempo de espectros: a_s, rho~_s, modos de cada campo y condensado."""
    modelo, p = leer_infos(d)
    fs, om, alpha, n = rg.variables_programa(modelo, p)
    if n != 2 or alpha != 0:
        sys.exit(f"{d}: el modelo de partículas libres está pensado para n = 2 (alpha = 0)")
    q, dk = p["q"], p["kIR"] * p.get("deltaKBin", 1.0)
    S = leer_tabla(os.path.join(d, "average_scale_factor.txt"))
    E = leer_tabla(os.path.join(d, "average_energies.txt"))
    ts = np.atleast_1d(np.loadtxt(os.path.join(d, "average_spectra_times.txt"), comments="#"))
    t_s = ts[-1]
    a_s = np.interp(t_s, S["t"], S["a"])
    rho_s = np.interp(t_s, E["t"], E["E_tot"])
    campos = []
    for i in (0, 1):
        b = leer_espectros(os.path.join(d, f"spectra_scalar_{i}.txt"))[len(ts) - 1]
        campos.append(dict(k=b[:, 0], Dphi=b[:, 1] * dk / b[:, 0], Dpi=b[:, 2] * dk / b[:, 0]))
    F0 = leer_tabla(os.path.join(d, "average_scalar_0.txt"))
    j = np.argmin(np.abs(F0["t"] - t_s))
    return dict(fs=fs, om=om, q=q, a_s=a_s, rho_s=rho_s, campos=campos,
                phi_hom2=F0["<phi>"][j] ** 2, pi_hom2=F0["<pi>"][j] ** 2)


def masas(var_phi, var_chi, q):
    """m~_eff^2 de campo medio: phi tiene la masa desnuda (1) más q<chi^2>; chi sólo q<phi^2>."""
    return 1 + q * var_chi, q * var_phi


def evolucion(M, a):
    """rho~(a) y p~(a) del modelo de partículas libres, para a >= a_s (a puede ser un array).

    Con las masas de campo medio, sum n omega cuenta dos veces la interacción U = (q/2)<phi^2><chi^2>
    (una vez en cada masa). La energía es rho = sum n omega - U y la presión p = sum n k^2/(3 a^2 omega)
    + U: con esas dos, d(a^3 rho)/d ln a = -3 a^3 p exactamente (la variación de las masas con a
    trabaja contra U). El condensado entra como el modo k = 0 de phi, con la misma masa efectiva.
    """
    q, a_s = M["q"], M["a_s"]
    ph, ch = M["campos"]
    k = (np.append(0.0, ph["k"]), ch["k"])
    Dphi = (np.append(M["phi_hom2"], ph["Dphi"]), ch["Dphi"])
    Dpi = (np.append(M["pi_hom2"], ph["Dpi"]), ch["Dpi"])
    # Números de partículas en a_s, con las masas de las varianzas medidas.
    var = [np.sum(D) for D in Dphi]
    m2 = masas(*var, q)
    n = [0.5 * (Dpi[i] + (k[i] ** 2 / a_s ** 2 + m2[i]) * Dphi[i]) / np.sqrt(k[i] ** 2 / a_s ** 2 + m2[i])
         for i in (0, 1)]
    rho, pre = [], []
    for x in np.append(a_s, a):  # a_s primero, para normalizar
        dil = (a_s / x) ** 3
        v = [var[0] * dil, var[1] * dil]
        for _ in range(100):  # punto fijo de las masas efectivas
            m2 = masas(*v, q)
            om = [np.sqrt(k[i] ** 2 / x ** 2 + m2[i]) for i in (0, 1)]
            v_new = [np.sum(n[i] / om[i]) * dil for i in (0, 1)]
            if all(abs(v_new[i] - v[i]) <= 1e-12 * v[i] for i in (0, 1)):
                break
            v = v_new
        U = 0.5 * q * v[0] * v[1]
        rho.append(sum(np.sum(n[i] * om[i]) for i in (0, 1)) * dil - U)
        pre.append(sum(np.sum(n[i] * k[i] ** 2 / (3 * x ** 2 * om[i])) for i in (0, 1)) * dil + U)
    rho, pre = np.array(rho), np.array(pre)
    # Se usa sólo la forma: la energía en a_s se normaliza a la total de la simulación.
    escala = M["rho_s"] / rho[0]
    return rho[1:] * escala, pre[1:] * escala


def reheating(M, TRH, g_star=rg.G_SM):
    """a_RD (rho = rho_RD), eps y la curva w(a) entre a_s y a_RD."""
    rho_RD = np.pi ** 2 * g_star / 30 * TRH ** 4 / (M["fs"] ** 2 * M["om"] ** 2)  # en unidades de programa
    if rho_RD >= M["rho_s"]:
        return dict(a_RD=M["a_s"], eps=1.0, a=np.array([M["a_s"]]), w=np.array([np.nan]))
    a = M["a_s"] * np.logspace(0, 40, 2001)
    rho, pre = evolucion(M, a)
    i = np.argmax(rho <= rho_RD)
    if rho[i] > rho_RD:
        sys.exit("No se llega a rho_RD: ampliá el rango de a")
    # interpolación en log
    la = np.interp(np.log(rho_RD), np.log(rho[i - 1:i + 1][::-1]), np.log(a[i - 1:i + 1][::-1]))
    a_RD = np.exp(la)
    eps = (M["rho_s"] / rho_RD) * (M["a_s"] / a_RD) ** 4
    return dict(a_RD=a_RD, eps=eps, a=a[:i + 1], w=(pre / rho)[:i + 1], rho_RD=rho_RD)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="*")
    ap.add_argument("--TRH", type=float, default=1e10)
    ap.add_argument("--out", default=os.path.join(RAIZ, "Notas", "figuras", "w_post_cuadraticos.pdf"))
    args = ap.parse_args()
    if args.dirs:
        casos = [(d, os.path.basename(d.rstrip("/"))) for d in args.dirs]
    else:
        casos = [(os.path.join(RAIZ, "build_corridas", f"{c}_N64_kIR3"), e) for c, e in CASOS_N64_CUADRATICOS]
    fig, ax = plt.subplots(figsize=(6.2, 3.6), constrained_layout=True)
    print(f"T_RH = {args.TRH:.3g} GeV")
    print(f"{'caso':>26} {'w(a_s)':>7} {'a(w<0.1)/a_s':>12} {'a_RD/a_s':>10} {'eps':>10} {'eps(w=0)':>10} {'eps(w=.25)':>10}")
    res = {}
    for i, (d, et) in enumerate(casos):
        M = modos_finales(d)
        R = reheating(M, args.TRH)
        w = R["w"]
        x = R["a"] / M["a_s"]
        a01 = x[np.argmax(w < 0.1)] if np.any(w < 0.1) else np.nan
        # Comparación con w constante (mismo rho_RD): eps = (rho_RD/rho_s)^((1-3w)/(3(1+w))).
        epsw = {w0: (R["rho_RD"] / M["rho_s"]) ** ((1 - 3 * w0) / (3 * (1 + w0))) for w0 in (0.0, 0.25)}
        print(f"{os.path.basename(d.rstrip('/')):>26} {w[0]:7.3f} {a01:12.3g} {R['a_RD'] / M['a_s']:10.3g} "
              f"{R['eps']:10.3e} {epsw[0.0]:10.3e} {epsw[0.25]:10.3e}")
        # Espectro de GWs de hoy con este eps (mismo tiempo de espectros que los modos). Para n = 2,
        # reescaleo_gws.epsilon(2, x) = x, así que pasar eps como a_e/a_RD da exactamente eps.
        g = leer_espectros(os.path.join(d, "spectra_energy_gws.txt"))[-1]
        f0 = rg.frecuencia_hoy(g[:, 0], M["a_s"], M["rho_s"], M["om"], M["fs"], 2, R["eps"])
        Om0 = rg.omega_gw_hoy(g[:, 1], 2, R["eps"])
        jp = int(np.argmax(Om0))
        banda = f0[Om0 >= Om0[jp] / 2]
        print(f"{'':>26} hoy: máximo h^2 Omega_GW,0 = {Om0[jp]:.2e} en {f0[jp]:.2e} Hz; "
              f"más de la mitad entre {banda.min():.2e} y {banda.max():.2e} Hz")
        ax.semilogx(x, w, color=CAT[i], lw=1.5, label=et)
        res[d] = dict(M=M, R=R, a01=a01, epsw=epsw, Om_max=Om0[jp], f_max=f0[jp], banda=(banda.min(), banda.max()))
    ax.axhline(1 / 3, color="#52514e", lw=0.8, ls="--")
    ax.axhline(0, color="#52514e", lw=0.8)
    ax.set(xlabel=r"$a/a_s$ (desde el final de la simulación hasta $a_{RD}$)", ylabel=r"$w = p/\rho$",
           title=rf"Partículas libres después de la simulación ($T_{{\rm RH}} = 10^{{{np.log10(args.TRH):.3g}}}$ GeV)", ylim=(-0.02, 0.36))
    ax.legend(fontsize=8)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fig.savefig(args.out)
    print(f"\nFigura: {args.out}")
    return res


if __name__ == "__main__":
    main()
