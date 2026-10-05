#!/usr/bin/env python3
"""Parámetros de los alpha-attractors cuárticos (k = 4) favorecidos por Ellis, Garcia, Olive y
Verner, PRD 113, 063571 (2026) [arXiv del paper: d35r-7bn8], traducidos a la convención de
gentmodel/genemodel (ver Notas/Notas.tex, Sec. "Parámetros cosmológicos para k = 4").

Uso: python3 code/parametros_cuarticos.py

Convención del paper (Eqs. 39-40), con alpha_att el alpha de la métrica de Kähler:
    E:  V = (3/4) lam M_P^4 (1 - exp(-sqrt(2/(3 alpha_att)) phi/M_P))^k
    T:  V = (3/4) lam M_P^4 tanh^k(phi/(sqrt(6 alpha_att) M_P))
Convención de los .h:  V = A |f(phi/M)|^n, así que
    n = k,  A = (3/4) lam M_P^4,  M = sqrt(3 alpha_att/2) M_P (E)  o  sqrt(6 alpha_att) M_P (T).

Para cada caso:
  1. phi_end: epsilon_H = 1, integrando el fondo homogéneo exacto (no la aproximación B11/B12).
  2. phi_*: N_* e-folds de slow-roll antes de phi_end (Eq. 18, integral numérica).
  3. lam: normalización A_s = V/(24 pi^2 eps_V M_P^4) en phi_* (A_s = 2.1e-9, como el paper).
  4. n_s, r en slow-roll (n_s = 1 - 6 eps + 2 eta, r = 16 eps), para compararlos con la
     Tabla II del paper (N_* = 55, cálculo exacto de perturbaciones).
  5. Condiciones iniciales de CosmoLattice con las fórmulas de Notas.tex: phi_0 (eps_V = 1),
     phi'_0 = -sqrt((sqrt(7/3)-1) V(phi_0)) y omega_* = sqrt(nA) M^(-n/2) phi_0^(n/2-1).
"""

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

MP = 2.435e18   # GeV, la misma masa de Planck reducida que los .in
AS = 2.1e-9     # amplitud escalar usada por el paper (Sec. I, Eq. 38)
K = 4

# (modelo, alpha_att, N_*) -- N_* de la Ec. (B15) del paper para k = 4 (independiente de T_RH):
# 55.7 (E) y 55.8 (T) para alpha_att = 1 con g_RH del MSSM.
CASOS = [("E", 1.0, 55.7), ("T", 1.0, 55.8), ("E", 5.0, 55.7)]


def Mesc(mod, al):
    """Escala M (en unidades de M_P) de la convención de los .h."""
    return np.sqrt(1.5 * al) if mod == "E" else np.sqrt(6 * al)


# Potencial con lam = 1 en unidades M_P = 1 (lam sólo fija la escala, no la forma).
def f(mod, x, M):
    return 1 - np.exp(-x / M) if mod == "E" else np.tanh(x / M)


def df(mod, x, M):
    return np.exp(-x / M) / M if mod == "E" else 1 / (M * np.cosh(x / M) ** 2)


def d2f(mod, x, M):
    if mod == "E":
        return -np.exp(-x / M) / M ** 2
    return -2 * np.tanh(x / M) / (M ** 2 * np.cosh(x / M) ** 2)


def V(mod, x, M, k=K):
    return 0.75 * f(mod, x, M) ** k


def dV(mod, x, M, k=K):
    return 0.75 * k * f(mod, x, M) ** (k - 1) * df(mod, x, M)


def d2V(mod, x, M, k=K):
    u = f(mod, x, M)
    return 0.75 * k * ((k - 1) * u ** (k - 2) * df(mod, x, M) ** 2 + u ** (k - 1) * d2f(mod, x, M))


def epsV(mod, x, M, k=K):
    return 0.5 * (dV(mod, x, M, k) / V(mod, x, M, k)) ** 2


def phi_end_exacto(mod, M, k=K):
    """Integra el fondo exacto (M_P = 1, lam = 1) desde eps_V = 1e-4 (plateau) hasta eps_H = 1."""
    def rhs(t, y):
        H = np.sqrt((V(mod, y[0], M, k) + y[1] ** 2 / 2) / 3)
        return [y[1], -3 * H * y[1] - dV(mod, y[0], M, k)]

    def ev(t, y):
        return 1.5 * y[1] ** 2 / (V(mod, y[0], M, k) + y[1] ** 2 / 2) - 1
    ev.terminal = True
    xi = brentq(lambda x: epsV(mod, x, M, k) - 1e-4, 0.5, 40 * M)
    Hi = np.sqrt(V(mod, xi, M, k) / 3)
    sol = solve_ivp(rhs, [0, 1e6], [xi, -dV(mod, xi, M, k) / (3 * Hi)], events=ev,
                    rtol=1e-11, atol=1e-14)
    return float(sol.y_events[0][0][0]), float(sol.y_events[0][0][1])


def calcular(mod, al, Nst, k=K):
    """Parámetros del caso (mod, alpha_K, N_*) con exponente k. Nst puede ser una función de
    (phi_*, phi_end, phi'_end) en unidades de M_P y lam = 1 (para k = 2, donde N_* depende de T_RH):
    se resuelve N(phi_*) = Nst(phi_*, ...) de forma autoconsistente."""
    M = Mesc(mod, al)
    xe, pe = phi_end_exacto(mod, M, k)

    def N(x):  # Eq. (18) del paper
        return quad(lambda y: V(mod, y, M, k) / dV(mod, y, M, k), xe, x)[0]
    objetivo = Nst if callable(Nst) else (lambda xs, xe, pe: Nst)
    xs = brentq(lambda x: N(x) - objetivo(x, xe, pe), xe, 30 * M)
    Nst = N(xs)
    eps, eta = epsV(mod, xs, M, k), d2V(mod, xs, M, k) / V(mod, xs, M, k)
    lam = AS * 24 * np.pi ** 2 * eps / V(mod, xs, M, k)
    ns, r = 1 - 6 * eps + 2 * eta, 16 * eps
    lam38 = 24 * al * np.pi ** 2 * AS / Nst ** 2

    # Convención de los .h, en GeV.
    A = 0.75 * lam * MP ** 4
    MG = M * MP
    if mod == "E":
        phi0 = MG * np.log(1 + k * MP / (np.sqrt(2) * MG))
        V0 = A * (1 - np.exp(-phi0 / MG)) ** k
    else:
        phi0 = MG / 2 * np.arcsinh(np.sqrt(2) * k * MP / MG)
        V0 = A * np.tanh(phi0 / MG) ** k
    pi0 = -np.sqrt((np.sqrt(7 / 3) - 1) * V0)
    om = np.sqrt(k * A) * MG ** (-k / 2) * phi0 ** (k / 2 - 1)
    Hinf = np.sqrt(lam * V(mod, xs, M, k) / 3) * MP
    return dict(M_MP=M, phi_end=xe * MP, phi_star=xs * MP, lam=lam, lam_eq38=lam38, ns=ns, r=r, N_star=Nst,
                rho_end=lam * (V(mod, xe, M, k) + pe ** 2 / 2) * MP ** 4,
                A=A, M=MG, phi0=phi0, pi0=pi0, omega_star=om, H_star=Hinf,
                lam_eff=4 * A / MG ** 4)


def chequeo_tabla_II():
    """n_s y r de slow-roll con alpha_att = 1, N_* = 55, contra la Tabla II del paper."""
    tabla = {"E": (0.9643, 0.00371), "T": (0.9636, 0.00390)}
    for mod in "ET":
        d = calcular(mod, 1.0, 55.0)
        print(f"Tabla II, {mod}-model k=4, N_*=55: n_s = {d['ns']:.4f} (paper {tabla[mod][0]}), "
              f"r = {d['r']:.5f} (paper {tabla[mod][1]})")


def monomial_control(Nst=55.0):
    """Monomial cuártico de control: V = A |phi|^4 con el lambda_eff del T-model (alpha_K = 1).

    No es un modelo del paper (phi^4 puro está excluido por r): sirve para aislar el efecto del
    plateau, porque cerca del mínimo tiene el mismo potencial de programa y, con el mismo q, el mismo
    acople físico g = sqrt(q lambda_eff). Condición inicial: epsilon_V = 1 (Notas.tex «Condiciones iniciales» y «Momento inicial»).
    n_s y r de phi^4 en slow-roll: phi_*^2 = 8(N_* + 1) m_p^2, n_s = 1 - 3/(N_*+1), r = 16/(N_*+1).
    """
    dT = calcular("T", 1.0, 55.8)
    lam_eff = dT["lam_eff"]
    A = lam_eff / 4
    phi0 = K * MP / np.sqrt(2)
    pi0 = -np.sqrt((np.sqrt(7 / 3) - 1) * A * phi0 ** K)
    om = np.sqrt(K * A) * phi0
    return dict(A=A, lam_eff=lam_eff, phi0=phi0, pi0=pi0, omega_star=om,
                ns=1 - 3 / (Nst + 1), r=16 / (Nst + 1))


IN_CUARTICOS = {("T", 1.0): "gentmodel_cuartico.in", ("E", 1.0): "genemodel_cuartico.in",
                ("E", 5.0): "genemodel_cuartico_aK5.in"}


def chequear_in():
    """Compara M, A, phi_0 y phi'_0 de los .in cuárticos con lo que calcula este script."""
    import os
    import re
    pars = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "models", "parameter-files")
    peor = 0.0
    revisados = []
    for mod, al, Nst in CASOS:
        ruta = os.path.join(pars, IN_CUARTICOS[(mod, al)])
        if not os.path.exists(ruta):  # el .in está en la rama de su modelo
            continue
        revisados.append(IN_CUARTICOS[(mod, al)])
        d = calcular(mod, al, Nst)
        txt = open(ruta, encoding="utf-8").read()
        leer = lambda k: float(re.search(rf"^{k}\s*=\s*(\S+)", txt, re.M).group(1))  # noqa: E731
        for k, v in (("M", d["M"]), ("A", d["A"]), ("initial_amplitudes", d["phi0"]),
                     ("initial_momenta", d["pi0"])):
            peor = max(peor, abs(leer(k) / v - 1))
    ruta = os.path.join(pars, "genmonomial_cuartico.in")
    d = monomial_control()
    txt = open(ruta, encoding="utf-8").read() if os.path.exists(ruta) else ""
    if txt:
        revisados.append("genmonomial_cuartico.in")
    for k, v in ([] if not txt else (("A", d["A"]), ("initial_amplitudes", d["phi0"]), ("initial_momenta", d["pi0"]))):
        peor = max(peor, abs(float(re.search(rf"^{k}\s*=\s*(\S+)", txt, re.M).group(1)) / v - 1))
    if not revisados:
        print(".in cuárticos: ninguno en esta rama (están en t-model, e-model y monomial)")
        return None
    print(f".in cuárticos ({', '.join(revisados)}) vs este script: máxima diferencia relativa = {peor:.1e}")
    assert peor < 1e-4
    return peor


if __name__ == "__main__":
    chequeo_tabla_II()
    for mod, al, Nst in CASOS:
        d = calcular(mod, al, Nst)
        print(f"\n{mod}-model, k = 4, alpha_att = {al:g}, N_* = {Nst}:")
        print(f"  n_s = {d['ns']:.4f}, r = {d['r']:.5f}, lambda = {d['lam']:.4e} "
              f"(Ec. 38: {d['lam_eq38']:.3e}), H_* = {d['H_star']:.3e} GeV")
        print(f"  phi_end = {d['phi_end']:.4e} GeV, phi_* = {d['phi_star']:.4e} GeV")
        print(f"  .in:  n = 4, M = {d['M']:.5e} GeV ({d['M_MP']:.4f} M_P), A = {d['A']:.5e} GeV^4")
        print(f"        initial_amplitudes = {d['phi0']:.5e}, initial_momenta = {d['pi0']:.4e}")
        print(f"        omega_* = {d['omega_star']:.4e} GeV, "
              f"lambda_eff = 4A/M^4 = {d['lam_eff']:.3e}")
    d = monomial_control()
    print(f"\nMonomial de control (lambda_eff del T-model): A = {d['A']:.5e}, phi_0 = {d['phi0']:.5e} GeV, "
          f"phi'_0 = {d['pi0']:.4e} GeV^2, omega_* = {d['omega_star']:.4e} GeV; "
          f"phi^4 con N_* = 55: n_s = {d['ns']:.4f}, r = {d['r']:.3f}")
    chequear_in()
