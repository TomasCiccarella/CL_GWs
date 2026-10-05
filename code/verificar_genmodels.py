#!/usr/bin/env python3
"""Verificación de los modelos genmonomial, gentmodel y genemodel (ver Notas/Notas.tex).

Uso: python3 code/verificar_genmodels.py

Tres bloques, cada uno independiente de CosmoLattice:
  chequear_derivadas      compara con sympy las derivadas del potencial de programa que
                          están escritas a mano en los .h contra la derivada simbólica.
  condiciones_iniciales   recalcula phi_0 (epsilon_V = 1), el momento inicial y omega_*
                          a partir de los parámetros de cada .in, y los compara con lo que
                          el .in tiene escrito. También integra el fondo exacto para ver
                          cuánto se aparta epsilon_V = 1 del final exacto (epsilon_H = 1).
  cotas_estabilidad       factor de escala a partir del cual el integrador (VV2) deja de
                          ser estable para n > 4 (Notas.tex, subsección de estabilidad).
"""

import os
import re

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PARS = os.path.join(ROOT, "models", "parameter-files")
MP = 2.435e18  # masa de Planck reducida en GeV, la misma que usan los .in de CosmoLattice
MODELOS = ("genmonomial", "gentmodel", "genemodel")


def leer_in(modelo):
    """Lee el .in de un modelo y devuelve {clave: [floats]} (sólo las claves numéricas)."""
    out = {}
    for linea in open(os.path.join(PARS, modelo + ".in"), encoding="utf-8"):
        linea = linea.split("#")[0].strip()
        m = re.match(r"(\w+)\s*=\s*(.+)", linea)
        if not m:
            continue
        try:
            out[m.group(1)] = [float(v) for v in m.group(2).split()]
        except ValueError:
            pass
    return out


# --- Potenciales físicos y sus formas cerradas ------------------------------------------

def V(modelo, f, p):
    n, A = p["n"], p["A"]
    if modelo == "genmonomial":
        return A * abs(f) ** n
    if modelo == "gentmodel":
        return A * abs(np.tanh(f / p["M"])) ** n
    return A * abs(1 - np.exp(-f / p["M"])) ** n


def dV(modelo, f, p):
    n, A = p["n"], p["A"]
    if modelo == "genmonomial":
        return A * n * np.sign(f) * abs(f) ** (n - 1)
    if modelo == "gentmodel":
        t = np.tanh(f / p["M"])
        return A * n / p["M"] * np.sign(t) * abs(t) ** (n - 1) / np.cosh(f / p["M"]) ** 2
    u = 1 - np.exp(-f / p["M"])
    return A * n / p["M"] * np.sign(u) * abs(u) ** (n - 1) * np.exp(-f / p["M"])


def phi0_epsilon1(modelo, p):
    """phi_0 tal que epsilon_V = 1 (Notas.tex, subsecciones 'Condiciones iniciales')."""
    n = p["n"]
    if modelo == "genmonomial":
        return n * MP / np.sqrt(2)
    if modelo == "gentmodel":
        return p["M"] / 2 * np.arcsinh(np.sqrt(2) * n * MP / p["M"])
    return p["M"] * np.log(1 + n * MP / (np.sqrt(2) * p["M"]))


def omega_star(modelo, p, f0):
    n, A = p["n"], p["A"]
    if modelo == "genmonomial":
        return np.sqrt(n * A) * f0 ** (n / 2 - 1)
    return np.sqrt(n * A) * p["M"] ** (-n / 2) * f0 ** (n / 2 - 1)


def momento_inicial(V0):
    """phi'_0 de slow-roll autoconsistente en epsilon_V = 1 (Notas.tex, 'Momento inicial').

    phi' = -V'/(3H) con H^2 = (V + phi'^2/2)/(3 m_p^2), y m_p^2 V'^2 = 2 V^2 en epsilon_V = 1,
    da phi'^2 = (sqrt(7/3) - 1) V(phi_0), igual para los tres modelos.
    """
    return -np.sqrt((np.sqrt(7 / 3) - 1) * V0)


def fin_exacto(modelo, p, om):
    """Integra el fondo homogéneo desde lo profundo de la inflación hasta epsilon_H = 1."""
    def rhs(t, y):
        f, fd = y
        H = np.sqrt((V(modelo, f, p) + fd ** 2 / 2) / 3) / MP
        return [fd, -3 * H * fd - dV(modelo, f, p)]

    def eps_H_menos_1(t, y):  # epsilon_H = (3/2) phi'^2 / rho
        return 1.5 * y[1] ** 2 / (V(modelo, y[0], p) + y[1] ** 2 / 2) - 1
    eps_H_menos_1.terminal = True

    f0 = phi0_epsilon1(modelo, p)
    fi = 3 * f0 if modelo == "genmonomial" else f0 + 3 * p["M"]
    Hi = np.sqrt(V(modelo, fi, p) / 3) / MP
    sol = solve_ivp(rhs, [0, 1e6 / om], [fi, -dV(modelo, fi, p) / (3 * Hi)],
                    events=eps_H_menos_1, rtol=1e-10, atol=1e-30, max_step=1 / om)
    return float(sol.y_events[0][0][0])


# --- Bloques de verificación ------------------------------------------------------------

def chequear_derivadas():
    """Derivadas del potencial de programa de los .h contra la derivada simbólica de sympy."""
    x, n, Mt = sp.symbols("x n Mt", positive=True)
    th, e = sp.tanh(x / Mt), sp.exp(-x / Mt)
    casos = {  # (V~, V~' del .h, V~'' del .h), rama x > 0
        "genmonomial": (x ** n / n, x ** (n - 1), (n - 1) * x ** (n - 2)),
        "gentmodel": ((Mt * th) ** n / n,
                      Mt ** (n - 1) * th ** (n - 1) / sp.cosh(x / Mt) ** 2,
                      Mt ** (n - 2) * th ** (n - 2) / sp.cosh(x / Mt) ** 2
                      * ((n - 1) / sp.cosh(x / Mt) ** 2 - 2 * th ** 2)),
        "genemodel": ((Mt * (1 - e)) ** n / n,
                      Mt ** (n - 1) * e * (1 - e) ** (n - 1),
                      Mt ** (n - 2) * e * (1 - e) ** (n - 2) * (n * e - 1)),
    }
    peor = 0.0
    for modelo, (Vt, d1, d2) in casos.items():
        for nv in (1.5, 2, 2.5, 3, 4, 6):
            for xv in (0.3, 1.0, 2.5):
                s = {n: nv, Mt: 0.7, x: xv}
                peor = max(peor,
                           abs(float((sp.diff(Vt, x) - d1).subs(s))),
                           abs(float((sp.diff(Vt, x, 2) - d2).subs(s))))
    print(f"derivadas de los .h vs sympy: máxima diferencia = {peor:.1e}")
    assert peor < 1e-12
    return peor


def condiciones_iniciales():
    """phi_0, phi'_0 y omega_* de cada .in, recalculados y comparados con lo escrito."""
    res = {}
    for modelo in MODELOS:
        d = leer_in(modelo)
        p = {k: d[k][0] for k in ("n", "A", "M", "q") if k in d}
        f0 = phi0_epsilon1(modelo, p)
        pi0 = momento_inicial(V(modelo, f0, p))
        om = omega_star(modelo, p, f0)
        eps = 0.5 * MP ** 2 * (dV(modelo, f0, p) / V(modelo, f0, p)) ** 2
        fe = fin_exacto(modelo, p, om)
        res[modelo] = dict(phi0=f0, pi0=pi0, omega_star=om, phi_fin_exacto=fe)
        print(f"{modelo}: phi_0 = {f0:.5e} GeV (en .in {d['initial_amplitudes'][0]:.5e}), "
              f"eps_V(phi_0) = {eps:.6f}")
        print(f"    phi'_0 = {pi0:.4e} GeV^2 (en .in {d['initial_momenta'][0]:.4e}), "
              f"omega_* = {om:.4e} GeV, epsilon_H = 1 exacto en phi = {fe:.4e} GeV")
    return res


def cotas_estabilidad(modelo="genmonomial"):
    """Factor de escala máximo estable para n > 4 (Notas.tex, 'Estabilidad numérica').

    Gradiente: k_max dt a^((2n-8)/(n+2)) < 2, con k_max = 2 sqrt(3)/dx, dx = 2 pi/(N kIR).
    Hija:      sqrt(q) dt a^(3(n-4)/(n+2)) < 2, tomando la envolvente phi~ ~ a^(-6/(n+2)).
    """
    d = leer_in(modelo)
    n, q, dt, N, kIR = (d[k][0] for k in ("n", "q", "dt", "N", "kIR"))
    if n <= 4:
        print(f"{modelo}: n = {n} <= 4, sin cota creciente con a")
        return None
    kmax = 2 * np.sqrt(3) * N * kIR / (2 * np.pi)
    a_grad = (2 / (kmax * dt)) ** ((n + 2) / (2 * n - 8))
    a_hija = (2 / (np.sqrt(q) * dt)) ** ((n + 2) / (3 * (n - 4)))
    print(f"{modelo}: n = {n}, q = {q}, dt = {dt}, N = {N:.0f}, kIR = {kIR}: "
          f"a_max(gradiente) = {a_grad:.0f}, a_max(hija) ~ {a_hija:.0f}")
    return dict(a_max_gradiente=a_grad, a_max_hija=a_hija)


if __name__ == "__main__":
    chequear_derivadas()
    condiciones_iniciales()
    cotas_estabilidad()
