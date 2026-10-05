#!/usr/bin/env python3
"""Parámetros de los alpha-attractors cuadráticos (k = 2) favorecidos por Ellis, Garcia, Olive y
Verner, PRD 113, 063571 (2026), traducidos a la convención de gentmodel/genemodel/genmonomial
(Notas/Notas.tex, «Parámetros cosmológicos para modelos cuadráticos»).

Uso: python3 code/parametros_cuadraticos.py [--TRH 1e10]

A diferencia de k = 4, para k = 2 el inflatón oscila como materia (w = 0) hasta el reheating, así
que N_* depende de T_RH (Ec. 27 del paper, con w_int = 0):

    N_* = 61.49 - ln(g_RH)/12 + (1/4) ln(V_*^2/(M_P^4 rho_end)) + (1/12) ln(rho_RH/rho_end),

con rho_RH = g_RH pi^2 T_RH^4/30 y rho_end la densidad en epsilon_H = 1 (fondo exacto). Como
lambda (y por lo tanto V_* y rho_end) depende de phi_*, N_* se resuelve en forma autoconsistente.
Se chequea contra las Ecs. B9-B10 del paper (alpha_K = 1) y contra «n_s ~ 0.961 para N_* = 50,
T_RH ~ 5e7 GeV» (Sec. VII.A).

Casos (los mismos alpha_K que en el caso cuártico): E y T con alpha_K = 1 y E con alpha_K = 5, con
T_RH = 1e10 GeV por defecto: es la cota de gravitinos, dentro de la ventana que Planck 2018 permite
para E (220-1e10 GeV) y T (2e4-1e10 GeV) con alpha_K = 1, y la que da el n_s más alto de esa ventana
(Fig. 2 del paper). Se agrega el monomial m^2 phi^2 / 2 de control (excluido por r) con la masa del
T-model.
"""

import argparse

import numpy as np

import parametros_cuarticos as pc

K = 2
G_RH = 915 / 4  # MSSM, como el paper para N_* (Ec. 30)
CASOS = [("E", 1.0), ("T", 1.0), ("E", 5.0)]


def N_estrella(mod, al, TRH, g_RH=G_RH):
    """Función objetivo N_*(phi_*, phi_end, phi'_end) de la Ec. 27 con w_int = 0 (unidades M_P)."""
    M = pc.Mesc(mod, al)
    rho_RH = g_RH * np.pi ** 2 / 30 * (TRH / pc.MP) ** 4

    def objetivo(xs, xe, pe):
        lam = pc.AS * 24 * np.pi ** 2 * pc.epsV(mod, xs, M, K) / pc.V(mod, xs, M, K)
        Vst = lam * pc.V(mod, xs, M, K)
        rho_end = lam * (pc.V(mod, xe, M, K) + pe ** 2 / 2)
        return (61.49 - np.log(g_RH) / 12 + 0.25 * np.log(Vst ** 2 / rho_end)
                + np.log(rho_RH / rho_end) / 12)
    return objetivo


def calcular(mod, al, TRH):
    return pc.calcular(mod, al, N_estrella(mod, al, TRH), k=K)


def chequeo_B9_B10():
    """N_* contra las aproximaciones analíticas B9 (E) y B10 (T) del paper, con alpha_K = 1."""
    for TRH in (1e2, 5e7, 1e10):
        for mod, c in (("E", 59.55), ("T", 59.67)):
            d = calcular(mod, 1.0, TRH)
            Nb = c
            for _ in range(20):  # B9/B10 son implícitas en N_*
                Nb = c - np.log(Nb) / 3 + np.log(TRH / pc.MP) / 3
            print(f"  {mod}, T_RH = {TRH:.0e} GeV: N_* = {d['N_star']:.2f} (Ec. B{9 if mod == 'E' else 10}: "
                  f"{Nb:.2f}), n_s = {d['ns']:.4f}")


def monomial_control(TRH):
    """V = A phi^2 con la masa del T-model alpha_K = 1 cerca del mínimo (m^2 = 2 A_T/M_T^2).
    phi_0 = sqrt(2) m_p (epsilon_V = 1). n_s y r de m^2 phi^2 en slow-roll, con el N_* del T-model."""
    dT = calcular("T", 1.0, TRH)
    m2 = 2 * dT["A"] / dT["M"] ** 2
    A = m2 / 2
    phi0 = K * pc.MP / np.sqrt(2)
    pi0 = -np.sqrt((np.sqrt(7 / 3) - 1) * A * phi0 ** K)
    Nst = dT["N_star"]
    return dict(A=A, m=np.sqrt(m2), phi0=phi0, pi0=pi0, omega_star=np.sqrt(m2), N_star=Nst,
                ns=1 - 2 / (Nst + 0.5), r=8 / (Nst + 0.5))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--TRH", type=float, default=1e10)
    TRH = ap.parse_args().TRH
    print("Chequeo de N_*(T_RH) contra el paper:")
    chequeo_B9_B10()
    for mod, al in CASOS:
        d = calcular(mod, al, TRH)
        print(f"\n{mod}-model, k = 2, alpha_K = {al:g}, T_RH = {TRH:.0e} GeV: N_* = {d['N_star']:.2f}")
        print(f"  n_s = {d['ns']:.4f}, r = {d['r']:.5f}, lambda = {d['lam']:.4e}, H_* = {d['H_star']:.3e} GeV, "
              f"rho_end^(1/4) = {d['rho_end'] ** 0.25:.3e} GeV")
        print(f"  .in:  n = 2, M = {d['M']:.5e} GeV ({d['M_MP']:.4f} M_P), A = {d['A']:.5e} GeV^4")
        print(f"        initial_amplitudes = {d['phi0']:.5e}, initial_momenta = {d['pi0']:.4e}")
        print(f"        omega_* = m = {d['omega_star']:.4e} GeV, phi_0/M = {d['phi0'] / d['M']:.3f}")
    d = monomial_control(TRH)
    print(f"\nMonomial de control (m del T-model): A = {d['A']:.5e} GeV^2, m = {d['m']:.4e} GeV, "
          f"phi_0 = {d['phi0']:.5e} GeV, phi'_0 = {d['pi0']:.4e} GeV^2; con N_* = {d['N_star']:.1f}: "
          f"n_s = {d['ns']:.4f}, r = {d['r']:.3f}")
