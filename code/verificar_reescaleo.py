#!/usr/bin/env python3
"""Verificación con sympy de las fórmulas de Notas.tex, Sec. "Análisis de las corridas", y
comparación numérica con la literatura.

Uso: python3 code/verificar_reescaleo.py

Bloques:
  sympy_autosimilar      w = (n-2)/(n+2) (virial) y amplitud phi ~ a^(-6/(n+2)).
  sympy_programa         continuidad y 2da ec. de Friedmann en variables de programa.
  sympy_reescaleo        f_0 y h^2 Omega_GW,0: cadena de factores de escala, eps, y su forma en
                         variables de programa (omega_*/f_*) para cada modelo.
  comparar_literatura    constantes numéricas contra Dufaux et al. 2007 (arXiv:0707.0875, Ecs.
                         45-47) y Figueroa & Torrenti 2017 (arXiv:1707.04533, Ecs. 2.24-2.27).
"""

import os
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reescaleo_gws as rg  # noqa: E402


def es_cero(expr):
    return sp.simplify(expr) == 0


def sympy_autosimilar():
    n, A, phi, a, c = sp.symbols("n A phi a c", positive=True)
    V = A * phi ** n
    # Virial para V ~ phi^n (promedio sobre un período): <phi V'> = <phidot^2>, y phi V' = n V.
    K = sp.Rational(1, 2) * n * V          # <phidot^2/2> = (n/2) <V>
    w = sp.simplify((K - V) / (K + V))
    assert es_cero(w - (n - 2) / (n + 2))
    # rho ~ a^(-3(1+w)) y rho ~ Phi^n  =>  Phi ~ a^(-3(1+w)/n) = a^(-6/(n+2)).
    expo = sp.simplify(-3 * (1 + w) / n)
    assert es_cero(expo + 6 / (n + 2))
    # 1 - 3w = 2(4-n)/(n+2): eps = (a_e/a_RD)^(1-3w) vale 1 para n = 4.
    assert es_cero(1 - 3 * w - 2 * (4 - n) / (n + 2))
    print("autosimilar: w = (n-2)/(n+2), Phi ~ a^(-6/(n+2)), 1-3w = 2(4-n)/(n+2)  [ok]")


def sympy_programa():
    """Ecuaciones de fondo en t~ (d t~ = omega_* a^(-alpha) dt), rho = rho~ f_*^2 omega_*^2."""
    t = sp.symbols("t", positive=True)
    al, om, fs, mp = sp.symbols("alpha omega_s f_s m_p", positive=True)
    a = sp.Function("a")(t)
    rho, p = sp.Function("rho")(t), sp.Function("p")(t)
    # Cosmología en tiempo cósmico:  H^2 = rho/(3 mp^2),  addot/a = -(rho+3p)/(6 mp^2),
    # rhodot = -3H(rho+p).
    H = sp.diff(a, t) / a
    addot = -a * (rho + 3 * p) / (6 * mp ** 2)
    adot2 = a ** 2 * rho / (3 * mp ** 2)
    # d/dt~ = (a^alpha/omega_*) d/dt
    D = lambda f: a ** al / om * sp.diff(f, t)  # noqa: E731
    ap = D(a)
    app = sp.expand(D(ap)).subs(sp.diff(a, t, 2), addot)
    app = app.subs(sp.diff(a, t) ** 2, adot2)
    rt, pt = rho / (fs ** 2 * om ** 2), p / (fs ** 2 * om ** 2)
    # 2da ec. de Friedmann en programa: a'' = a^(2alpha+1) (f_*/m_p)^2 [(2alpha-1) rho~ - 3 p~]/6
    app_notas = a ** (2 * al + 1) * (fs / mp) ** 2 * ((2 * al - 1) * rt - 3 * pt) / 6
    assert es_cero(app - app_notas)
    # 1ra ec. de Friedmann (la que chequea CosmoLattice): a'^2 = a^(2alpha+2) (f_*/m_p)^2 rho~/3
    assert es_cero((ap ** 2).subs(sp.diff(a, t) ** 2, adot2)
                   - a ** (2 * al + 2) * (fs / mp) ** 2 * rt / 3)
    # Continuidad: d rho~/dt~ + 3 (a'/a)(rho~ + p~) = 0, sin depender de alpha.
    rhodot = -3 * H * (rho + p)
    cont = D(rt).subs(sp.diff(rho, t), rhodot) + 3 * ap / a * (rt + pt)
    assert es_cero(cont)
    print("programa: a'' = a^(2a+1)(f*/mp)^2[(2a-1)rho~-3p~]/6, a'^2 = a^(2a+2)(f*/mp)^2 rho~/3,"
          " rho~' = -3(a'/a)(rho~+p~)  [ok]")


def sympy_reescaleo():
    k, a_e, a_RD, a_0, rho_e, T0, g, gs, g_s0, w = sp.symbols(
        "k a_e a_RD a_0 rho_e T_0 g_* g_s* g_s0 w", positive=True)
    # Entre a_e y a_RD: rho ~ a^(-3(1+w)). En RD: rho = (pi^2/30) g T^4, y entropía conservada:
    # a T g_s^(1/3) = cte  =>  a_RD/a_0 = (g_s0/g_s*)^(1/3) T_0/T_RD.
    rho_RD = rho_e * (a_e / a_RD) ** (3 * (1 + w))
    T_RD = (30 * rho_RD / (sp.pi ** 2 * g)) ** sp.Rational(1, 4)
    a0_sobre_aRD = (gs / g_s0) ** sp.Rational(1, 3) * T_RD / T0
    ae_sobre_a0 = (a_e / a_RD) / a0_sobre_aRD
    eps = (a_e / a_RD) ** (1 - 3 * w)
    # f_0 = k/(2 pi a_0) con a_i = 1:
    f0 = k / (2 * sp.pi) / a_e * ae_sobre_a0
    f0_notas = (k / (a_e * rho_e ** sp.Rational(1, 4)) * eps ** sp.Rational(1, 4)
                * T0 * (sp.pi ** 2 * g / 30) ** sp.Rational(1, 4) * (g_s0 / gs) ** sp.Rational(1, 3)
                / (2 * sp.pi))
    assert es_cero(sp.powsimp(sp.expand_power_base(f0 / f0_notas, force=True), force=True) - 1)
    # Amplitud: rho_GW ~ a^-4 => Omega_GW,0 = Omega_GW,e rho_e (a_e/a_0)^4 / rho_c0.
    # rho_c0 = rho_gamma,0/Omega_gamma,0 con rho_gamma,0 = (pi^2/30) 2 T_0^4.
    Og = sp.symbols("Omega_gamma", positive=True)
    rho_c0 = (sp.pi ** 2 / 30) * 2 * T0 ** 4 / Og
    amp = rho_e * ae_sobre_a0 ** 4 / rho_c0
    amp_notas = Og * g / 2 * (g_s0 / gs) ** sp.Rational(4, 3) * eps
    assert es_cero(sp.powsimp(sp.expand_power_base(amp / amp_notas, force=True), force=True) - 1)
    # Con g = g_s (y g_0 en lugar de g_s0) se recupera Omega_rad (g_*/g_0)^(-1/3) de Dufaux et al.
    g0 = sp.symbols("g_0", positive=True)
    Orad = Og * g0 / 2
    assert es_cero(sp.simplify((amp_notas.subs({gs: g, g_s0: g0}) / eps) - Orad * (g / g0) ** sp.Rational(-1, 3)))
    # Variables de programa: k = k~ omega_*, rho_e = rho~_e f_*^2 omega_*^2
    kt, rt, om, fs = sp.symbols("ktilde rhotilde omega_s f_s", positive=True)
    comb = sp.simplify((k / rho_e ** sp.Rational(1, 4)).subs({k: kt * om, rho_e: rt * fs ** 2 * om ** 2}))
    assert es_cero(comb - kt * sp.sqrt(om / fs) / rt ** sp.Rational(1, 4))
    # omega_*/f_* para cada modelo (Notas.tex): con n = 4 no depende de phi_0.
    n, A, M, phi0, lam = sp.symbols("n A M phi_0 lambda", positive=True)
    om_mono = sp.sqrt(n * A) * phi0 ** (n / 2 - 1)
    om_TE = sp.sqrt(n * A) * M ** (-n / 2) * phi0 ** (n / 2 - 1)
    assert es_cero((om_mono / phi0).subs(n, 4) - 2 * sp.sqrt(A))
    assert es_cero((om_TE / phi0).subs(n, 4) - sp.sqrt(4 * A / M ** 4))
    assert es_cero((sp.sqrt(lam) * phi0 / phi0) - (om_mono / phi0).subs({n: 4, A: lam / 4}))
    print("reescaleo: f_0 y h^2 Omega_GW,0 de la cadena de factores de escala; forma con "
          "k~ (omega_*/f_*)^(1/2)/rho~^(1/4); omega_*/f_* = sqrt(lambda_eff) para n = 4  [ok]")


def comparar_literatura():
    res = {}
    # Dufaux et al. 2007, Ec. (45)-(47): f = k/(a_j rho_j^(1/4)) * 4e10 Hz y, para lambda phi^4,
    # Omega h^2 = 9.3e-6 x (...), con Omega_rad h^2 = 4.3e-5 y g_*/g_0 = 100, sin distinguir g de g_s.
    g_d = 100 * rg.G0
    Cf_d = (rg.T0 * (np.pi ** 2 * g_d / 30) ** 0.25 * (rg.G0 / g_d) ** (1 / 3) * rg.GEV_A_HZ)
    Orad = rg.omega_gamma_h2() * rg.G0 / 2
    CO_d = Orad * (g_d / rg.G0) ** (-1 / 3)
    res["C_f_dufaux"] = Cf_d
    res["C_Omega_dufaux"] = CO_d
    res["Omega_rad_h2"] = Orad
    print(f"Dufaux et al. (g_*/g_0 = 100, g = g_s): C_f = {Cf_d:.3e} Hz (paper: 4e10), "
          f"C_Omega = {CO_d:.3e} (paper: 9.3e-6 con Omega_rad h^2 = 4.3e-5; acá {Orad:.3e})")
    print(f"    reescalando nuestro Omega_rad h^2 al del paper: {CO_d * 4.3e-5 / Orad:.3e}")
    # Figueroa & Torrenti 2017, Ec. (2.26): f = eps^(1/4) k/rho_i^(1/4) x 8e9 Hz con
    # rho_0 ~ 2e-15 eV^4 y (g_s,o/g_s,RD)^(1/3) (g_o/g_RD)^(-1/4) ~ 1.
    rho0_eV = 2e-15
    f_FT = rho0_eV ** 0.25 * 1e-9 * rg.GEV_A_HZ   # rho_0^(1/4)/(2 pi), en Hz
    res["C_f_FT_sin_g"] = f_FT
    print(f"Figueroa & Torrenti (Ec. 2.26): rho_0^(1/4)/(2pi) = {f_FT:.3e} Hz; el paper escribe 8e9 Hz,"
          f" un factor {f_FT / 8e9:.2f} ~ 2pi menor")
    # Ec. (2.27): h^2 Omega_rad (g_o/g_RD)^(1/3) -> con g_RD = 100:
    CO_FT = Orad * (rg.G0 / 100) ** (1 / 3)
    res["C_Omega_FT_g100"] = CO_FT
    res["cociente_FT"] = rg.C_Omega(100) / CO_FT
    print(f"Figueroa & Torrenti (Ec. 2.27), g_RD = 100: h^2 Omega_rad (g_o/g_RD)^(1/3) = {CO_FT:.3e};"
          f" nuestro C_Omega(100) = {rg.C_Omega(100):.3e}; el cociente {rg.C_Omega(100) / CO_FT:.3f} es "
          f"(g_s0/g_0)^(4/3) = {(rg.GS0 / rg.G0) ** (4 / 3):.3f} (ellos toman g_s = g hoy)")
    res["C_f_SM"] = rg.C_f()
    res["C_Omega_SM"] = rg.C_Omega()
    res["omega_gamma_h2"] = rg.omega_gamma_h2()
    print(f"Nuestro, Modelo Estándar (g_* = g_s* = 106.75): C_f = {rg.C_f():.4e} Hz, "
          f"C_Omega = {rg.C_Omega():.4e}, h^2 Omega_gamma = {rg.omega_gamma_h2():.4e}")
    return res


if __name__ == "__main__":
    sympy_autosimilar()
    sympy_programa()
    sympy_reescaleo()
    comparar_literatura()
