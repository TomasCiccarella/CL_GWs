"""Reescaleo a hoy de los espectros de CosmoLattice (ver Notas/Notas.tex, Sec. "Reescaleo a hoy").

Lo usan code/analisis_cosmolattice.py y code/verificar_reescaleo.py.

Notación (la de las notas):
    k~, rho~, t~      variables de programa: k = k~ omega_*, rho = rho~ f_*^2 omega_*^2
    a_e, rho_e        factor de escala y densidad total al final de la simulación (a_i = 1)
    n                 potencia del potencial cerca del mínimo, V ~ |phi|^n
    w = (n-2)/(n+2)   ecuación de estado de las oscilaciones autosimilares, que suponemos vale
                      desde a_e hasta a_RD (comienzo de la dominación por radiación)
    eps = (a_e/a_RD)^(1-3w) = (a_e/a_RD)^(2(4-n)/(n+2))

Resultados:
    f_0 = k~ (omega_*/f_*)^(1/2) / (a_e rho~_e^(1/4)) * eps^(1/4) * C_f(g_*, g_s*)
    h^2 Omega_GW,0 = C_Omega(g_*, g_s*) * eps * Omega_GW,e(k)
con
    C_f     = T_0 (pi^2 g_*/30)^(1/4) (g_s0/g_s*)^(1/3) / (2 pi), en Hz
    C_Omega = h^2 Omega_gamma,0 (g_*/2) (g_s0/g_s*)^(4/3)
y Omega_GW,e(k) = (1/rho_tot) d rho_GW/d log k, que es lo que escribe CosmoLattice en
spectra_energy_gws.txt (nota técnica de GWs de CosmoLattice, Sec. 4).
"""

import numpy as np

# --- Constantes -------------------------------------------------------------------------
MP = 2.435e18                    # GeV, masa de Planck reducida (la de CosmoLattice y los .in)
HBAR = 6.582119569e-25           # GeV s
KB = 8.617333262e-14             # GeV / K
T0 = 2.7255 * KB                 # GeV, temperatura del CMB hoy (Fixsen 2009)
H100 = 100 / 3.0856775814913673e19 * HBAR   # GeV: 100 km/s/Mpc en unidades naturales
G0 = 2 + 7 / 8 * 6 * (4 / 11) ** (4 / 3)    # grados de libertad relativistas hoy (3.363)
GS0 = 2 + 7 / 8 * 6 * (4 / 11)              # grados de libertad de entropía hoy (3.909)
G_SM = 106.75                    # Modelo Estándar a T >> 100 GeV
GEV_A_HZ = 1 / (2 * np.pi * HBAR)  # f [Hz] = E [GeV] * GEV_A_HZ  (f = omega/2pi, omega = E/hbar)


def omega_gamma_h2():
    """h^2 Omega_gamma,0 = (pi^2/15) T_0^4 / (3 m_p^2 H_100^2)."""
    return np.pi ** 2 / 15 * T0 ** 4 / (3 * MP ** 2 * H100 ** 2)


def w_autosimilar(n):
    """Ecuación de estado promedio de un condensado oscilando en V ~ |phi|^n (virial)."""
    return (n - 2) / (n + 2)


def epsilon(n, a_e_sobre_a_RD=1.0):
    """eps = (a_e/a_RD)^(1-3w), con w = (n-2)/(n+2) entre el final de la simulación y a_RD.

    Para n = 4 vale exactamente 1, sea cual sea a_RD.
    """
    return a_e_sobre_a_RD ** (1 - 3 * w_autosimilar(n))


def C_f(g_star=G_SM, g_s_star=None):
    """Constante de la frecuencia de hoy, en Hz por GeV^(0) (ver docstring del módulo)."""
    g_s_star = g_star if g_s_star is None else g_s_star
    return T0 * (np.pi ** 2 * g_star / 30) ** 0.25 * (GS0 / g_s_star) ** (1 / 3) * GEV_A_HZ


def C_Omega(g_star=G_SM, g_s_star=None):
    """Constante de la amplitud de hoy: h^2 Omega_gamma,0 (g_*/2) (g_s0/g_s*)^(4/3)."""
    g_s_star = g_star if g_s_star is None else g_s_star
    return omega_gamma_h2() * g_star / 2 * (GS0 / g_s_star) ** (4 / 3)


def frecuencia_hoy(k_prog, a_e, rho_prog_e, omega_star, f_star, n, a_e_sobre_a_RD=1.0,
                   g_star=G_SM, g_s_star=None):
    """Frecuencia de hoy [Hz] de un modo comóvil k~ (en unidades de programa, con a_i = 1).

    f_0 = k~ (omega_*/f_*)^(1/2) / (a_e rho~_e^(1/4)) * eps^(1/4) * C_f
    """
    return (np.asarray(k_prog) * np.sqrt(omega_star / f_star) / (a_e * rho_prog_e ** 0.25)
            * epsilon(n, a_e_sobre_a_RD) ** 0.25 * C_f(g_star, g_s_star))


def omega_gw_hoy(omega_gw_e, n, a_e_sobre_a_RD=1.0, g_star=G_SM, g_s_star=None):
    """h^2 Omega_GW hoy a partir del espectro de CosmoLattice al final de la simulación."""
    return C_Omega(g_star, g_s_star) * epsilon(n, a_e_sobre_a_RD) * np.asarray(omega_gw_e)


# --- Variables de programa de cada modelo (Notas.tex, subsecciones 'Variables computacionales')

def variables_programa(modelo, p):
    """Devuelve (f_*, omega_*, alpha, n) a partir de los parámetros del .in / .infos.

    p: dict con 'initial_amplitudes' (lista) y los parámetros del modelo.
    """
    f_star = p["initial_amplitudes"][0]
    if modelo == "lphi4":
        return f_star, np.sqrt(p["lambda"]) * f_star, 1.0, 4.0
    n, A = p["n"], p["A"]
    alpha = 3 * (n - 2) / (n + 2)
    if modelo == "genmonomial":
        return f_star, np.sqrt(n * A) * f_star ** (n / 2 - 1), alpha, n
    if modelo in ("gentmodel", "genemodel"):
        return f_star, np.sqrt(n * A) * p["M"] ** (-n / 2) * f_star ** (n / 2 - 1), alpha, n
    raise ValueError(f"modelo desconocido: {modelo} (pasá --fstar --omegastar --alpha --n)")


def omega_sobre_f(modelo, p):
    """omega_*/f_*, la combinación que aparece en f_0. Para n = 4 no depende de phi_0:
    sqrt(lambda) (lphi4), sqrt(4A) (genmonomial), sqrt(4A)/M^2 = sqrt(lambda_eff) (T/E)."""
    f_star, om, _, _ = variables_programa(modelo, p)
    return om / f_star
