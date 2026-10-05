#!/usr/bin/env python3
"""Análisis de una corrida de CosmoLattice (fondo, campos, espectros, errores, GWs, energías).

Uso:
    python3 code/analisis_cosmolattice.py DIRECTORIO_DE_LA_CORRIDA [opciones]

El directorio es el 'outputfile' de la corrida: tiene que contener el <modelo>.infos que escribe
CosmoLattice y los average_*.txt / spectra_*.txt. Los gráficos y un resumen.txt se guardan en
DIRECTORIO/analisis/ (o en --out).

Opciones:
    --gstar G          g_* al comienzo de la dominación por radiación (106.75, Modelo Estándar)
    --gsstar G         g_s* (por defecto = g_*)
    --a_e_sobre_a_RD X a_e/a_RD: cuánto se expande el universo entre el final de la simulación y la
                       dominación por radiación, con w = (n-2)/(n+2). Irrelevante para n = 4.
    --TRH T            en vez de lo anterior: temperatura de reheating en GeV; a_e/a_RD sale de
                       rho_e y de rho_RD = (pi^2 g_*/30) T^4 (para n = 2, Notas.tex «Reescaleo a hoy»)
    --w_post W         ecuación de estado entre el final de la simulación y a_RD (por defecto la
                       autosimilar (n-2)/(n+2)); por ejemplo, el <w> medido al final de la corrida
    --modelo M         forzar el modelo (por defecto, el nombre del .infos)
    --fstar, --omegastar, --alpha, --n
                       forzar las variables de programa (para modelos que el script no conoce)

Qué se grafica (detalles y fórmulas en Notas/Notas.tex, Sec. "Análisis de las corridas"):
    01_fondo.png             a(t~), H físico, ecuación de estado w, y rho a^{3(1+w_n)}
    02_campos.png            <phi~> reescalado por a^{6/(n+2)} y rms de cada campo
    03_espectros_campos.png  Delta_phi~(k~) a cada tiempo, con la frecuencia de hoy arriba
    04_errores.png           Friedmann, continuidad, a'' (2da ec. de Friedmann), cola UV, GWs
    05_gws.png               Omega_GW(k~, t~) y h^2 Omega_GW(f) hoy
    06_energias.png          cada término de la energía (absoluto y fracción del total)
"""

import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reescaleo_gws as rg  # noqa: E402

# Paleta categórica de referencia (orden fijo) y rampa secuencial de un solo tono para el tiempo.
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
TINTA, TINTA2 = "#0b0b0b", "#52514e"
plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 150, "font.size": 10, "axes.titlesize": 11,
    "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2, "ytick.color": TINTA2,
    "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6, "lines.linewidth": 1.6,
    "legend.frameon": False, "axes.spines.top": False, "axes.spines.right": False,
})


def colores_tiempo(m):
    """m colores de una rampa azul (claro = temprano, oscuro = tardío)."""
    return plt.cm.Blues(np.linspace(0.35, 1.0, max(m, 1)))


# --- Lectura ----------------------------------------------------------------------------

def leer_infos(d):
    """Parámetros de la corrida desde <modelo>.infos. Devuelve (modelo, dict)."""
    infos = glob.glob(os.path.join(d, "*.infos"))
    if not infos:
        sys.exit(f"No hay *.infos en {d}: ¿es el directorio de salida de CosmoLattice?")
    p = {}
    for linea in open(infos[0], encoding="utf-8"):
        partes = linea.split()
        if len(partes) < 2 or partes[0].endswith(":"):
            continue
        vals = []
        for v in partes[1:]:
            try:
                vals.append(float(v))
            except ValueError:
                vals = None
                break
        if vals:
            p[partes[0]] = vals if len(vals) > 1 or partes[0].startswith("initial_") else vals[0]
    return os.path.basename(infos[0])[:-len(".infos")], p


def leer_tabla(ruta):
    """Tabla average_*.txt -> dict {columna: array}."""
    with open(ruta, encoding="utf-8") as fh:
        cab = fh.readline().lstrip("#").split()
    datos = np.atleast_2d(np.loadtxt(ruta, comments="#"))
    return {c: datos[:, i] for i, c in enumerate(cab)}


def leer_espectros(ruta):
    """spectra_*.txt -> lista (un elemento por tiempo) de arrays [k~, valor(es)..., multiplicidad]."""
    bloques, actual = [], []
    for linea in open(ruta, encoding="utf-8"):
        linea = linea.strip()
        if not linea:
            if actual:
                bloques.append(np.array(actual))
                actual = []
        elif linea[0].isdigit() and ":" in linea.split()[0]:
            continue  # cabecera "1:binCentralValue 2:valuesAverage ..."
        else:
            actual.append([float(x) for x in linea.split()])
    if actual:
        bloques.append(np.array(actual))
    return bloques


# --- Cantidades derivadas ---------------------------------------------------------------

def presion(E):
    """p~ a partir de average_energies: escalares p = K - G/3; gauge p = rho/3; potencial -V."""
    p = np.zeros_like(E["t"])
    for c, v in E.items():
        if c.startswith("E^kin_U1") or c.startswith("E^grad_U1") or c.endswith("SU2"):
            p += v / 3
        elif c.startswith("E^kin"):
            p += v
        elif c.startswith("E^grad"):
            p -= v / 3
        elif c.startswith("Vpot"):
            p -= v
    return p


def nombre_energia(c):
    """Etiqueta legible para una columna de average_energies.txt."""
    if c == "rho_GW":
        return "GWs"
    for pre, nom in (("E^kin_scal", "cinética"), ("E^grad_scal", "gradiente")):
        if c.startswith(pre):
            i = c[len(pre):]
            return f"{nom} " + ("$\\phi$" if i == "0" else f"$\\chi_{{{i}}}$")
    if c.startswith("Vpot_term_"):
        return f"potencial, término {c[len('Vpot_term_'):]}"
    return c


def promedio_movil(t, y, ventana):
    """Promedio sobre una ventana de ancho 'ventana' en t (para promediar oscilaciones).

    Donde la ventana no entra completa (los bordes) devuelve NaN, para no inventar valores.
    """
    dt = np.median(np.diff(t))
    m = max(int(round(ventana / dt)), 1)
    out = np.full(len(y), np.nan)
    if m >= len(y):
        return out
    val = np.convolve(y, np.ones(m) / m, mode="valid")
    out[m // 2: m // 2 + len(val)] = val
    return out


def positivo(v):
    """v con los ceros exactos como NaN: en escala log, un cero (p. ej. el error de Friedmann en
    t~ = 0, o una coincidencia exacta en punto flotante) estiraría el eje hasta 1e-300."""
    return np.where(v > 0, v, np.nan)


def cruces_por_cero(t, x):
    """Tiempos (interpolados linealmente) en que x cambia de signo."""
    i = np.where(np.diff(np.sign(x)) != 0)[0]
    return t[i] - x[i] * (t[i + 1] - t[i]) / (x[i + 1] - x[i])


def derivada(t, y):
    """dy/dt con un spline cúbico (error O(dt^4), mejor que diferencias centradas)."""
    from scipy.interpolate import CubicSpline
    return CubicSpline(t, y)(t, 1)


def periodo_programa(n):
    """Período de oscilación de phi~ con V~ = |phi~|^n/n y amplitud 1, en tiempo de programa.

    T = 4 int_0^1 dx / sqrt(2 (1 - x^n)/n) (con a^alpha constante, que es la elección de alpha).
    """
    from scipy.integrate import quad
    return 4 * quad(lambda x: 1 / np.sqrt(2 * (1 - x ** n) / n), 0, 1)[0]


# --- Análisis ---------------------------------------------------------------------------

def analizar(d, args):
    modelo, p = leer_infos(d)
    if args.modelo:
        modelo = args.modelo
    if args.fstar and args.omegastar and args.alpha is not None and args.n:
        fs, om, alpha, n = args.fstar, args.omegastar, args.alpha, args.n
    else:
        fs, om, alpha, n = rg.variables_programa(modelo, p)
    out = args.out or os.path.join(d, "analisis")
    os.makedirs(out, exist_ok=True)
    kIR, N = p["kIR"], p["N"]
    w_n = rg.w_autosimilar(n)

    S = leer_tabla(os.path.join(d, "average_scale_factor.txt"))
    E = leer_tabla(os.path.join(d, "average_energies.txt"))
    C = leer_tabla(os.path.join(d, "average_energy_conservation.txt"))
    campos = sorted(glob.glob(os.path.join(d, "average_scalar_*.txt")))
    F = [leer_tabla(c) for c in campos]
    t, a, ap, Hp = S["t"], S["a"], S["aDot"], S["H"]
    rho, pr = E["E_tot"], presion(E)
    w = pr / rho
    # Período de <phi~>: medido de los cruces por cero (depende de la amplitud), o el teórico
    # para amplitud 1 si no hay suficientes oscilaciones.
    tc = cruces_por_cero(F[0]["t"], F[0]["<phi>"])
    T_osc = 2 * np.mean(np.diff(tc[-5:])) if len(tc) >= 3 else periodo_programa(n)
    w_prom = promedio_movil(t, w, T_osc / 2)  # w oscila con el doble de la frecuencia de phi
    if len(tc) >= 5:  # promedio sobre los dos últimos períodos completos
        w_final = np.mean(w[(t >= tc[-5]) & (t < tc[-1])])
    else:
        w_final = np.mean(w[t > t[-1] - T_osc])
    H_fis = Hp * om * a ** (-alpha)            # GeV

    t_esp = np.atleast_1d(np.loadtxt(os.path.join(d, "average_spectra_times.txt"), comments="#"))
    esp_campos = [leer_espectros(os.path.join(d, f"spectra_scalar_{i}.txt"))
                  for i in range(len(campos))]
    ruta_gw = os.path.join(d, "spectra_energy_gws.txt")
    hay_gw = os.path.exists(ruta_gw)
    esp_gw = leer_espectros(ruta_gw) if hay_gw else []
    EGW = leer_tabla(os.path.join(d, "average_energies_gws.txt")) if hay_gw else None
    m = min(len(t_esp), *(len(e) for e in esp_campos))
    t_esp = t_esp[:m]

    # Final de la simulación: a_e y rho~_e (para el reescaleo a hoy).
    a_e, rho_e = a[-1], rho[-1]
    a_esp = np.interp(t_esp, t, a)
    rho_esp = np.interp(t_esp, t, rho)

    # Ecuación de estado entre el final de la simulación y a_RD: la autosimilar (n-2)/(n+2) o la
    # que se pida con --w_post. reescaleo_gws la toma de n, así que se pasa el n equivalente.
    w_post = w_n if args.w_post is None else args.w_post
    n_post = n if args.w_post is None else 2 * (1 + w_post) / (1 - w_post)

    def a_RD(rho_x):  # a_x/a_RD: fijo (--a_e_sobre_a_RD) o desde T_RH (--TRH), con rho~ -> GeV^4
        if args.TRH:
            return rg.a_sobre_a_RD(rho_x * fs ** 2 * om ** 2, args.TRH, n_post, args.gstar)
        return args.a_e_sobre_a_RD

    def f_hoy(k, a_x=a_e, rho_x=rho_e):
        return rg.frecuencia_hoy(k, a_x, rho_x, om, fs, n_post, a_RD(rho_x), args.gstar, args.gsstar)

    def Om_hoy(Om_e, rho_x=rho_e):
        return rg.omega_gw_hoy(Om_e, n_post, a_RD(rho_x), args.gstar, args.gsstar)
    factor_f = f_hoy(1.0)  # Hz por unidad de k~ al final

    resumen = [f"Corrida: {os.path.abspath(d)}",
               f"modelo = {modelo}, n = {n:g}, alpha = {alpha:g}, w_n = (n-2)/(n+2) = {w_n:.4f}",
               f"f_* = {fs:.5e} GeV, omega_* = {om:.5e} GeV, omega_*/f_* = {om / fs:.4e}",
               f"N = {N:g}, kIR = {kIR:g}, k~_max = {kIR * N * np.sqrt(3) / 2:.3g}, "
               f"t~ final = {t[-1]:g}, a_e = {a_e:.5g}, rho~_e = {rho_e:.5e}",
               f"g_* = {args.gstar}, g_s* = {args.gsstar or args.gstar}, "
               + (f"T_RH = {args.TRH:.3g} GeV (w = {w_post:.3g} hasta entonces), " if args.TRH else "")
               + f"a_e/a_RD = {a_RD(rho_e):.4g} -> eps = {rg.epsilon(n_post, a_RD(rho_e)):.4g}",
               f"rho_e^(1/4) = {(rho_e * fs ** 2 * om ** 2) ** 0.25:.4e} GeV",
               f"f_0 = {factor_f:.4e} Hz x k~   (C_f = {rg.C_f(args.gstar, args.gsstar):.4e} Hz)",
               f"h^2 Omega_GW,0 = {Om_hoy(1.0):.4e}"
               f" x Omega_GW,e",
               f"<w> en los últimos dos períodos = {w_final:.4f} "
               f"(autosimilar: {w_n:.4f}); período de <phi~> en t~: T = {T_osc:.4f}"]

    # ---------------------------------------------------------------- 01 fondo
    fig, ax = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    ax[0, 0].plot(t, a, color=CAT[0])
    ax[0, 0].set(xlabel=r"$\tilde t$", ylabel=r"$a$", title="Factor de escala")
    ax[0, 1].semilogy(t, H_fis, color=CAT[0])
    ax[0, 1].set(xlabel=r"$\tilde t$", ylabel=r"$H$ [GeV]", title=r"Hubble físico, $H = \omega_* a^{-\alpha}\, a'/a$")
    ax[1, 0].plot(t, w, color=CAT[0], lw=0.6, alpha=0.35, label="instantánea")
    ax[1, 0].plot(t, w_prom, color=CAT[0], label=r"promedio en $T/2$")
    ax[1, 0].axhline(w_n, color=TINTA2, ls="--", lw=1, label=rf"$(n-2)/(n+2) = {w_n:.3f}$")
    ax[1, 0].set(xlabel=r"$\tilde t$", ylabel=r"$w = \tilde p/\tilde\rho$", title="Ecuación de estado")
    ax[1, 0].legend(loc="lower right")
    esc = rho * a ** (3 * (1 + w_n))
    ax[1, 1].plot(t, esc / esc[0], color=CAT[0])
    ax[1, 1].set(xlabel=r"$\tilde t$", ylabel=r"$\tilde\rho\, a^{3(1+w_n)} / (\cdot)_0$",
                 title=r"$\rho \propto a^{-3(1+w_n)}$ (constante si es autosimilar)")
    fig.savefig(os.path.join(out, "01_fondo.png"))
    plt.close(fig)

    # ---------------------------------------------------------------- 02 campos
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    ax[0].plot(F[0]["t"], F[0]["<phi>"] * a ** (6 / (n + 2)), color=CAT[0])
    ax[0].set(xlabel=r"$\tilde t$", ylabel=r"$a^{6/(n+2)} \langle\tilde\phi\rangle$",
              title=r"Inflatón homogéneo reescalado por su envolvente autosimilar")
    for i, Fi in enumerate(F):
        nombre = r"\tilde\phi" if i == 0 else rf"\tilde\chi_{{{i}}}" if len(F) > 2 else r"\tilde\chi"
        ax[1].semilogy(Fi["t"], np.maximum(Fi["rms(phi)"], 1e-300), color=CAT[i % 8],
                       label=rf"rms$({nombre})$")
    ax[1].semilogy(F[0]["t"], np.abs(F[0]["<phi>"]) + 1e-300, color=CAT[0], lw=0.5, alpha=0.4,
                   label=r"$|\langle\tilde\phi\rangle|$")
    ax[1].set(xlabel=r"$\tilde t$", ylabel="amplitud (unidades de programa)",
              title="Fluctuaciones de cada campo")
    ax[1].legend()
    fig.savefig(os.path.join(out, "02_campos.png"))
    plt.close(fig)

    # ---------------------------------------------------------------- 03 espectros de campos
    cols = colores_tiempo(m)
    fig, ax = plt.subplots(1, len(esp_campos), figsize=(5.4 * len(esp_campos), 4.4),
                           constrained_layout=True, squeeze=False)
    for i, bloques in enumerate(esp_campos):
        for j in range(m):
            b = bloques[j]
            ok = b[:, 1] > 0
            ax[0, i].loglog(b[ok, 0], b[ok, 1], color=cols[j], lw=1.1)
        nombre = r"\tilde\phi" if i == 0 else r"\tilde\chi"
        ax[0, i].set(xlabel=r"$\tilde k$", ylabel=rf"$\Delta_{{{nombre}}}(\tilde k)$",
                     title=f"Espectro de scalar_{i} ({m} tiempos, claro = temprano)")
        sec = ax[0, i].secondary_xaxis("top", functions=(lambda k: k * factor_f, lambda f: f / factor_f))
        sec.set_xlabel(r"$f_0$ hoy [Hz]")
    sm = plt.cm.ScalarMappable(cmap=plt.cm.Blues, norm=plt.Normalize(t_esp[0], t_esp[-1]))
    fig.colorbar(sm, ax=ax[0, -1], label=r"$\tilde t$")
    fig.savefig(os.path.join(out, "03_espectros_campos.png"))
    plt.close(fig)

    # ---------------------------------------------------------------- 04 errores
    fig, ax = plt.subplots(2, 2, figsize=(10.5, 7.5), constrained_layout=True)
    ax[0, 0].semilogy(C["t"], positivo(np.abs(C["rel_diff_friedmann"])), color=CAT[0])
    ax[0, 0].set(xlabel=r"$\tilde t$", ylabel=r"$|(\mathcal{L}-\mathcal{R})/(\mathcal{L}+\mathcal{R})|$",
                 title=r"1ra ec. de Friedmann: $a'^2$ vs $a^{2\alpha+2}(f_*/m_p)^2\tilde\rho/3$")
    # Continuidad: d rho~/dt~ + 3 (a'/a) (rho~ + p~) = 0 (vale para cualquier alpha).
    drho = derivada(t, rho)
    res_cont = np.abs(drho + 3 * Hp * (rho + pr)) / np.abs(3 * Hp * rho)
    # 2da de Friedmann: a'' = a^{2alpha+1} (f_*/m_p)^2 [(2alpha-1) rho~ - 3 p~]/6.
    app_num = derivada(t, ap)
    pref = a ** (2 * alpha + 1) * (fs / rg.MP) ** 2 / 6
    app_teo = pref * ((2 * alpha - 1) * rho - 3 * pr)
    res_app = np.abs(app_num - app_teo) / (pref * (np.abs((2 * alpha - 1) * rho) + 3 * np.abs(pr)))
    interior = slice(3, -3)  # el spline es peor en los bordes
    t_in, res_cont, res_app = t[interior], res_cont[interior], res_app[interior]
    ax[0, 1].semilogy(t_in, positivo(res_cont), color=CAT[0], label="continuidad")
    ax[0, 1].semilogy(t_in, positivo(res_app), color=CAT[1], label=r"$a''$ (2da Friedmann)")
    ax[0, 1].set(xlabel=r"$\tilde t$", ylabel="residuo relativo",
                 title=f"Derivadas por spline (limitadas por tOutputFreq = {np.median(np.diff(t)):.3g})")
    ax[0, 1].legend()
    colas = {}
    for i, bloques in enumerate(esp_campos):
        colas[f"scalar_{i}"] = [b[-1, 1] / b[:, 1].max() if b[:, 1].max() > 0 else np.nan
                                for b in bloques[:m]]
    if hay_gw:
        colas["GWs"] = [b[-1, 1] / b[:, 1].max() if b[:, 1].max() > 0 else np.nan for b in esp_gw[:m]]
    for i, (nom, c) in enumerate(colas.items()):
        ax[1, 0].semilogy(t_esp, c, "o-", ms=4, color=CAT[i % 8], label=nom)
    ax[1, 0].set(xlabel=r"$\tilde t$", ylabel=r"$\Delta(\tilde k_{\max}) / \max \Delta$",
                 title="Resolución UV: cola del espectro respecto del pico")
    ax[1, 0].legend()
    if hay_gw:
        # Integral del espectro en log k contra rhoGW/rho que escribe CosmoLattice.
        dk = kIR * p.get("deltaKBin", 1.0)
        integ = np.array([np.sum(b[:, 1] * dk / b[:, 0]) for b in esp_gw[:m]])
        ref = np.interp(t_esp, EGW["t"], EGW["rhoGW_over_rho"])
        ok = ref > 0
        ax[1, 1].semilogy(t_esp[ok], positivo(np.abs(integ[ok] / ref[ok] - 1)), "o-", ms=4, color=CAT[0],
                          label=r"$|\sum \Omega_{GW}\,\Delta k/k \,/\, (\rho_{GW}/\rho) - 1|$")
        # Modos de GWs fuera del horizonte: k~/(a^{1-alpha} H~) < 1 no redshiftean como radiación.
        kH = kIR / (a_esp ** (1 - alpha) * np.interp(t_esp, t, Hp))
        ax[1, 1].semilogy(t_esp, 1 / kH, "s--", ms=4, color=CAT[1],
                          label=r"$aH/k_{\rm IR}$ (debe ser $<1$)")
        ax[1, 1].legend()
        chequeo = f"{abs(integ[-1] / ref[-1] - 1):.2e}" if ref[-1] > 0 else "sin GWs todavía"
        resumen.append(f"GWs: |integral/rhoGW_over_rho - 1| final = {chequeo}; "
                       f"k_IR/(aH) final = {kH[-1]:.3g}")
    ax[1, 1].set(xlabel=r"$\tilde t$", title="Chequeos de GWs")
    fig.savefig(os.path.join(out, "04_errores.png"))
    plt.close(fig)
    resumen += [f"Friedmann: max |rel_diff| = {np.max(np.abs(C['rel_diff_friedmann'])):.2e}",
                f"Continuidad: mediana residuo = {np.median(res_cont):.2e}; "
                f"a'': mediana residuo = {np.median(res_app):.2e}"]
    resumen += [f"Cola UV final {nom}: {c[-1]:.2e}" for nom, c in colas.items()]

    # ---------------------------------------------------------------- 05 GWs
    if hay_gw:
        fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)
        for j in range(m):
            b = esp_gw[j]
            ok = b[:, 1] > 0
            if ok.any():
                ax[0].loglog(b[ok, 0], b[ok, 1], color=cols[j], lw=1.1)
                ax[1].loglog(f_hoy(b[ok, 0], a_esp[j], rho_esp[j]),
                             Om_hoy(b[ok, 1], rho_esp[j]),
                             color=cols[j], lw=1.1)
        ax[0].set(xlabel=r"$\tilde k$", ylabel=r"$\Omega_{GW}(\tilde k, \tilde t) = \frac{1}{\rho}\frac{d\rho_{GW}}{d\log k}$",
                  title="Espectro de GWs en la simulación")
        b = esp_gw[m - 1]
        fh = f_hoy(b[:, 0])
        Oh = Om_hoy(b[:, 1])
        jp = int(np.argmax(Oh))
        ax[1].plot(fh[jp], Oh[jp], "o", ms=8, color=CAT[1], mec="white", mew=2,
                   label=rf"pico final: $f = {fh[jp]:.2e}$ Hz, $h^2\Omega = {Oh[jp]:.2e}$")
        ax[1].set(xlabel=r"$f_0$ [Hz]", ylabel=r"$h^2\Omega_{GW,0}(f)$",
                  title=r"Hoy (cada curva, como si la producción terminara en ese $\tilde t$)")
        ax[1].legend(loc="lower right")
        fig.colorbar(sm, ax=ax[1], label=r"$\tilde t$")
        fig.savefig(os.path.join(out, "05_gws.png"))
        plt.close(fig)
        resumen.append(f"Pico de GWs hoy: f = {fh[jp]:.4e} Hz (k~ = {b[jp, 0]:.4g}), "
                       f"h^2 Omega_GW = {Oh[jp]:.4e}; rhoGW/rho final = {EGW['rhoGW_over_rho'][-1]:.4e}")
        np.savetxt(os.path.join(out, "gws_hoy.txt"), np.column_stack([b[:, 0], fh, b[:, 1], Oh]),
                   header="k~  f_0[Hz]  Omega_GW,e  h2_Omega_GW,0  (último tiempo de espectros)")

    # ---------------------------------------------------------------- 06 energías
    # Arriba: valores instantáneos. Abajo: fracciones promediadas sobre un período de <phi~>
    # (cinética y potencial del inflatón se intercambian dentro de cada oscilación).
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    componentes = [c for c in E if c not in ("t", "E_tot")]
    if hay_gw:
        rgw = np.interp(t, EGW["t"], EGW["rhoGW"])
    for i, c in enumerate(componentes + (["rho_GW"] if hay_gw else [])):
        v = np.abs(rgw if c == "rho_GW" else E[c])
        col, nom = CAT[i % 8], nombre_energia(c)
        ax[0].semilogy(t, np.where(v > 0, v, np.nan), color=col, lw=1.1, label=nom)
        fr = promedio_movil(t, v / rho, T_osc)
        ax[1].semilogy(t, np.where(fr > 0, fr, np.nan), color=col, label=nom)
    ax[0].semilogy(t, rho, color=TINTA, lw=1, ls="--", label="total")
    ax[0].set(xlabel=r"$\tilde t$", ylabel=r"$\tilde\rho_i = \rho_i/(f_*^2\omega_*^2)$",
              title="Energía de cada término (instantánea)")
    ax[1].set(xlabel=r"$\tilde t$", ylabel=r"$\langle\rho_i\rangle_T / \rho_{\rm tot}$",
              title=f"Fracción del total, promediada en un período (T = {T_osc:.3g})")
    ax[1].legend(fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.savefig(os.path.join(out, "06_energias.png"))
    plt.close(fig)

    texto = "\n".join(resumen)
    open(os.path.join(out, "resumen.txt"), "w", encoding="utf-8").write(texto + "\n")
    print(texto)
    print(f"\nGráficos en {out}")
    return dict(f_pico=fh[jp] if hay_gw else None, Omega_pico=Oh[jp] if hay_gw else None,
                factor_f=factor_f, a_e=a_e, rho_e=rho_e, friedmann=np.max(np.abs(C["rel_diff_friedmann"])))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("directorio")
    ap.add_argument("--out")
    ap.add_argument("--gstar", type=float, default=rg.G_SM)
    ap.add_argument("--gsstar", type=float, default=None)
    ap.add_argument("--a_e_sobre_a_RD", type=float, default=1.0)
    ap.add_argument("--TRH", type=float, help="temperatura de reheating [GeV]; reemplaza a --a_e_sobre_a_RD")
    ap.add_argument("--w_post", type=float, help="w entre el final de la simulación y a_RD (por defecto (n-2)/(n+2))")
    ap.add_argument("--modelo")
    ap.add_argument("--fstar", type=float)
    ap.add_argument("--omegastar", type=float)
    ap.add_argument("--alpha", type=float)
    ap.add_argument("--n", type=float)
    analizar(ap.parse_args().directorio, ap.parse_args())


if __name__ == "__main__":
    main()
