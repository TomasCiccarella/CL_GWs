# Bibliografía

Bibliografía útil para el proyecto, ordenada por tema. Para cada trabajo: qué aporta y dónde se usa en el repositorio. Las referencias en BibTeX están en [`referencias.bib`](referencias.bib) (las claves son las de la primera columna).

**Los PDFs no se versionan**: el repositorio es público y la mayoría están bajo derechos de autor (libros y versiones publicadas). Van en `Bibliografia/pdf/`, que está en `.gitignore`. Para bajar los que están en arXiv:

```bash
Bibliografia/descargar.sh        # baja a Bibliografia/pdf/ los que faltan (arXiv)
```

Los libros y la nota técnica de CosmoLattice hay que conseguirlos aparte.

## Preheating y resonancia paramétrica

| clave | trabajo | qué aporta | dónde se usa |
|---|---|---|---|
| `Kofman1997` | Kofman *et al.* (1997), «Towards the Theory of Reheating After Inflation», Phys.Rev.D56:3258-3295,1997. [arXiv:hep-ph/9704452](https://arxiv.org/abs/hep-ph/9704452), [doi:10.1103/PhysRevD.56.3258](https://doi.org/10.1103/PhysRevD.56.3258) | Teoría del preheating: resonancia paramétrica ancha en un universo en expansión, el parámetro $q$ y la resonancia estocástica para $m^2\phi^2$. | base de los casos cuadráticos |
| `Greene1997` | Greene *et al.* (1997), «Structure of Resonance in Preheating after Inflation», Phys.Rev.D56:6175-6192,1997. [arXiv:hep-ph/9705347](https://arxiv.org/abs/hep-ph/9705347), [doi:10.1103/PhysRevD.56.6175](https://doi.org/10.1103/PhysRevD.56.6175) | Estructura de bandas de resonancia para $\lambda\phi^4 + g^2\phi^2\chi^2$, donde la expansión se elimina con variables conformes. | base de los casos cuárticos ($q = g^2/\lambda$) |
| `Amin2015` | Amin *et al.* (2015), «Nonperturbative Dynamics Of Reheating After Inflation: A Review», Int. J. Mod. Phys. D 24 (2015): 1530003. [arXiv:1410.3808](https://arxiv.org/abs/1410.3808), [doi:10.1142/S0218271815300037](https://doi.org/10.1142/S0218271815300037) | Revisión de la dinámica no perturbativa del reheating: resonancia, fragmentación, oscilones, termalización. | referencia general |
| `Figueroa2017a` | Figueroa y Torrenti (2017), «Parametric Resonance in the Early Universe - A Fitting Analysis», JCAP02(2017)001. [arXiv:1609.05197](https://arxiv.org/abs/1609.05197), [doi:10.1088/1475-7516/2017/02/001](https://doi.org/10.1088/1475-7516/2017/02/001) | Ajustes de la resonancia paramétrica en simulaciones en la red para potenciales de tipo potencia: tiempos y fracciones de energía en función de $q$. | elegir $q$ y $\tilde t_{\max}$ |
| `Mishra2024` | Mishra (2024), «Cosmic Inflation: Background dynamics, Quantum fluctuations and Reheating». [arXiv:2403.10606](https://arxiv.org/abs/2403.10606) | Notas de curso: dinámica de fondo de la inflación, fluctuaciones y reheating. | referencia general |
| `Barman2024` | Barman *et al.* (2024), «Resonant Reheating». [arXiv:2404.16090](https://arxiv.org/abs/2404.16090) | Reheating resonante: producción de partículas no perturbativa y su efecto sobre $T_{\text{RH}}$. | contexto del reheating |
| `Elia2025` | Elía *et al.* (2025), «Hydrodynamic models of reheating». [arXiv:2510.08685](https://arxiv.org/abs/2510.08685) | Modelos hidrodinámicos del reheating (grupo de la UBA). | contexto del reheating y de $w$ |

## Ecuación de estado después del preheating

| clave | trabajo | qué aporta | dónde se usa |
|---|---|---|---|
| `Podolsky2006` | Podolsky *et al.* (2006), «Equation of state and Beginning of Thermalization After Preheating», Phys.Rev. D73 (2006) 023501. [arXiv:hep-ph/0507096](https://arxiv.org/abs/hep-ph/0507096), [doi:10.1103/PhysRevD.73.023501](https://doi.org/10.1103/PhysRevD.73.023501) | Simulaciones de $m^2\phi^2/2 + g^2\phi^2\chi^2/2$: $w$ salta de 0 a $\sim 0.2$--$0.3$ al final del preheating y no llega a $1/3$; un inflatón masivo termina dominando como materia. | Notas.tex, «Por qué $\langle w\rangle$ no es 0» |
| `Lozanov2017` | Lozanov y Amin (2017), «The Equation of State and Duration to Radiation Domination After Inflation», Phys. Rev. Lett. 119, 061301 (2017). [arXiv:1608.01213](https://arxiv.org/abs/1608.01213), [doi:10.1103/PhysRevLett.119.061301](https://doi.org/10.1103/PhysRevLett.119.061301) | Ecuación de estado después de la fragmentación del inflatón para $V \propto |\phi|^{2n}$: $w \to 0$ para $n = 1$ y $w \to 1/3$ para $n$ mayores. | Notas.tex, «Qué pasa después de la simulación» |

## Ondas gravitacionales del preheating

| clave | trabajo | qué aporta | dónde se usa |
|---|---|---|---|
| `Dufaux2007` | Dufaux *et al.* (2007), «Theory and Numerics of Gravitational Waves from Preheating after Inflation», Phys.Rev.D76:123517,2007. [arXiv:0707.0875](https://arxiv.org/abs/0707.0875), [doi:10.1103/PhysRevD.76.123517](https://doi.org/10.1103/PhysRevD.76.123517) | Formalismo y numérica de las GWs del preheating; reescaleo a hoy (Sec. III.D) y caso $\lambda\phi^4$ con $q = 120$. | Notas.tex, «Reescaleo a hoy» y «Comparación con la literatura»; $q = 120$ de los .in cuárticos |
| `GarciaBellido2007` | Garcia-Bellido y Figueroa (2007), «A stochastic background of gravitational waves from hybrid preheating», Phys.Rev.Lett.98:061302,2007. [arXiv:astro-ph/0701014](https://arxiv.org/abs/astro-ph/0701014), [doi:10.1103/PhysRevLett.98.061302](https://doi.org/10.1103/PhysRevLett.98.061302) | Fondo de GWs del preheating híbrido (taquiónico). | contexto |
| `Figueroa2017b` | Figueroa y Torrenti (2017), «Gravitational wave production from preheating: parameter dependence», JCAP10(2017)057. [arXiv:1707.04533](https://arxiv.org/abs/1707.04533), [doi:10.1088/1475-7516/2017/10/057](https://doi.org/10.1088/1475-7516/2017/10/057) | GWs del preheating: dependencia en $q$ y en el potencial, y reescaleo a hoy con $\epsilon$. | Notas.tex, «Comparación con la literatura» |
| `Gross2024` | Gross *et al.* (2024), «Gravitational Wave Production During Reheating: From the Inflaton to Primordial Black Holes». [arXiv:2411.04189](https://arxiv.org/abs/2411.04189) | GWs producidas durante el reheating, del inflatón a agujeros negros primordiales. | contexto |
| `Das2025` | Das y Revankar (2025), «Preheating and gravitational waves in large-field hilltop inflation», The European Physical Journal Special Topics, 2025. [arXiv:2508.07442](https://arxiv.org/abs/2508.07442), [doi:10.1140/epjs/s11734-025-01863-x](https://doi.org/10.1140/epjs/s11734-025-01863-x) | Preheating y GWs en inflación hilltop de campo grande, con simulaciones en la red. | comparación con otro tipo de modelo |

## CosmoLattice

| clave | trabajo | qué aporta | dónde se usa |
|---|---|---|---|
| `Figueroa2021` | Figueroa *et al.* (2021), «The art of simulating the early Universe -- Part I», JCAP 04 (2021) 035. [arXiv:2006.15122](https://arxiv.org/abs/2006.15122), [doi:10.1088/1475-7516/2021/04/035](https://doi.org/10.1088/1475-7516/2021/04/035) | «The art of simulating the early Universe», Parte I: discretización, integradores y variables de programa en que se basa CosmoLattice. | variables de programa, reescaleo temporal $\alpha$ |
| `Figueroa2023` | Figueroa *et al.* (2023), «CosmoLattice», Comput.Phys.Commun. 283 (2023) 108586. [arXiv:2102.01031](https://arxiv.org/abs/2102.01031), [doi:10.1016/j.cpc.2022.108586](https://doi.org/10.1016/j.cpc.2022.108586) | Artículo del código CosmoLattice (manual de uso). | cita del código |
| `Figueroa2024` | Figueroa *et al.* (2024), «Present and future of CosmoLattice», Rep. Prog. Phys. 87 (2024) 094901. [arXiv:2312.15056](https://arxiv.org/abs/2312.15056), [doi:10.1088/1361-6633/ad616a](https://doi.org/10.1088/1361-6633/ad616a) | Estado y planes de CosmoLattice. | referencia general |
| `BaezaBallesteros2025` | Baeza-Ballesteros *et al.* (2025), «The art of simulating the early Universe. Part II. Non-canonical cases & gravitational waves». [arXiv:2512.15627](https://arxiv.org/abs/2512.15627) | «The art of simulating the early Universe», Parte II: casos no canónicos y GWs. | módulo de GWs |
| `CLTechNoteGW` | J. Baeza-Ballesteros, D. G. Figueroa y N. Loayza, «CosmoLattice Technical Note II: Gravitational Waves» (2022, corregida en 2023), [cosmolattice.net](https://cosmolattice.net) | Nota técnica II de CosmoLattice: cómo se calculan y normalizan los espectros de GWs (no está en arXiv; ver cosmolattice.net). | Notas.tex, «Qué escribe CosmoLattice» |
| `BaezaBallesteros2026` | Baeza-Ballesteros *et al.* (2026), «CosmoLattice 2.0». [arXiv:2607.24978](https://arxiv.org/abs/2607.24978) | CosmoLattice v2.0: acoplamientos no mínimos, axiones y nuevas condiciones iniciales. | versión nueva del código (este repositorio usa la incluida en include/) |

## Modelos de inflación y observaciones

| clave | trabajo | qué aporta | dónde se usa |
|---|---|---|---|
| `Ellis2026` | Ellis *et al.* (2026), «Constraints on Attractor Models of Inflation and Reheating from Planck, BICEP/Keck, ACT DR6, and SPT-3G Data», Phys. Rev. D **113**, 063571. [arXiv:2510.18656](https://arxiv.org/abs/2510.18656), [doi:10.1103/d35r-7bn8](https://doi.org/10.1103/d35r-7bn8) | $\alpha$-attractors E y T generalizados frente a Planck, BICEP/Keck, ACT DR6 y SPT-3G; $N_*(T_{\text{RH}})$ (Ec. 27, App. B). | parámetros de todos los .in cuárticos y cuadráticos |
| `Saini2024` | Saini y Nautiyal (2024), «Observational constraints on $α$-Starobinsky inflation». [arXiv:2409.05615](https://arxiv.org/abs/2409.05615) | Cotas observacionales sobre $\alpha$-Starobinsky. | comparación de cotas en $\alpha_K$ |
| `Ketov2025` | Ketov (2026), «On Legacy of Starobinsky Inflation», In: Open Issues in Gravitation and Cosmology, Springer, 2026, ISBN:978-3-032-15702-7. [arXiv:2501.06451](https://arxiv.org/abs/2501.06451) | Revisión del modelo de Starobinsky. | contexto del E-model con $\alpha_K = 1$ |
| `Baumann2009` | Baumann (2009), «TASI Lectures on Inflation». [arXiv:0907.5424](https://arxiv.org/abs/0907.5424) | Notas TASI de inflación: slow-roll, perturbaciones, $n_s$ y $r$. | referencia general |

## Fondos estocásticos de GWs y reescaleo a hoy

| clave | trabajo | qué aporta | dónde se usa |
|---|---|---|---|
| `Caprini2018` | Caprini y Figueroa (2018), «Cosmological Backgrounds of Gravitational Waves». [arXiv:1801.04268](https://arxiv.org/abs/1801.04268), [doi:10.1088/1361-6382/aac608](https://doi.org/10.1088/1361-6382/aac608) | Revisión de los fondos cosmológicos de GWs, incluidos los del preheating. | referencia general |
| `Christensen2019` | Christensen (2019), «Stochastic Gravitational Wave Backgrounds», Reports on Progress in Physics, Vol. 82, 016903 (2019). [arXiv:1811.08797](https://arxiv.org/abs/1811.08797), [doi:10.1088/1361-6633/aae6b5](https://doi.org/10.1088/1361-6633/aae6b5) | Revisión de los fondos estocásticos de GWs y su detección. | referencia general |
| `Renzini2022` | Renzini *et al.* (2022), «Stochastic Gravitational-Wave Backgrounds: Current Detection Efforts and Future Prospects». [arXiv:2202.00178](https://arxiv.org/abs/2202.00178) | Estado actual y perspectivas de detección de fondos estocásticos. | contexto observacional |
| `Watanabe2006` | Watanabe y Komatsu (2006), «Improved Calculation of the Primordial Gravitational Wave Spectrum in the Standard Model», Phys.Rev. D73 (2006) 123515. [arXiv:astro-ph/0604176](https://arxiv.org/abs/astro-ph/0604176), [doi:10.1103/PhysRevD.73.123515](https://doi.org/10.1103/PhysRevD.73.123515) | Efecto de $g_*(T)$ y $g_{s*}(T)$ sobre la transferencia de las GWs hasta hoy. | constantes $C_f$ y $C_\Omega$ del reescaleo |
| `Ringwald2022` | Ringwald y Tamarit (2022), «Revealing the Cosmic History with Gravitational Waves». [arXiv:2203.00621](https://arxiv.org/abs/2203.00621), [doi:10.1103/PhysRevD.106.063027](https://doi.org/10.1103/PhysRevD.106.063027) | Cómo los fondos de GWs permiten reconstruir la historia térmica (incluido el reheating). | contexto de la dependencia en $T_{\text{RH}}$ |
| `Moore2015` | Moore *et al.* (2015), «Gravitational-wave sensitivity curves», Classical & Quantum Gravity, 32(1):015014 (2015). [arXiv:1408.0740](https://arxiv.org/abs/1408.0740), [doi:10.1088/0264-9381/32/1/015014](https://doi.org/10.1088/0264-9381/32/1/015014) | Curvas de sensibilidad de detectores de GWs en $h^2\Omega_{GW}$. | comparar las señales con detectores |

## Detección de GWs de alta frecuencia (MHz–GHz)

| clave | trabajo | qué aporta | dónde se usa |
|---|---|---|---|
| `Aggarwal2021` | Aggarwal *et al.* (2021), «Challenges and Opportunities of Gravitational Wave Searches at MHz to GHz Frequencies», Living Reviews in Relativity volume 24, Article number: 4 (2021). [arXiv:2011.12414](https://arxiv.org/abs/2011.12414), [doi:10.1007/s41114-021-00032-5](https://doi.org/10.1007/s41114-021-00032-5) | Revisión de búsquedas de GWs entre MHz y GHz: fuentes (incluido el preheating) y detectores. | dónde caen nuestras señales ($10^7$--$10^9$ Hz) |
| `Domcke2023` | Domcke (2023), «Electromagnetic high-frequency gravitational wave detection». [arXiv:2306.04496](https://arxiv.org/abs/2306.04496) | Detección electromagnética de GWs de alta frecuencia. | contexto observacional |
| `Berlin2022` | Berlin *et al.* (2021), «Detecting High-Frequency Gravitational Waves with Microwave Cavities». [arXiv:2112.11465](https://arxiv.org/abs/2112.11465), [doi:10.1103/PhysRevD.105.116011](https://doi.org/10.1103/PhysRevD.105.116011) | Detección de GWs de alta frecuencia con cavidades de microondas. | contexto observacional |
| `Herman2023` | Herman *et al.* (2023), «Electromagnetic Antennas for the Resonant Detection of the Stochastic Gravitational Wave Background», Phys. Rev. D 108, 124009 (2023). [arXiv:2203.15668](https://arxiv.org/abs/2203.15668), [doi:10.1103/PhysRevD.108.124009](https://doi.org/10.1103/PhysRevD.108.124009) | Antenas electromagnéticas para la detección resonante del fondo estocástico. | contexto observacional |
| `Gatti2024` | Gatti *et al.* (2024), «Cavity Detection of Gravitational Waves: Where Do We Stand?», Phys. Rev. D 110, 023018 (2024). [arXiv:2403.18610](https://arxiv.org/abs/2403.18610), [doi:10.1103/PhysRevD.110.023018](https://doi.org/10.1103/PhysRevD.110.023018) | Estado de la detección de GWs con cavidades. | contexto observacional |
| `Mohanty2026` | Mohanty *et al.* (2026), «Testing the Starobinsky model of inflation with resonant cavities», Phys. Rev. D 113, 063574 (2026). [arXiv:2503.06858](https://arxiv.org/abs/2503.06858), [doi:10.1103/ckhr-ffkg](https://doi.org/10.1103/ckhr-ffkg) | Cómo testear el modelo de Starobinsky con cavidades resonantes. | directamente relevante para el E-model cuadrático con $\alpha_K = 1$ |
| `Amaral2026` | Amaral *et al.* (2026), «Global detector network to search for high-frequency gravitational waves (GravNet): conceptual design». [arXiv:2603.24645](https://arxiv.org/abs/2603.24645) | Diseño conceptual de GravNet, red global de detectores de alta frecuencia. | contexto observacional |
| `Abac2026` | Abac *et al.* (2026), «The Science of the Einstein Telescope», JCAP 03 (2026) 081. [arXiv:2503.12263](https://arxiv.org/abs/2503.12263), [doi:10.1088/1475-7516/2026/03/081](https://doi.org/10.1088/1475-7516/2026/03/081) | Ciencia del Einstein Telescope. | contexto observacional |

## Libros

| clave | libro | para qué |
|---|---|---|
| `Maggiore2018` | M. Maggiore, *Gravitational Waves, Vol. 2: Astrophysics and Cosmology*, Oxford University Press (2018). | GWs cosmológicas, fondos estocásticos, preheating. |
| `Dodelson2020` | S. Dodelson y F. Schmidt, *Modern Cosmology*, 2.ª ed., Academic Press (2020). | cosmología de fondo, entropía y $g_*$. |
| `Peter2009` | P. Peter y J.-P. Uzan, *Primordial Cosmology*, Oxford University Press (2009). | inflación y reheating. |
| `Parker2009` | L. Parker y D. Toms, *Quantum Field Theory in Curved Spacetime*, Cambridge University Press (2009). | producción de partículas en un universo en expansión. |

Además, en la carpeta local de bibliografía del usuario hay trabajos que no se incluyeron porque no tratan temas del proyecto: defectos cósmicos (arXiv:1212.5458), transiciones de fase en LISA (arXiv:2403.03723), materia oscura freeze-in (arXiv:2506.08106), señales de nueva física en NANOGrav (arXiv:2306.16219), un libro escaneado sin identificar y diapositivas de un journal club.
