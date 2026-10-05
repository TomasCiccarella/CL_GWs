#ifndef GENMONOMIAL_H // Usual macro guard to prevent multiple inclusion
#define GENMONOMIAL_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Generalized monomial potential V(phi) = A|phi|^n, for a general real power n.
// See Notas/Notas.tex (la sección «Modelo monomial» (rama monomial)) for the derivation of the program variables and initial conditions,
// and Figueroa, Florio, Torrenti & Valkenburg, arXiv:2006.15122, Sec. 8, for the general-p rescaling.

#include "CosmoInterface/cosmointerface.h"

namespace TempLat
{
  /////////
  // Model name and number of fields
  /////////

  struct ModelPars : public TempLat::DefaultModelPars {
    static constexpr size_t NScalars = 2;
    // Inflaton (0) + daughter field (1).
    static constexpr size_t NPotTerms = 2;
    // Inflaton potential term + interaction term.

    using FloatType = double;
  };

#define MODELNAME genmonomial
  // This should match the name of this file.

  template <class R> using Model = MakeModel(R, ModelPars);

  class MODELNAME : public Model<MODELNAME>
  {
  private:
    FloatType n, A, q, g;
    // n: power of the monomial. A: amplitude (V = A|phi|^n). q: dimensionless resonance
    // parameter (V_int = 0.5*g^2*phi^2*chi^2, with q = g^2*phi0^(4-n)/(n*A)). g is derived from q.

  public:
    MODELNAME(ParameterParser &parser, RunParameters<FloatType> &runPar, auto toolBox)
        : Model<MODELNAME>(parser, runPar.getLatParams(), toolBox, runPar.dt,
                           STRINGIFY(MODELLABEL))
    {
      /////////
      // Independent parameters of the model and initial homogeneous components of the fields
      /////////

      n = parser.get<FloatType>("n");
      A = parser.get<FloatType>("A");
      q = parser.get<FloatType>("q");

      fldS0 = parser.get<FloatType, 2>("initial_amplitudes");
      piS0 = parser.get<FloatType, 2>("initial_momenta", {0, 0});

      /////////
      // Rescaling for program variables
      /////////

      // alpha = 3(n-2)/(n+2) keeps the inflaton oscillation period approximately constant in
      // program time (self-similar monomial oscillations, see Notas.tex «Reescaleo temporal»). This
      // reproduces alpha=0 for n=2 (tanh2) and alpha=1 for n=4 (lphi4).
      alpha = 3 * (n - 2) / (n + 2);
      // For n > 4 (alpha > 1) the gradient term (~a^(2 alpha - 2)) and the daughter-field mass
      // (~q a^(6(n-4)/(n+2))) grow in program time, so a fixed dt eventually goes unstable; see
      // Notas.tex «Estabilidad numérica para n > 4» and code/verificar_genmodels.py::cotas_estabilidad.

      fStar = fldS0[0];
      // omegaStar is the physical oscillation frequency Omega(phi) = sqrt(V'(phi)/phi) evaluated
      // at phi=fStar (see Notas.tex «Variables de programa»); this is what makes alpha above actually keep the
      // oscillation period constant in program time. Reduces to sqrt(lambda)*fStar for n=4,
      // A=lambda/4, matching lphi4.h exactly.
      omegaStar = sqrt(n * A) * pow(fStar, n / 2 - 1);

      g = sqrt(q) * omegaStar / fStar;
      // g is the physical coupling, kept only for reference/diagnostics; the lattice dynamics
      // only use q directly.

      setInitialPotentialAndMassesFromPotential();
    }

    /////////
    // Program potential
    /////////

    auto potentialTerms(Tag<0>) const // Inflaton potential energy: |phi_tilde|^n / n
    {
      return pow(abs(fldS(0_c)), n) / n;
    }
    auto potentialTerms(Tag<1>) const // Interaction energy: (1/2) q phi_tilde^2 chi_tilde^2
    {
      return FloatType(0.5) * q * pow<2>(fldS(0_c) * fldS(1_c));
    }

    /////////
    // Derivatives of the program potential
    /////////

    auto potDeriv(Tag<0>) // d V / d phi_tilde
    {
      return fldS(0_c) * pow(abs(fldS(0_c)), n - 2) + q * fldS(0_c) * pow<2>(fldS(1_c));
    }

    auto potDeriv(Tag<1>) // d V / d chi_tilde
    {
      return q * fldS(1_c) * pow<2>(fldS(0_c));
    }

    /////////
    // Second derivatives of the program potential
    /////////

    auto potDeriv2(Tag<0>) // d^2 V / d phi_tilde^2
    {
      return (n - 1) * pow(abs(fldS(0_c)), n - 2) + q * pow<2>(fldS(1_c));
    }

    auto potDeriv2(Tag<1>) // d^2 V / d chi_tilde^2
    {
      return q * pow<2>(fldS(0_c));
    }
  };
} // namespace TempLat

#endif // GENMONOMIAL_H
