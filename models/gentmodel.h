#ifndef GENTMODEL_H // Usual macro guard to prevent multiple inclusion
#define GENTMODEL_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Generalized alpha-attractor T-model, V(phi) = A|tanh(phi/M)|^n, for a general real
// power n (generalizes tanh2.h, which is the fixed n=2 case). See Notas/Notas.tex (Sec. 2) for the
// derivation of the program variables and initial conditions, and Figueroa, Florio, Torrenti &
// Valkenburg, arXiv:2006.15122, Sec. 8, for the same model (there called p instead of n).

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

#define MODELNAME gentmodel
  // This should match the name of this file.

  template <class R> using Model = MakeModel(R, ModelPars);

  class MODELNAME : public Model<MODELNAME>
  {
  private:
    FloatType n, A, M, q, g, Mtilde;
    // n: power. A: plateau-normalization amplitude (V = A|tanh(phi/M)|^n). M: attractor scale.
    // q: dimensionless resonance parameter (V_int = 0.5*g^2*phi^2*chi^2, with
    // q = g^2*phi0^(4-n)*M^n/(n*A)). g is derived from q.
    // Mtilde = M/fStar is the dimensionless attractor scale in program units.

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
      M = parser.get<FloatType>("M");
      q = parser.get<FloatType>("q");

      fldS0 = parser.get<FloatType, 2>("initial_amplitudes");
      piS0 = parser.get<FloatType, 2>("initial_momenta", {0, 0});

      /////////
      // Rescaling for program variables
      /////////

      // Same alpha as the generalized monomial model: during reheating the field oscillates in
      // the region |phi| << M, where V ~ (A/M^n)|phi|^n, so the same self-similar-oscillation
      // argument applies (see Notas.tex Sec. 2.2).
      alpha = 3 * (n - 2) / (n + 2);
      // For n > 4 (alpha > 1) the gradient term (~a^(2 alpha - 2)) and the daughter-field mass
      // (~q a^(6(n-4)/(n+2))) grow in program time, so a fixed dt eventually goes unstable; see
      // Notas.tex Sec. 1.6 and code/verificar_genmodels.py::cotas_estabilidad.

      fStar = fldS0[0];
      // omegaStar is the physical oscillation frequency Omega(phi) = sqrt(V_local'(phi)/phi)
      // evaluated at phi=fStar, using the *local* monomial behaviour near phi=0 (effective
      // amplitude A/M^n) rather than the full bounded tanh^n term -- this is what makes alpha
      // above actually keep the oscillation period constant in program time; it reduces exactly
      // to tanh2.h's omegaStar for n=2 (with A=Lambda4/2). See Notas.tex Sec. 2.2.
      omegaStar = sqrt(n * A) * pow(M, -n / 2) * pow(fStar, n / 2 - 1);

      Mtilde = M / fStar;

      g = sqrt(q) * omegaStar / fStar;
      // g is the physical coupling, kept only for reference/diagnostics; the lattice dynamics
      // only use q directly.

      setInitialPotentialAndMassesFromPotential();
    }

    /////////
    // Program potential
    /////////

    auto potentialTerms(Tag<0>) const // |Mtilde*tanh(phi_tilde/Mtilde)|^n / n
    {
      return pow(abs(Mtilde * tanh(fldS(0_c) / Mtilde)), n) / n;
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
      return pow(Mtilde, n - 1) * tanh(fldS(0_c) / Mtilde) *
                 pow(abs(tanh(fldS(0_c) / Mtilde)), n - 2) / pow<2>(cosh(fldS(0_c) / Mtilde)) +
             q * fldS(0_c) * pow<2>(fldS(1_c));
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
      return pow(Mtilde, n - 2) * pow(abs(tanh(fldS(0_c) / Mtilde)), n - 2) /
                 pow<2>(cosh(fldS(0_c) / Mtilde)) *
                 ((n - 1) / pow<2>(cosh(fldS(0_c) / Mtilde)) - 2 * pow<2>(tanh(fldS(0_c) / Mtilde))) +
             q * pow<2>(fldS(1_c));
    }

    auto potDeriv2(Tag<1>) // d^2 V / d chi_tilde^2
    {
      return q * pow<2>(fldS(0_c));
    }
  };
} // namespace TempLat

#endif // GENTMODEL_H
