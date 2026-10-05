#ifndef GENEMODEL_H // Usual macro guard to prevent multiple inclusion
#define GENEMODEL_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Generalized alpha-attractor E-model, V(phi) = A|1 - exp(-phi/M)|^n, for a general
// real power n. See Notas/Notas.tex (la sección «E-models» (rama e-model)) for the derivation of the program variables and
// initial conditions. Unlike the T-model, this potential is NOT symmetric under phi -> -phi: it
// has a plateau at phi -> +infinity but grows steeply (exponentially, then like |phi|^n once
// n-th-powered) for phi -> -infinity. This is the expected, physical shape of E-models; it is not
// an implementation issue, but it does make the post-inflationary oscillations asymmetric.

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

#define MODELNAME genemodel
  // This should match the name of this file.

  template <class R> using Model = MakeModel(R, ModelPars);

  class MODELNAME : public Model<MODELNAME>
  {
  private:
    FloatType n, A, M, q, g, Mtilde;
    // n: power. A: plateau-normalization amplitude (V = A|1-exp(-phi/M)|^n). M: attractor scale.
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

      // Same alpha as the generalized monomial and T-model: near phi=0 (the region that matters
      // right after inflation, before the field feels the asymmetry of the potential) this model
      // also reduces to V ~ (A/M^n)|phi|^n. See Notas.tex «E-models», variables computacionales.
      alpha = 3 * (n - 2) / (n + 2);
      // For n > 4 (alpha > 1) the gradient term (~a^(2 alpha - 2)) and the daughter-field mass
      // (~q a^(6(n-4)/(n+2))) grow in program time, so a fixed dt eventually goes unstable; see
      // Notas.tex «Estabilidad numérica para n > 4» and code/verificar_genmodels.py::cotas_estabilidad.

      fStar = fldS0[0];
      // Same local-monomial-based omegaStar as the T-model (identical formula: both models share
      // the same near-origin behaviour). See Notas.tex «E-models», variables computacionales.
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

    auto potentialTerms(Tag<0>) const // |Mtilde*(1 - exp(-phi_tilde/Mtilde))|^n / n
    {
      return pow(abs(Mtilde * (FloatType(1) - exp(-fldS(0_c) / Mtilde))), n) / n;
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
      return pow(Mtilde, n - 1) * exp(-fldS(0_c) / Mtilde) *
                 (FloatType(1) - exp(-fldS(0_c) / Mtilde)) *
                 pow(abs(FloatType(1) - exp(-fldS(0_c) / Mtilde)), n - 2) +
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
      return pow(Mtilde, n - 2) * exp(-fldS(0_c) / Mtilde) *
                 pow(abs(FloatType(1) - exp(-fldS(0_c) / Mtilde)), n - 2) *
                 (n * exp(-fldS(0_c) / Mtilde) - FloatType(1)) +
             q * pow<2>(fldS(1_c));
    }

    auto potDeriv2(Tag<1>) // d^2 V / d chi_tilde^2
    {
      return q * pow<2>(fldS(0_c));
    }
  };
} // namespace TempLat

#endif // GENEMODEL_H
