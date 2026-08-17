# From R_rs to concentrations: the GIOP inversion as applied here

What the inversion actually solves, which quantities are **fitted**, which are
**assumed**, and which are neither. Companion to `THEORY_SCALED_MEAN.md` (which covers
how the R_rs being inverted was formed) and `LOC1_GIOP_FINDINGS.md` (the verdict on
LOC1). The full port-level theory, every equation mapped to the upstream MATLAB line, is
in the sibling repo's `giop_python/THEORY.md`.

Upstream MATLAB: <https://github.com/kelseybisson/GIOP> @ `ef9b93f`.
Our port: <https://github.com/tmittal22/giop-workbench>. Reference: Werdell et al. (2013),
*Generalized ocean color inversion model for retrieving marine inherent optical
properties*, Applied Optics 52(10), 2019–2037, doi:10.1364/AO.52.002019.

---

## 1. The chain, in four equations

**(1) Above water to below water.** R_rs is measured in air; the radiative transfer is
below the surface. Lee et al. (2002), QAA:

> r_rs = R_rs / (0.52 + 1.7 R_rs)

**(2) The AOP–IOP relation.** Gordon et al. (1988), quadratic in the single variable u:

> r_rs = g₀ u + g₁ u²,  with g₀ = 0.0949, g₁ = 0.0794
> u = b_b / (a + b_b)

**This is the only thing R_rs measures.** Everything after it is library, not
measurement. §3 below is the consequence.

**(3) The decomposition.** Total absorption and backscatter are split into a known water
term plus three constituents, each written as a **fixed spectral shape times a free
amplitude**:

> a(λ) = a_w(λ) + M_φ·a*_φ(λ) + M_dg·exp[−S_dg(λ − 443)]
> b_b(λ) = b_bw(λ) + M_bp·(443/λ)^η

**(4) Inversion.** Find (M_dg, M_bp, M_φ) minimising the misfit between the modelled and
the measured r_rs. Three free amplitudes; here, 301 wavelengths.

## 2. What is fitted, what is assumed, and the naming trap

| symbol | role | in this analysis |
|---|---|---|
| **M_dg** | fitted amplitude | free |
| **M_bp** | fitted amplitude | free |
| **M_φ** | fitted amplitude | free |
| S_dg | shape parameter of a_dg | GIOP-DC **assumes** 0.018 nm⁻¹. We also fit it — see §4 |
| η | shape parameter of b_bp | derived per spectrum by QAA, or fitted — §4 |
| a*_φ(λ) | shape of phytoplankton absorption | Bricaud et al. (1998), selected by an OC4 chlorophyll seed; Ciotti & Bricaud (2006) as the size-parameterised alternative |
| a_w, b_bw | pure water | tabulated, never fitted |

> **`adg443` ≡ M_dg and `bbp443` ≡ M_bp, exactly.** Both shapes are normalised to 1 at
> 443 nm by construction — exp[−S_dg·0] = 1 and (443/443)^η = 1 — so the fitted amplitude
> *is* the coefficient at 443 nm. The same is **not** true of phytoplankton: Bricaud's
> a*_φ is chlorophyll-*specific* absorption, so a*_φ(443) = 0.05394, and
> a_φ(443) = M_φ × 0.05394 ≠ M_φ. M_φ carries units of chlorophyll; the other two carry
> units of m⁻¹.

## 3. Why the answer is weaker than it looks

Equation (2) shows R_rs constrains **u(λ) = b_b/(a+b_b)** and nothing else. Splitting u
into a and b_b, and then a into three additive smooth curves, is done entirely by the
prescribed shapes. Two consequences that show up in the LOC1 figures:

- **a_dg and a_φ trade off directly.** Both rise toward the blue, so an increase in one
  can be absorbed by a decrease in the other at nearly constant χ². Across the 12
  per-scan fits their correlation is **r = −0.92** (`giop7_covariance.png`).
- **A band ratio and a magnitude carry different errors.** OC4 is a ratio, so it is
  *exactly invariant* to an overall scaling of R_rs; the amplitudes are not. This is why
  the 1.7 % shape uncertainty and the 11 % amplitude spread from `THEORY_SCALED_MEAN.md`
  must be propagated separately — nearly all of the amplitude term lands on b_bp.

**Report u(λ) alongside any retrieval.** It carries only the AOP–IOP operator and the
air–water transfer, both independently validated, and no assumption about what is in the
water.

## 4. Freeing the shapes, and the χ² treatment

`giop(..., fit_shapes=True)` also fits S_dg and η. It is solved as a **profile**: at each
trial (S_dg, η) the three amplitudes are solved exactly by the bounded trust-region
solver, and only the two shape parameters are searched. The configured shapes are always
among the starts, so the returned cost is **≤ the fixed-shape cost by construction** —
freeing parameters cannot make the optimum worse, since the fixed-shape solution is a
point inside the free search box.

> That property is a test, not a remark. A previous joint 5-D simplex version violated it
> by 33× and the failure was written up as a result about the data. See
> `LOC1_GIOP_FINDINGS.md` §2f and `giop_python/tests/test_fit_shapes_guard.py`.

**How to grade the arms.** χ² against the *measured* per-band uncertainty (column 3 of
`FINAL_Rrs.csv`), not against an assumed noise floor:

> χ² = Σ_λ [(R_rs^model(λ) − R_rs^obs(λ)) / σ(λ)]²,  ν = n − 3

Three cautions, all of which bind on real data:

1. **Grade each solver on the objective it actually minimised.** `inv='bounded'` with a
   σ minimises a weighted sum of squares; `inv='fmin'` minimises an unweighted one.
   Comparing across them on the wrong objective produces a spurious "nesting violation".
2. **χ² ranks; it does not assign probabilities.** The residual here has lag-1
   autocorrelation **0.9964** across 301 bands — one smooth systematic curve, not noise.
   Formal Δχ² intervals assume independent residuals and are meaningless in that regime.
   Under an AR(1) reading, n_eff ≈ 0.5.
3. **Inflate errors before weighting, and say that you did.** With χ²_ν,min ≫ 1 the
   misfit is model inadequacy, so raw exp(−Δχ²/2) puts all weight on one arm. Inflating
   by ŝ = χ²_ν,min (Avni 1976) is the standard treatment — and on LOC1 it *still*
   collapses to one arm, which is itself the finding: the assumption sweep is a set of
   rejected models, not an uncertainty band.

## 5. Where the model stops

**Bricaud's a*_φ table is defined 400–700 nm.** GIOP therefore cannot be run redward of
700 nm at all, and the inversion refuses rather than extrapolating. The instrument
records to 2500 nm, so on turbid water the two clearest features — the ~700 nm
particulate peak and the ~810 nm water-absorption shoulder — are outside the inversion
entirely. On LOC1 that is where the residual is worst (a −20σ notch at 690 nm, identical
in all 12 scans).

**GIOP is a Case-1-calibrated model.** Both OC4 and the Bricaud a*_φ come from waters
where phytoplankton co-varies with everything else. Sediment-dominated water violates
that, and a blue/green ratio there reads suspended mineral as chlorophyll.

For water like this the honest next step is a turbid-water inversion that extends past
700 nm (QAA-turbid, or `titanspec`'s registry with a mineral component), not a better-tuned
GIOP.

## References

- Werdell, P.J., et al. (2013), Applied Optics 52(10), 2019–2037, doi:10.1364/AO.52.002019
- Gordon, H.R., et al. (1988), *A semianalytic radiance model of ocean color*, JGR 93(D9),
  10909–10924, doi:10.1029/JD093iD09p10909
- Lee, Z., Carder, K.L., Arnone, R.A. (2002), *Deriving inherent optical properties from
  water color: a multiband quasi-analytical algorithm for optically deep waters*, Applied
  Optics 41(27), 5755–5772, doi:10.1364/AO.41.005755
- Bricaud, A., Morel, A., Babin, M., Allali, K., Claustre, H. (1998), *Variations of light
  absorption by suspended particles with chlorophyll a concentration*, JGR 103(C13),
  31033–31044, doi:10.1029/98JC02712
- Ciotti, A.M., Bricaud, A. (2006), *Retrievals of a size parameter for phytoplankton and
  spectral light absorption by colored detrital matter*, Limnology and Oceanography:
  Methods 4(7), 237–253, doi:10.4319/lom.2006.4.237
- Avni, Y. (1976), *Energy spectra of X-ray clusters of galaxies*, ApJ 210, 642–646,
  doi:10.1086/154870
