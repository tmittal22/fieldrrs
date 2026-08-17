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

### 2.1 The seed chain — how η, S_dg and a*_φ are actually computed

Three inputs are neither fitted nor measured: they are **computed from the spectrum
itself by published parameterisations, before the inversion runs**. Figure
`giop2_seed_chain.png` shows all three, evaluated on this spectrum. Everything below
depends on one number, the subsurface blue/green ratio

> x = r_rs(443) / r_rs(555),  **= 0.30077 on the LOC1 mean**

**η — QAA v5** (`model.eta_qaa`, the `eta='qaa'` default; note it uses *subsurface*
r_rs):

> **η = 2.0 · [ 1 − 1.2 · exp(−0.9 x) ]**   → **η = 0.1692** here

The alternative is `eta='gsm'`, a constant 1.03373 — six times larger. Across the 12
per-scan LOC1 fits the QAA formula returns η = 0.127–0.215, so it is stable within the
station but sits far below the GSM constant, which is what a low blue/green ratio does.

**S_dg — four options, and they disagree by 71 %:**

| option | formula | on LOC1 |
|---|---|---|
| **0.018** | GIOP-DC default, a fixed constant | 0.0180 |
| `'qaa'` | 0.015 + 0.002 / (0.6 + x) | 0.01722 |
| `'obpg'` | clip[ 0.015 + 0.0038·log₁₀(R_rs(412)/R_rs(555)), 0.01, 0.02 ] | **0.01205** |
| `'gsm'` | constant | 0.02061 |
| — | **fitted here** (§4) | **0.0113** |

> ⚠ `'qaa'` uses *subsurface* r_rs while `'obpg'` uses *above-surface* R_rs. That
> asymmetry is real, is preserved from the upstream MATLAB deliberately, and is recorded
> in `giop_python/PORTING_NOTES.md` D3.
>
> **This is corroboration, not just disagreement.** GIOP's own `'obpg'` parameterisation
> gives 0.01205, within **7 %** of the value χ² independently prefers (0.0113), while the
> GIOP-DC default 0.018 is **59 % high**. Two routes that share no machinery — a published
> band-ratio formula and our χ² minimisation — agree that this water has a shallower CDOM
> slope than the default assumes.

**a*_φ — Bricaud et al. (1998)**, and this is the *only* thing the OC4 number does:

> **a*_φ(λ) = A_φ(λ) · chl^( E_φ(λ) − 1 )**,  rescaled so a*_φ(442) = 0.055 m² mg⁻¹

A_φ and E_φ are tabulated per wavelength; the chlorophyll enters as a **power-law
exponent**, so the seed changes the *shape*, not the amplitude. It is **frozen at the
seed and never iterated** — GIOP does not re-derive chl from its own retrieval. The shape
change is not cosmetic:

| seed chl | a*_φ(443) | a*_φ(443)/a*_φ(675) |
|---|---|---|
| 0.5 | 0.05472 | 2.349 |
| 2 | 0.05435 | 1.781 |
| **9.84 (OC4 here)** | **0.05394** | **1.296** |
| 20 | 0.05376 | 1.125 |
| 50 | 0.05352 | 0.937 |

The blue-to-red contrast of the phytoplankton shape moves by a factor 2.5 across the
plausible seed range while its 443 nm value barely moves — packaging: higher chlorophyll
means larger, more pigment-packed cells, flattening the blue peak. **That is why sweeping
the OC4 input in `giop8_assumption_free.png` changes the answer at all**, and why M_φ
ranges 0 → 41 under that sweep.

**Where OC4 itself comes from** (`empirical.get_oc`, NASA OC v6, ported from
`get_oc.m`): a 4th-order polynomial in one band ratio,

> X = log₁₀[ max(R_rs 443, 490, 510) / R_rs(555) ],  log₁₀ chl = Σ aᵢ Xⁱ

On LOC1 the numerator band selected is **510 nm**, X = −0.1878, **chl = 9.84 mg m⁻³**.
A blue/green ratio in sediment-dominated water reads suspended mineral as chlorophyll,
which is the Case-1/Case-2 problem in §5 — so the seed that sets the a*_φ shape is itself
the least trustworthy number in the chain.

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

### 4.1 "Free" is not assumption-free — the OC4 seed survives it

`fit_shapes=True` frees S_dg and η. It does **not** free a*_φ: `_invert_shapes` builds
the phytoplankton eigenvector **once** from the chlorophyll seed and holds it fixed
inside the profile. So the OC4 number — a Case-1 band ratio, on Case-2 water — still sets
the phytoplankton shape in every "free" fit reported here.

The genuinely maximal arm (`_maxfree` in `make_giop_figures.py`) profiles the seed as
well, over 10 Bricaud seeds **and** 6 Ciotti S_f values, and reports the best. On the
LOC1 mean spectrum:

| arm | free parameters | M_φ | a_dg(443) | b_bp(443) | S_dg | η | χ²_ν |
|---|---|---|---|---|---|---|---|
| constrained (GIOP-DC) | 3 amplitudes | 11.47 | 1.254 | 0.0836 | 0.018 *fixed* | 0.169 *QAA* | **74.5** |
| free (still OC4-seeded) | + S_dg, η | 2.26 | 0.779 | 0.0430 | 0.01176 | −1.000 | **18.1** |
| **max freedom** (Bricaud seed 3) | + a*_φ family/seed | 1.68 | 0.781 | 0.0416 | 0.01144 | −1.000 | **17.2** |

Two things follow, and they point in opposite directions:

- **The CDOM slope matters far more than the seed.** Releasing S_dg and η took χ²_ν from
  74.5 to 18.1; additionally releasing the a*_φ family took it only to 17.2. So the a*_φ
  prescription — the assumption most obviously wrong on Case-2 water — is **not** what
  the misfit is made of.
- **a_dg and b_bp are stable across the two free arms** (0.779 → 0.781, and 0.0430 →
  0.0416, i.e. 0.2 % and 3 %). Once the shapes are free, those two stop caring about the
  seed. **M_φ does not**: 2.26 → 1.68, and it was 11.47 constrained.

### 4.2 "If you know chl from OC4, don't you know the shape *and* the amplitude?"

Yes — and GIOP deliberately throws the amplitude away. Bricaud gives the **absolute**
absorption

> a_φ(λ) = A_φ(λ) · chl^E_φ(λ)

GIOP divides by chl to get the **specific** a*_φ(λ) = A_φ(λ)·chl^(E_φ−1), renormalises it
to 0.055 at 442 nm, and then fits a free amplitude M_φ on top. Because of that
normalisation M_φ carries units of chlorophyll, so if the seed were believed you would
simply set **M_φ = chl** and phytoplankton would cost zero free parameters. GIOP does not,
because inheriting a Case-1 band ratio is the thing a semi-analytical inversion exists to
avoid.

That makes **M_φ vs its own seed a real internal-consistency test**, and GIOP never runs
it — the docstring is explicit that the seed is frozen and not iterated. Sweeping the
seed and asking where M_φ = seed (panel 7 of `giop10_final_result.png`) gives **two**
roots:

| root | chl | stability |
|---|---|---|
| lower | 1.22 | **unstable** — this is where M_φ has collapsed onto its 0 bound; iterate away and you never return |
| upper | **11.18** | **stable** — iterating the seed converges here |

> ⚠ Taking the *first* sign change reports the artefact. The first version of this
> analysis did exactly that and would have published "self-consistent chl = 1.2 against
> OC4 9.8", a factor-8 contradiction that does not exist.

**The stable fixed point is 11.18 against OC4's 9.84 — agreement to 14 %.** So in the
constrained configuration GIOP and OC4 *are* mutually consistent, and the earlier
"factor-2.4 internal contradiction" is dead twice over (once by going hyperspectral,
once by this).

But that consistency is a property of the **constrained** model, which fits at χ²_ν = 74.
The arms that actually fit want **M_φ ≈ 1.7–2.3, roughly 5× less phytoplankton than OC4
reports**. Two Case-1 relations agreeing with each other is not evidence about the water;
it is evidence that they share a calibration.

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
