# The amplitude-normalised mean

Full derivation of `fieldrrs.rrs.scaled_mean`: what it estimates, why a plain average is
the wrong summary of replicate R_rs spectra, how the estimator follows from the model,
what it is equivalent to, and when it lies.

---

## 1. The observation that motivates it

At LOC1 (Kotzebue, 2026-08-16, 12 water scans over 17 minutes) the scan-to-scan variation
in R_rs is almost entirely **multiplicative**. Three measurements say so:

| quantity | value |
|---|---|
| variance carried by the leading SVD mode | **98.3 %** |
| relative scatter, raw | 11.2 % |
| relative scatter after normalising each spectrum at 555 nm | **3.0 %** |
| correlation of spectral distance with time gap | r = +0.03 (p = 0.81) |

So the spectra differ by a scale factor, not by shape, and not as a function of time. The
footprint is 0.16 m across on a turbid nearshore surface: each scan samples a different
patch with a different amount of the same suspended material.

That is a statement about **concentration**, not about **water type**, and the two should
not be reported with one error bar.

---

## 2. Model

Index scans by *i* = 1…*n* and wavelength by *j* = 1…*m*. Write the measured spectra as

$$R_i(\lambda_j) \;=\; c_i\,S(\lambda_j) \;+\; \varepsilon_{ij}
\tag{1}$$

- **S(λ)** the common **shape**, the quantity that identifies the water. Defined up to a
  scale, which §4 fixes.
- **c_i** the per-scan **amplitude**, proportional to the amount of scattering material in
  that footprint. This is a real environmental variable, not an error.
- **ε_ij** measurement noise, assumed zero-mean and independent of *i*.

Equation (1) says the data matrix **X** (n × m) is **rank one plus noise**. Everything
below is the consequence.

---

## 3. Why the plain mean misreports the uncertainty

The plain average is

$$\bar R(\lambda_j) = \frac{1}{n}\sum_i R_i(\lambda_j) = \bar c\,S(\lambda_j) + \bar\varepsilon_j
\tag{2}$$

which is an unbiased estimate of the shape up to the constant $\bar c$. **The mean itself
is fine.** The problem is the uncertainty attached to it. The band-wise standard deviation
of the replicates is

$$\mathrm{sd}_j \;=\; S(\lambda_j)\,\mathrm{sd}(c) \;+\; O(\varepsilon)
\tag{3}$$

so the *relative* scatter is

$$\frac{\mathrm{sd}_j}{\bar R_j} \;=\; \frac{\mathrm{sd}(c)}{\bar c}
\qquad\text{at every wavelength.}
\tag{4}$$

The amplitude scatter is reproduced identically in every band and reads as though the
spectrum were uncertain by that much *everywhere*. It is not. Under (1) the **shape** is
known far better than (4) suggests, and (4) is really a measurement of how much the
concentration varied.

This is why the LOC1 error bar looked like 11.7 % when the shape is determined to 1.7 %.

---

## 4. The estimator

Fit (1) by least squares. Rescale each spectrum onto a common shape and ask for the
scaling and the shape that agree best:

$$J(S, \{a_i\}) \;=\; \sum_{i=1}^{n}\sum_{j=1}^{m} w_j\,\bigl(a_i R_{ij} - S_j\bigr)^2
\tag{5}$$

with optional band weights $w_j \ge 0$ (default 1). Here $a_i \approx 1/c_i$ is the factor
that takes scan *i* onto the shape.

(5) is bilinear, so it has no closed form in both arguments at once but a closed form in
each. Alternate.

**Step A — amplitudes, holding S fixed.** Differentiate (5) with respect to $a_i$:

$$\frac{\partial J}{\partial a_i} = 2\sum_j w_j R_{ij}\,(a_i R_{ij} - S_j) = 0$$

$$\boxed{\;a_i \;=\; \frac{\langle R_i, S\rangle_w}{\langle R_i, R_i\rangle_w}\;}
\qquad \langle u,v\rangle_w \equiv \sum_j w_j u_j v_j
\tag{6}$$

**Step B — shape, holding the amplitudes fixed.** Differentiate with respect to $S_j$:

$$\frac{\partial J}{\partial S_j} = -2\sum_i w_j (a_i R_{ij} - S_j) = 0$$

$$\boxed{\;S_j \;=\; \frac{1}{n}\sum_i a_i R_{ij}\;}
\tag{7}$$

which is just the mean of the rescaled spectra, independent of the weights.

**Step C — fix the scale.** (5) is invariant under $S \to \kappa S$, $a_i \to \kappa a_i$,
so without a constraint the iteration drifts (and $S=0$, $a=0$ is a global minimum). Impose

$$\frac{1}{n}\sum_i a_i = 1
\tag{8}$$

by dividing both $\{a_i\}$ and, implicitly, the resulting $S$ by the mean of the
amplitudes. This choice keeps **S in the original physical units of sr⁻¹** and makes it
the spectrum *at the mean amplitude of the station*, which is the quantity you want to
report. Any other normalisation (unit norm, unit value at a reference band) is a rescaling
of the same answer.

Iterate A → B → C until $\max_j |S_j^{(k+1)} - S_j^{(k)}| / \max_j|S_j^{(k)}| < \text{tol}$.
On the LOC1 data this converges in **4 iterations**; on synthetic pure-amplitude data, in 2.

### Why it converges

Each of Steps A and B is the exact minimiser of (5) in its argument, so $J$ is
non-increasing and bounded below by 0, hence convergent. Step C is a reparameterisation
along the invariant direction and does not change $J$. This is block coordinate descent on
a bilinear objective; it converges to a stationary point, and for a rank-1-dominant matrix
the stationary point is the global one, reached from the plain mean as the starting value.

---

## 5. What it is equivalent to

### 5.1 The rank-1 SVD

The standard rank-1 approximation minimises the residual on the *unscaled* data,

$$J_{\rm SVD}(S,\{c_i\}) = \sum_{ij}\bigl(R_{ij} - c_i S_j\bigr)^2
\tag{9}$$

whose solution is the leading singular triplet, $c_i S_j = \sigma_1 u_{i1} v_{j1}$. The
difference from (5) is which side carries the scaling, hence the denominator in (6):

| | amplitude estimate |
|---|---|
| this estimator, eq. (6) | $a_i = \langle R_i, S\rangle / \langle R_i, R_i\rangle$ |
| rank-1 SVD, eq. (9) | $c_i = \langle R_i, S\rangle / \langle S, S\rangle$ |

(5) weights a bright scan *less* once rescaled, since rescaling amplifies its noise too;
(9) is the maximum-likelihood fit for iid Gaussian noise on $R$. **(9) is the more
principled estimator for model (1).**

Measured on LOC1, they are indistinguishable in practice:

| | residual shape scatter, 450–700 nm |
|---|---|
| alternating, eq. (5)–(8) | **1.747 %** |
| rank-1 SVD, eq. (9) | **1.739 %** |

with the two means agreeing to **2.0 %** worst case and the amplitudes to 1.5 %, and the
leading mode carrying **99.93 %** of the total energy.

**`scaled_mean` implements (5)–(8) rather than the SVD** for one reason: `fieldrrs` has no
third-party dependencies, so that it can ship as a single Windows executable, and (5)–(8)
is a few dot products in the standard library while (9) needs LAPACK. The cost of that
choice is 0.5 % on the residual, which is measured above rather than assumed. Use the SVD
if numpy is already in your stack; you will get the same answer.

### 5.2 Multiplicative scatter correction

In near-infrared spectroscopy the same problem — replicate spectra of one substance
differing by a multiplicative path-length or particle-size factor — is treated by
**Multiplicative Scatter Correction** (Geladi, MacDougall & Martens 1985,
doi:10.1366/0003702854248656), and its extension EMSC. MSC regresses each spectrum against
a reference (usually the mean) and divides out the slope. Equations (6)–(7) are MSC with
the reference **iterated to self-consistency** rather than fixed at the first mean, and
with the additive offset term omitted, since an additive offset in R_rs is a residual
skylight or glint error that belongs in ρ and Δ, not here.

---

## 6. What the two outputs mean

$$R_{rs}(\lambda) \;=\; \underbrace{S(\lambda)}_{\text{shape} \;\pm\; \sigma_S}
\;\times\; \underbrace{\bar c \,(1 \pm \sigma_c)}_{\text{amplitude}}$$

- **σ_S — shape scatter.** Residual relative scatter after rescaling. This is the
  uncertainty on *what the water is*. **Band ratios inherit this.**
- **σ_c — amplitude scatter.** Spread of $1/a_i$. This is *how much of it there is*, and it
  is a property of the water body, not of the instrument. **Absolute magnitudes inherit
  this.**

At LOC1: **σ_S = 1.7 %** over 450–700 nm, **σ_c = 11.4 %**. A plain average reports a
single 11.7 % band and so understates a ratio product by a factor 7 while being about
right for an absolute one.

### Reporting σ_S honestly

A *relative* scatter is meaningless where the signal approaches zero, and R_rs does so
below ~430 nm and above ~800 nm. Quoting σ_S over the full 400–900 nm range gives 7.1 %
for the same data, entirely from those tails. `scaled_mean` returns the all-band figure and
`analyse_location` prints **both**, with the core-band value as the headline and the
full-range value beside it, so the choice of window is visible rather than flattering.

---

## 7. When this is the wrong tool

The estimator assumes model (1). It will happily return an answer when (1) is false, and
the answer will be a fiction: a shape that is an average of genuinely different waters,
with a small σ_S implying confidence that is not there.

**The diagnostic is σ_S itself, compared against the unscaled spread.** If rescaling does
not reduce the scatter substantially, the spectra differ in shape and the mean shape is not
a meaningful object. Stated in the docstring and reproduced here:

| situation | σ_S vs raw | what to do |
|---|---|---|
| one water, varying load | falls sharply (11.2 → 1.7 %) | use it; report both terms |
| two water types mixed | barely falls | do not average; cluster first |
| a spectrally structured error (e.g. an ρ or sky-subtraction residual) | falls partly | fix the error; it is not amplitude |

Two further limits:

- **It cannot separate a genuine concentration change from a multiplicative instrument
  error**, because both are exactly a scale factor. At LOC1 the panel replicates bound the
  instrument at 0.6 %, which is 18× below the 11.2 % observed, so the attribution to water
  is safe *there* — it is not automatic.
- **It assumes the noise is not itself multiplicative.** If it were, Step A would be
  biased; the near-machine-precision recovery in §8 shows this is not an issue at the
  signal levels here.

---

## 8. Verification

**A control that must succeed exactly.** Five spectra built as $a\cdot S$ with
$a = 0.8, 0.9, 1.0, 1.15, 1.3$ and no noise:

```
shape_cv   = 4.5e-14 %          (machine precision)
amplitudes = 0.800, 0.900, 1.000, 1.150, 1.300   recovered exactly
iterations = 2
```

**A control that must fail.** If the input spectra differ in shape rather than amplitude,
σ_S must *not* collapse — otherwise the estimator would be manufacturing agreement. This
is the property the table in §7 rests on and is what the diagnostic in the docstring tells
the user to check.

**On real data.** LOC1, 12 scans: converges in 4 iterations, σ_S = 1.7 % against a raw
11.2 %, and agrees with the independent rank-1 SVD to 2.0 % (§5.1).

---

## 9. Where it sits in the pipeline

The scaled mean is the **last** step, after everything that could put a spurious
multiplicative factor into the data has been removed:

1. `rrs_three_scan` — R_rs per scan, with ρ evaluated at that scan's own view angle
   (`rho_at_angle`), so pointing does not masquerade as amplitude.
2. `match_by_angle` — one sky per water, mirrored geometry, within one location and one
   panel-reference block.
3. **`scaled_mean`** — split what remains into shape and amplitude.

Applying it before step 1 would fold the ρ(θ) systematic into the amplitude term and
report a real instrument effect as water variability.

---

## References

- Geladi, MacDougall & Martens (1985), *Linearization and scatter-correction for
  near-infrared reflectance spectra of meat*, Applied Spectroscopy 39, 491.
  doi:10.1366/0003702854248656
- Mobley (1999), *Estimation of the remote-sensing reflectance from above-surface
  measurements*, Applied Optics 38, 7442. doi:10.1364/AO.38.007442
- Quan & Fry (1995), *Empirical equation for the index of refraction of seawater*,
  Applied Optics 34, 3477. doi:10.1364/AO.34.003477
