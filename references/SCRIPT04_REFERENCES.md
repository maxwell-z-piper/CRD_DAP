# Script 04 citation record

Retained 2026-10-24. The implementation and its departures are described in
METHODS.md section 7. The [pipeline-wide index](README.md) links the other
papers and records their roles; reusable citation metadata is in
[script04.bib](script04.bib).

## Curvature penalty

[Sellwood & Zánmar Sánchez (2010), Quantifying non-circular streaming motions
in disc galaxies](https://doi.org/10.1111/j.1365-2966.2010.16430.x), MNRAS 404,
1733–1744. [Full preprint](https://arxiv.org/html/0912.5493v2).

Appendix A defines the squared-second-difference penalty and its normalization.
Section 4.4, Table 2 and footnote 6 discuss the empirical coefficient and its
sensitivity. This does not provide a universally calibrated strength.
Our central-origin triplet omission, fixed C across extents, two-component
spectral-cube objective and human-reviewed strength sweep are explicit project
choices. Borrowing the penalty does not validate the full CRD procedure.

## Observed brightness proxy

[DiGiorgio Zanger et al. (2024), The Strength of Bisymmetric Modes in
SDSS-IV/MaNGA Barred Galaxy Kinematics](https://arxiv.org/html/2407.11908v1).

Section III, equations 4–6 and the following paragraph, uses intensity-weighted
PSF convolution and substitutes observed brightness for intrinsic brightness.
Our common observed RH3 proxy is analogous; it does not separate the disks'
intrinsic photometry. We do not adopt that paper's full barred-disk model,
Bayesian sampler or data exclusions.

## Ring-model framework

[López-Cobá, Lin & Sánchez (2024), XookSuut: A Bayesian Tool for Modeling
Circular and Non-Circular Flows on 2D Velocity Maps](https://doi.org/10.22201/ia.01851101p.2024.60.01.03),
Revista Mexicana de Astronomía y Astrofísica 60, 19–39.
[Preprint](https://arxiv.org/abs/2110.05095).

This is a reference for the non-parametric radial-ring framework.
Script 04 is not an invocation of XookSuut and does not reproduce its posterior
inference machinery. The legacy BibTeX key `LopezCoba2023XookSuut`
is retained; the published citation year is 2024.
