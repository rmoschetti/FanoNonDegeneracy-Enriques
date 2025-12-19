# Part 4: Examples from the paper

This folder contains the computational examples discussed in the paper. Unlike Parts 1–3, these scripts are not meant to be run as fully automated pipelines producing a fixed dataset. They are intended to be explored **interactively**, executing them step-by-step while following the in-code comments.

## `SingularMukaiModels.sage` (Example 5.9)

**Dependency note.** This script calls

```python
load('cnd_finder.py')
```

The file `cnd_finder.py` is **not** included in this repository. It is part of *CndFinder* [MRS22b] and can be obtained from the upstream source.

Example: an Enriques surface with finite automorphism group (type **KVI**) admitting a Mukai model which is

- **non-degenerate** (given by an isotropic sequence of 9 half-fibers), and
- **singular** (it contracts a curve).

Since the Mukai polarization arises from a length-9 isotropic sequence, its orthogonal complement in $\mathrm{Num}(Y)$ is a $(-2)$-class; in this example it is effective, hence represented by a smooth rational curve that is contracted by the Mukai polarization.

## `SmoothMukaiSingularFano.sage` (Example 5.10)

Example: the $(\tau,\overline{\tau})$-generic Enriques surface $Y_{158}$, for which every Fano model is singular (since $\mathrm{Fnd}(Y_{158}) = 9$), but there exists a Mukai model which is

- **non-degenerate** (given by an isotropic sequence of 9 half-fibers), and
- **smooth** (it does not contract any curve).

The key point is that, for a non-degenerate length-9 sequence, the only $(-2)$-class that can be contracted is the orthogonal class. In this case, the script verifies that the orthogonal $(-2)$-class is **not** effective, hence no curve is contracted by the Mukai polarization.


## Bibliography

[MRS22b] Moschetti, Riccardo; Rota, Franco; Schaffler, Luca. *SageMath code CndFinder*. Available at:
[https://github.com/rmoschetti/CNDFinder](https://github.com/rmoschetti/CNDFinder)
