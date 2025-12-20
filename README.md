# Enriques surfaces with non-generic non-degeneracy

This repository contains the computational data supporting the paper *“Enriques surfaces with non-generic non-degeneracy”* by Riccardo Moschetti, Franco Rota, and Luca Schaffler.

The main contribution is the completion of the program initiated in [MRS24a]: computing the non-degeneracy invariant for the realizable $(\tau,\overline{\tau})$-generic Enriques surfaces (cf. [BS22]). The computations concern the two invariants

- **Fnd** (Fano non-degeneracy),
- **Mnd** (Mukai non-degeneracy).

Isotropic vectors and smooth rational curves are represented as vectors in $\mathrm{Num}(Y)$ with coordinates with respect to a fixed $E_{10}$-basis, as in [BS22].

## Software and file formats

- This project uses **SageMath 9.3**.
- Outputs are saved both as plain text (`.txt`) and as Python-pickled compressed files (`.pickle`) via the `pickle` module.

*Practical note.* For convenience, we recommend downloading (or cloning) the full repository, since several scripts assume the relative folder structure and local file paths. Users running SageMath through a notebook interface (e.g. Jupyter) can still use the .sage scripts by copying/pasting their contents into a notebook cell (or by running them via a terminal Sage session).

## About the [SB20] dataset and regenerated inputs

The original dataset from [SB20] is **not included** in this repository (to avoid duplicating external data). In particular, the JSON input files derived from [SB20] are *not* shipped.

A SageMath script `Parser.sage` (currently written to be largely Python-compatible by importing `from sage.all import *`) is provided to parse the [SB20] files (in GAP format) and generate the JSON inputs used by the computations. A reader who wants to cross-check [SB20] against the present computations can regenerate all JSON files via `Parser.sage` and compare the resulting databases.

**Exception:** in `4_Examples` there is one file included from the parsed data, because it is required to run the examples in that folder. This is documented in the corresponding folder `README.md`.

## Repository structure

For clarity, the project is split into independent parts. **Each folder is standalone** and contains its own `README.md` describing inputs/outputs and how to run the corresponding computations.

- **Part 1 (`1_LowerBounds`)**: Lower bounds for **Fnd** and **Mnd** (generalizing the data in [MRS24b] to these invariants).
- **Part 2 (`2_FiniteAutomorphisms`)**: Fano and Mukai polarizations for Enriques surfaces with finite automorphism group; uses the code [MRS22b] associated to [MRS22a].
- **Part 3 (`3_UpperBounds`)**: Upper bounds for **Fnd** and **Mnd**.
- **Part 4 (`4_Examples`)**: Additional data illustrating specific behaviours of **Fnd** and **Mnd** (Example 5.9 and Example 5.10 in the paper).

## Bibliography

[BS22] Brandhorst, Simon and Shimada, Ichiro. *Automorphism groups of certain Enriques surfaces*. Found. Comput. Math. 22 (2022).

[MRS22a] Moschetti, Riccardo; Rota, Franco; Schaffler, Luca. *A computational view on the non-degeneracy invariant for Enriques surfaces*. Experimental Mathematics. Published online: 29 August 2022.

[MRS22b] Moschetti, Riccardo; Rota, Franco; Schaffler, Luca. *SageMath code CndFinder*. Available at:
https://github.com/rmoschetti/CNDFinder

[MRS24a] Moschetti, Riccardo; Rota, Franco; Schaffler, Luca. *The non-degeneracy invariant of Brandhorst and Shimada’s families of Enriques surfaces*. Documenta Mathematica. Published online: 25 September 2024.

[MRS24b] Moschetti, Riccardo; Rota, Franco; Schaffler, Luca. Data for *“A computational view on the non-degeneracy invariant for Enriques surfaces”*. Available at:
https://github.com/rmoschetti/CND-TauTauBar-Enriques

[SB20] Shimada, Ichiro and Brandhorst, Simon. *Automorphism groups of certain Enriques surfaces — computational data*. Available at:
https://zenodo.org/record/4327019
