# Part 1: Lower bounds for the Fnd and the Mnd invariants

This folder generalizes the data contained in the project [MRS24b] to the invariants **Fnd** and **Mnd**. For each realizable $(\tau,\overline{\tau})$-generic Enriques surface we provide two sequences of isotropic vectors:
- one giving a lower bound for the **Fano non-degeneracy** invariant (**Fnd**),
- one giving a lower bound for the **Mukai non-degeneracy** invariant (**Mnd**).

## Summary of the files contained in this folder

- `LowerBoundsFanoMukai.pickle` (optionally also `LowerBoundsFanoMukai.txt`): a dictionary containing the isotropic sequences giving the lower bounds for **Fnd** and **Mnd**.
  - The `.pickle` file is the primary serialized format used by the checker script.
  - The `.txt` file, when present, stores the same data in JSON text format and can be imported in SageMath.
- `CheckDataLowerBounds.sage`: a SageMath script which loads `LowerBoundsFanoMukai.pickle` and performs automatic checks to ensure that the sequences stored in `LowerBoundsFanoMukai.*` are consistent and satisfy the required conditions.

## Data structure of `LowerBoundsFanoMukai.*`

`LowerBoundsFanoMukai.*` contains a dictionary indexed by the number of $(\tau,\overline{\tau})$-generic Enriques surfaces. The corresponding value for each key is a dictionary structured as follows:

- `no`: the number identifying the $(\tau,\overline{\tau})$-generic Enriques surface as listed in [Table 1, BS22]. This typically coincides with the key of the outer dictionary.
- `exists`: a boolean flag; if `False`, the corresponding surface is marked as non-existing and is skipped by the checker.
- `FanoPolarization`: a dictionary containing the sequence giving a lower bound for a Fano polarization:
  - `Fnd`: the value of **Fnd** corresponding to this sequence (it coincides with the length of both `IsotropicSequenceDescription` and `IsotropicSequenceVectors`).
  - `IsotropicSequenceDescription`: a sequence of *descriptions* referring to the lists `IsotropicVectors` and `Automorphisms` of the current Enriques surface (see examples below).
  - `IsotropicSequenceVectors`: the resulting sequence, as a list of vectors expressed in the $E_{10}$-basis for $\mathrm{Num}(Y)$.
- `MukaiPolarization`: a dictionary containing the sequence giving a lower bound for a Mukai polarization:
  - `Mnd`: the value of **Mnd** corresponding to this sequence (it coincides with the length of both `IsotropicSequenceDescription` and `IsotropicSequenceVectors`).
  - `IsotropicSequenceDescription`: as above.
  - `IsotropicSequenceVectors`: as above.
- `IsotropicVectors`: dictionary containing the coefficients of the isotropic vectors in `Ellstemp` used for computing the lower bounds. The indices reflect the original ordering of the vectors in the data of [SB20].
- `Automorphisms`: dictionary containing the matrices of the automorphisms in `HHH` used for computing the lower bounds. The indices reflect the original ordering of the data of [SB20].

### Description format

Each entry of `IsotropicSequenceDescription` is a list describing how to obtain the corresponding vector in `IsotropicSequenceVectors`.

Example:
> `['I27']`

This denotes the isotropic vector named `I27` in the dictionary `IsotropicVectors` (coming from `Ellstemp` in [SB20]).

Example:
> `['I37', 'H4']`

This denotes the vector obtained from `IsotropicVectors['I37']` by right-multiplication by the automorphism `Automorphisms['H4']` (row-vector convention $v \mapsto v \cdot H4$).

## Running `CheckDataLowerBounds.sage`

Run from within this folder:
- `sage CheckDataLowerBounds.sage`

The script iterates over indices `1..184` and prints a per-index log. Depending on the entry, it either reports that nothing is checked (if `exists == False`), or it stops at the first failed check and prints an error message, or it reports success.

For each isotropic sequence $(v_i)$, the script checks:

- That the length of the sequences coincides with the claimed lower bound for the corresponding invariant.
- The correctness of the reconstruction from `IsotropicSequenceDescription` to `IsotropicSequenceVectors`.
- The condition $v_i \cdot v_j = 1- \delta_{ij}$.
- If `Fnd = 10`, that $\sum v_i$ is divisible by 3.
- If `Mnd = 9`, that $\sum v_i$ is divisible by 2.
- If `Fnd = 9`, that $\sum v_i$ is **not** divisible by 2 (hence it can be extended to a Fano polarization).
- If `Fnd = Mnd \le 8`, that the data for the Fano and Mukai polarizations coincide.

## Bibliography

[BS22] Brandhorst, Simon and Shimada, Ichiro; *Automorphism groups of certain Enriques surfaces*, Found. Comput. Math. 22 (2022).

[MRS24b] Riccardo Moschetti, Franco Rota, and Luca Schaffler. Data for *"A computational view on the non-degeneracy invariant for Enriques surfaces"*. Available at: https://github.com/rmoschetti/CND-TauTauBar-Enriques.

[SB20] Shimada, Ichiro and Brandhorst, Simon; *Automorphism groups of certain Enriques surfaces - computational data*, https://zenodo.org/record/4327019
