# Part 3: Upper bounds for Fnd and Mnd

This folder contains the computational data and consistency checks used to establish the **upper bounds** for **Fnd** and **Mnd** in the $(\tau,\overline{\tau})$-generic cases appearing in the main theorems.

## Summary of the contents

- `TTB_Fano/` - data used to prove an upper bound for **Fnd** for the following $(\tau,\overline{\tau})$-generic Enriques surfaces:
  $Y_{84}, Y_{85}, Y_{121}, Y_{122}, Y_{123}, Y_{143}, Y_{144}, Y_{158}, Y_{159}, Y_{171}, Y_{176}$.
- `TTB_Mukai/` - data used to prove an upper bound for **Mnd** for:
  $Y_{143}, Y_{144}, Y_{171}$.
- `checkUpperBoundsFnd.sage` - performs consistency checks on the data contained in `TTB_Fano/`.
- `checkUpperBoundsMnd.sage` - performs consistency checks on the data contained in `TTB_Mukai/`.

## Data format in `TTB_Fano/`

For each surface $Y_n$ listed above, the corresponding file contains a dictionary with the following keys:

- `RatsUsed`: dictionary containing the coefficient vectors of the **smooth rational curves** (from `Ratstemp`) that are used in the verification. The indices reflect the ordering in the data of [SB20].
- `AutUsed`: dictionary containing the matrices of the **automorphisms** (from `HHH`) that are used in the verification. The indices reflect the ordering in the data of [SB20].
- `ListDelta`: a dictionary indexed by fano polarizations produced in Part 2. The corresponding value is a list of smooth rational curves written as elements of `Ratstemp`, possibly multiplied by automorphisms in `HHH`. For example:
  `(30, 18, 38, 60, 55, 50, 38, 27, 16, 5): [['R1'], ['R3', 'H2']]`

## Checks performed by `checkUpperBoundsFnd.sage`

For each file in `TTB_Fano/`, the script checks:

- For each element of `ListDelta`, the intersection product of each listed curve with the corresponding Fano polarization is zero.
- Coverage of the polarizations that must be checked for the corresponding finite-automorphism type (this requires data from Part 2).
- For each Fano polarization `fp`, the list `ListDelta[fp]` contains **at least** $10-\mathrm{ExpectedFnd}$ curves with intersection zero with `fp`.

The script uses a hardcoded dictionary `DataToCheck`, containing:
- a correspondence between the $(\tau,\overline{\tau})$-generic surfaces under consideration and the finite-automorphism Enriques types (to be cross-checked against Table 1.1 of [BS22a] and Table 1.2 of [BS22b]);
- the wall needed to move to the second chamber (either `None`, or a wall vector defining the required reflection);
- the claimed upper-bound value for **Fnd**.

## Data format in `TTB_Mukai/`

For each surface $Y_{143}, Y_{144}, Y_{171}$, the corresponding file contains a dictionary with three keys:

- `RatsUsed`: dictionary of smooth rational curves used in the verification (as above).
- `AutUsed`: dictionary of automorphisms used in the verification (as above).
- `ListDelta`: a dictionary indexed by Mukai polarizations produced in Part 2. For each key, the value is a dictionary containing:
  - `DataType`: either `Two curves` or `Sequence`;
  - `Data`: data used to certify the upper bound for the corresponding Mukai polarization, depending on `DataType`.

If `DataType = 'Two curves'`, then `Data` is a list of exactly two smooth rational curves (each expressed via `Ratstemp` and possibly `HHH`) which satisfy the required intersection condition(s), similarly to the `TTB_Fano/` case.

If `DataType = 'Sequence'`, then `Data` is a dictionary containing:
- `SaturatedEightSequence`: a non-degenerate isotropic sequence of length 8 (built from `Ellstemp` and `HHH`);
- `9th-vector`: a vector extending it to a length-9 isotropic sequence;
- `9th-RatCurve`: a rational curve (expressed via `Ratstemp` and possibly `HHH`) such that `9th-vector` equals `9th-RatCurve` plus one element of `SaturatedEightSequence`.

## Checks performed by `checkUpperBoundsMnd.sage`

The script is analogous to `checkUpperBoundsFnd.sage`, with additional checks depending on `DataType`.

- Coverage of the polarizations that must be checked for the corresponding finite-automorphism type (this requires data from Part 2, in `FiniteAut_Mukai`).

If `DataType = 'Two curves'`:
- check that `Data` contains exactly two elements;
- verify that each of the two curves has intersection zero with the corresponding Mukai polarization.

If `DataType = 'Sequence'`:
- check $v_i \cdot v_j = 1 - \delta_{ij}$ for `SaturatedEightSequence` and for the extended length-9 sequence obtained by adding `9th-vector`;
- check that the resulting length-9 sequence yields the intended Mukai polarization;
- check that `9th-vector` is obtained as `9th-RatCurve` plus one element of `SaturatedEightSequence`.

If you read until here we owe you a drink of your choice.

## Running the scripts

From within this folder:
- `sage checkUpperBoundsFnd.sage`
- `sage checkUpperBoundsMnd.sage`

## Bibliography

[BS22a] Brandhorst, Simon and Shimada, Ichiro. *Automorphism groups of certain Enriques surfaces*. Found. Comput. Math. 22 (2022).

[BS22b] Brandhorst, Simon and Shimada, Ichiro. *Borcherds’ method for Enriques surfaces*. Michigan Math. J. 71 (2022).

[SB20] Shimada, Ichiro and Brandhorst, Simon. *Automorphism groups of certain Enriques surfaces — computational data*. Available at:
https://zenodo.org/record/4327019
