
# Part 2: Fano and Mukai polarizations for Enriques surfaces with finite automorphism group.

This part of the project contains the data relative to Fano and Mukai polarizations for the Enriques surfaces with finite automorphism group. It consists of the following folders and files.

 - FiniteAut_Fano (folder) - a folder containing 14 files (A text file and a pickle file for each of the seven Enriques surfaces with finite automorphism group). Each file contains the list of all the Fano polarizations expressed as a coefficient vector with respect to a fixed $E_{10}$-basis for $\mathrm{Num}(Y)$, as done in [BS22a].
 -  FiniteAut_Mukai (folder) - a folder containing 14 files (A text file and a pickle file for each of the seven Enriques surfaces with finite automorphism group). Each file contains the list of all the Mukai polarizations expressed as a coefficient vector with respect to a fixed $E_{10}$-basis for $\mathrm{Num}(Y)$, as done in [BS22a].
 - GenerateFanoAndMukai_FiniteAut.py - a SageMath code which uses [MRS22] to construct the lists provided in the above folders.

** Algorithm of GenerateFanoAndMukai_FiniteAut.py **

The strategy for finding the polarizations is to compute all possible isotropic sequences, selecting either the ones giving **Fano** or  **Mukai** polarizations. Moreover, we also determine their equivalence classes under the action of the surface’s automorphism group by using the generators of the automorphism groups of these Enriques listed in L10L26compdata.txt of [BS22b].


The computation proceeds as follows:
 - Retrieve input data for the chosen Enriques surface type from the dataset [BS22b] consisting of smooth rational curves and the generators of the automorphism groups. To check consistency between data sets, match the curves obtained from [BS22b] with the walls of an induced chamber for a suitable $(\tau,\overline{\tau})$-generic Enriques surface, as given in [SB20].
 - Use *CndFinder* of [MRS22] to obtain all saturated isotropic sequences from the curves.
 - Enumerate and extend these sequences using smooth rational curves to reach length 10 (resp. 9) for Fano (resp. Mukai) non-degeneracy.
 - Save raw results as lists of vectors.
 - Apply automorphisms to identify equivalence classes.
 - Save a list of one representative per equivalence class.


**Dependency note.** This script calls

The file `cnd_finder.py` is **not** included in this repository. It is part of *CndFinder* [MRS22] and can be obtained from the upstream source.
The folder `BSDATA` is **not** included in this repository. It can be generated starting by [SB20] and the file Parser.sage.
The file `L10L26compdata.txt` is **not** included in this repository. It can be found in the computational data of the paper [BS22b].

**Running the computation**
1. Set parameters at the beginning of the file:
   - `NumberTTB` – TTB number corresponding to the Enriques surface (see the code comments for the mapping from `K`-types).
   - `NumberK` – label for the surface type (`"I"`, `"II"`, …).
   - `ResultType` – `"Fano"` or `"Mukai"`.
2. Ensure the following are available:
   - BSData files: `BSData/E<TTB>.txt`
   - *CndFinder* script: `cnd_finder.py`
   - Automorphism group data: `L10L26compdata.txt` 
3. Run in SageMath:
   ```bash
   sage <script_name>.sage

**Bibliography**
[BS22a] Brandhorst, Simon and Shimada, Ichiro; Automorphism groups of certain Enriques surfaces, Found. Comput. Math. 22, 2022.

[BS22b] Brandhorst, Simon and Shimada, Ichiro; Borcherds’ Method for Enriques Surfaces, Michigan Math. J. 71, 2022. Computational data available at:
https://home.hiroshima-u.ac.jp/ichiro-shimada/ComputationData.html#L10L26

[MRS22] Riccardo Moschetti, Franco Rota, and Luca Schaffler. SageMath code CndFinder. Available at: https://github.com/rmoschetti/CNDFinder

[SB20] Shimada, Ichiro and Brandhorst, Simon; Automorphism groups of certain Enriques surfaces - computational data, https://zenodo.org/record/4327019# CND-TauTauBar-Enriques
