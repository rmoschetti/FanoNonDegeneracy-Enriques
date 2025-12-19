# SingularMukaiModel.sage
#
# Goal (interactive): exhibit a Mukai polarization on the finite-automorphism Enriques surface KVI
# (here encoded as K6) which is non-degenerate (length-9 isotropic sequence) but singular,
# i.e. it has a (-2)-curve in its orthogonal complement that is contracted (intersection 0).

# Requires CndFinder (cnd_finder.py) from [MRS22b]; see README file.
load('cnd_finder.py')
logging.getLogger().setLevel(logging.ERROR)

# Intersection matrix of the curve configuration (20 curves) for K6.
# This is the input to CndFinder: it encodes pairings of the curves in the curve-index basis.
IntersectionMatrixK6 = matrix([
     [-2, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2],
     [1, -2, 1, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 1, -2, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0],
     [0, 0, 1, -2, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0],
     [0, 0, 0, 1, -2, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0],
     [1, 0, 0, 0, 1, -2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0],
     [0, 1, 0, 0, 1, 0, -2, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 1, -2, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 1, 0, 0, 1, 0, 1, -2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0],
     [1, 0, 0, 1, 0, 0, 0, 1, 0, -2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 2, 0, 0, -2, 1, 1, 0, 0, 0, 1, 1, 1, 1],
     [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 1, 1, 0, 1, 1, 0],
     [0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 1, -2, 0, 1, 1, 1, 0, 0, 1],
     [0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, -2, 1, 1, 1, 1, 1, 1],
     [0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 1, 1, -2, 1, 0, 1, 0, 1],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 1, 1, 1, -2, 1, 0, 1, 0],
     [0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, -2, 0, 1, 1],
     [0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, -2, 1, 1],
     [0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, -2, 0],
     [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, -2]
])

# A Num-basis for CndFinder: 10 vectors (rank 10) expressed in the curve-index space (length 20).
# Here they are chosen as standard basis vectors corresponding to the listed curves.
BasisNumK6 = [
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # R1
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # R2
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # R3
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # R4
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # R5
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # R7
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # R11
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],  # R12
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],  # R14
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],  # R17
]

# Expected execution time (as recorded by the authors): ~212 seconds on an i7-7700HQ.
CndResult = CndFinder(IntersectionMatrixK6, BasisNumK6)

for sat in CndResult['SaturatedSequences']:
    # We select saturated sequences of length 9 of a specific "type" produced by CndFinder.
    # For these, non-degeneracy is built into the construction (CndFinder classification).
    if sat['Found cnd'] == 9 and sat['Type'] == '6 x (1 A1HF + 1 A2F + 1 A5F), 3 x (1 A3F + 1 D5F)':
        for iseq, seq in enumerate(sat['ListSequencesVector']):

            # seq is a list of 9 vectors in the curve-index space (length 20).
            SeqVector = matrix(QQ, seq)

            # Mukai polarization associated to a length-9 isotropic sequence: (sum v_i)/2.
            MukaiPol = (matrix(QQ, [1,1,1,1,1,1,1,1,1]) * SeqVector) / 2

            # Express the Mukai polarization in the chosen Num-basis (BasisNumK6) and check integrality.
            CoefficientsMukaiPol = vector(MukaiPol * IntersectionMatrixK6 * transpose(matrix(BasisNumK6)))
            if all(t in ZZ for t in CoefficientsMukaiPol) == False:
                print("Error: working on a length-9 sequence which does not give a Mukai polarization")
                continue

            # Search for a (-2)-curve among the 20 curves which is contracted by this Mukai polarization,
            # i.e. has intersection 0 with MukaiPol.
            for irat in range(20):
                Vrat = [1 if i == irat else 0 for i in range(20)]  # curve R_{irat+1} in curve-index space
                Int = matrix(QQ, Vrat) * IntersectionMatrixK6 * transpose(MukaiPol)

                # Int is a 1x1 matrix; Int == 0 detects orthogonality (contracted curve).
                if (Int == 0):
                    print(f"The Mukai polarization with coefficients {CoefficientsMukaiPol} is zero against the curve R_{irat+1}\n")
