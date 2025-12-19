from sage.all import *

# Compatibility: SageMath 9.3 ships with pickle, but the fallback covers older environments.
try:
    import pickle
except:
    from pickle5 import pickle


"""
This script validates the serialized lower-bound data stored in `LowerBoundsFanoMukai.pickle`.

For each Enriques surface index in 1..184, it checks:
  - internal consistency of stored sequence lengths vs. reported invariants (Fnd/Mnd),
  - that each stored isotropic sequence vector equals the result of applying the recorded
    automorphisms to the recorded starting isotropic vector,
  - divisibility constraints in the extreme cases (Fnd=10, Mnd=9, Fnd=9),
  - agreement of Fano/Mukai data in the small-length regime (<= 8) when Fnd==Mnd.

It prints a per-index log to stdout. It does not write output files.
"""

with open('LowerBoundsFanoMukai.pickle', 'rb') as f:
    DicLowerBounds = pickle.load(f)


def CheckLength(item):
    """Sanity-check that sequence lengths match the stored invariant values."""
    if item['FanoPolarization']['Fnd'] != len(item['FanoPolarization']['IsotropicSequenceDescription']):
        return False
    if item['FanoPolarization']['Fnd'] != len(item['FanoPolarization']['IsotropicSequenceVectors']):
        return False
    if item['MukaiPolarization']['Mnd'] != len(item['MukaiPolarization']['IsotropicSequenceDescription']):
        return False
    if item['MukaiPolarization']['Mnd'] != len(item['MukaiPolarization']['IsotropicSequenceVectors']):
        return False
    return True


def CheckComputation(item, Type='FanoPolarization'):
    """
    Recompute each isotropic vector sequence from its "description".

    Convention expected in the pickle:
      - Description[0] indexes an entry in item['IsotropicVectors']
      - Description[1:], each entry indexes an automorphism in item['Automorphisms']
      - The sequence vector is obtained by applying these automorphisms in order.
    """
    for iSeq, Seq in enumerate(item[Type]['IsotropicSequenceVectors']):
        Description = item[Type]['IsotropicSequenceDescription'][iSeq]

        # Start from the isotropic vector specified by Description[0].
        StartingVector = vector(item['IsotropicVectors'][Description[0]])

        # Apply the automorphisms encoded by Description[1:].
        AutomorphismList = [matrix(item['Automorphisms'][Description[t]]) for t in range(1, len(Description))]
        for mAut in AutomorphismList:
            # Note: row-vector convention (v * A).
            StartingVector = StartingVector * mAut

        # Stored vectors are lists; StartingVector is a Sage vector.
        if list(StartingVector) != Seq:
            print(iSeq, list(StartingVector), Seq)
            return False

    return True


def CheckIsotropic(seq):
    """
    Check that `seq` is an isotropic sequence for the lattice with Gram matrix GramL10:
      - self-pairings = 0
      - pairings between distinct vectors = 1

    (Currently unused in the main loop; kept as an auxiliary check.)
    """
    GramL10 = matrix([
        [-2, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, -2, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, -2, 1, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, -2, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, -2, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, -2, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, -2, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, -2, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, -2, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, -2],
    ])

    for iv in range(len(seq)):
        for iw in range(iv, len(seq)):
            pairing = (matrix(seq[iv]) * GramL10 * transpose(matrix(seq[iw])))[0,0]
            if iv==iw and pairing != 0:
                return False
            if iv!=iw and pairing != 1:
                return False
    return True


def CheckDivisibility(seq, Den=3):
    """Check whether (sum(seq) / Den) lies in ZZ^10 (componentwise)."""
    Res = sum([vector(QQ, t) for t in seq]) / Den
    return all(t in ZZ for t in Res)


for index in range(1, 185):
    print(f"Check Enriques no. {index}")

    if index not in DicLowerBounds:
        print(f"\tError: index not in the dictionary")
        continue

    if DicLowerBounds[index]['exists'] == False:
        print(f"\tNothing to check: the Enriques does not exists.\n")
        continue

    if CheckLength(DicLowerBounds[index]) == False:
        print(f"\tError: lenght of sequences does not coincide with Fnd / Mnd reported in the dictionary")
        continue

    if CheckComputation(DicLowerBounds[index], Type='FanoPolarization') == False:
        print(f"\tError: computations for the IsotropicSequenceVectors from IsotropicSequenceDescription of the Fano Polarization wrong")
        continue

    if CheckComputation(DicLowerBounds[index], Type='MukaiPolarization') == False:
        print(f"\tError: computations for the IsotropicSequenceVectors from IsotropicSequenceDescription of the Mukai Polarization wrong")
        continue


    if CheckIsotropic(DicLowerBounds[index]['FanoPolarization']['IsotropicSequenceVectors'])==False:
        print(f"\tError: isotropicity condition failed for the Fano polarization sequence")
        continue

    if CheckIsotropic(DicLowerBounds[index]['MukaiPolarization']['IsotropicSequenceVectors'])==False:
        print(f"\tError: isotropicity condition failed for the Mukai polarization sequence")
        continue


    # Divisibility constraints in the extremal-length cases.
    if DicLowerBounds[index]['FanoPolarization']['Fnd'] == 10:
        if CheckDivisibility(DicLowerBounds[index]['FanoPolarization']['IsotropicSequenceVectors'], Den=3) == False:
            print(f"\tError: the sequence for the Fano polarization gives a vector not divisible by three")
            continue

    if DicLowerBounds[index]['MukaiPolarization']['Mnd'] == 9:
        if CheckDivisibility(DicLowerBounds[index]['MukaiPolarization']['IsotropicSequenceVectors'], Den=2) == False:
            print(f"\tError: the sequence for the Mukai polarization gives a vector not divisible by two")
            continue

    if DicLowerBounds[index]['FanoPolarization']['Fnd'] == 9:
        if CheckDivisibility(DicLowerBounds[index]['FanoPolarization']['IsotropicSequenceVectors'], Den=2) == True:
            print(f"\tError: the sequence for the Fano polarization of length 9 gives a vector divisible by two")
            continue

    # If Fnd==Mnd and <=8, the script expects the *same* sequence data for Fano and Mukai.
    if DicLowerBounds[index]['FanoPolarization']['Fnd'] == DicLowerBounds[index]['MukaiPolarization']['Mnd']:
        if DicLowerBounds[index]['FanoPolarization']['Fnd'] <= 8:
            if (
                DicLowerBounds[index]['FanoPolarization']['IsotropicSequenceDescription'] != DicLowerBounds[index]['MukaiPolarization']['IsotropicSequenceDescription']
                or DicLowerBounds[index]['FanoPolarization']['IsotropicSequenceVectors'] != DicLowerBounds[index]['MukaiPolarization']['IsotropicSequenceVectors']
            ):
                print(f"\tError: Fano and Mukai of length <= 8 non coinciding")
                continue

    print(f"\tCondluded successfully without errors.\n")
