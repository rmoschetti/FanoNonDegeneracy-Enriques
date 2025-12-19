TTB_Mukai_folder_path = 'TTB_Mukai'
FiniteAut_folder_path = '../2_FiniteAutomorphisms/FiniteAut_Mukai'

# Hardcoded correspondence for the Mukai upper-bound checks:
# - CorrespondingK: finite-automorphism type used to import the expected set of Mukai polarizations
# - WallSecondChamber: if not None, also include polarizations reflected across this wall
DataToCheck = {
    '143Mukai.pickle': {'CorrespondingK': '2', 'WallSecondChamber': [2, 1, 2, 3, 3, 3, 2, 2, 1, 0]},
    '144Mukai.pickle': {'CorrespondingK': '3', 'WallSecondChamber': None},
    '171Mukai.pickle': {'CorrespondingK': '3', 'WallSecondChamber': None},
}

try:
    import pickle
except:
    from pickle5 import pickle


def getVectorFromDescription(itemDescription, curvesList, automorphismList):
    """
    Evaluate a description of the form:
      ['R<i>'] or ['R<i>', 'H<j>', ...]
    meaning: start from curvesList['R<i>'] and right-multiply by the listed automorphisms.

    Conventions:
      - curvesList: dictionary of vectors (E10 coordinates)
      - automorphismList: dictionary of matrices
      - row-vector convention: v -> v * A
    """
    StartingVector = vector(curvesList[itemDescription[0]])
    AutomorphismList = [matrix(automorphismList[itemDescription[t]]) for t in range(1, len(itemDescription))]
    for mAut in AutomorphismList:
        StartingVector = StartingVector * mAut
    return list(StartingVector)


# Gram matrix for L10 (E10 basis), used for intersection computations.
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


def checkIntersection(mp, ListRat, RatsUsed, AutUsed):
    """
    For a fixed Mukai polarization mp, verify that every curve in ListRat has intersection 0 with mp.
    """
    for rat in ListRat:
        VectCurve = getVectorFromDescription(rat, RatsUsed, AutUsed)
        if (vector(mp) * GramL10 * vector(VectCurve) != 0):
            return False
    return True


def myintProduct(x, y):
    return vector(x) * GramL10 * vector(y)


def Reflection(V, W):
    """Reflect V across the wall W (in the L10 pairing)."""
    C = -2 * myintProduct(V, W) / myintProduct(W, W)
    return [V[i] + C * W[i] for i in range(10)]


def checkMukaiPolarization(thisData, CorrespondingK, WallSecondChamber):
    """
    Compare the set of polarizations (keys of ListDelta) with those imported from Part 2
    for the finite-automorphism type K{CorrespondingK}. Optionally include the reflected
    set if a second chamber is required.
    """
    SetCurrentMp = thisData['ListDelta'].keys()
    with open(FiniteAut_folder_path + '/K' + CorrespondingK + '_Mukai.pickle', 'rb') as f:
        ListMpImportedFromKondo = pickle.load(f)

    SetMpImportedFromKondo = set([tuple(v) for v in ListMpImportedFromKondo])

    if WallSecondChamber != None:
        for imp in ListMpImportedFromKondo:
            SetMpImportedFromKondo.add(tuple(Reflection(imp, WallSecondChamber)))

    return SetCurrentMp == SetMpImportedFromKondo


def CheckIsotropic(seq):
    """
    Check isotropicity for a sequence {v_i} in L10:
      v_i·v_i = 0 and v_i·v_j = 1 for i != j.
    """
    for iv in range(len(seq)):
        for iw in range(iv, len(seq)):
            pairing = (matrix(seq[iv]) * GramL10 * transpose(matrix(seq[iw])))[0, 0]
            if iv == iw and pairing != 0:
                return False
            if iv != iw and pairing != 1:
                return False
    return True


def CheckNinthCurve(Seq, VectorCurve, NinthVector):
    """
    Check that NinthVector can be written as (VectorCurve + s) for some s in Seq.
    (Here VectorCurve is the vector of '9th-RatCurve'.)
    """
    for s in Seq:
        if (vector(s) + vector(VectorCurve) == vector(NinthVector)):
            return True
    return False


for key in DataToCheck:
    print(f"Checking TTB {key}")
    with open(TTB_Mukai_folder_path + '/' + key, 'rb') as f:
        thisData = pickle.load(f)

    if checkMukaiPolarization(thisData, DataToCheck[key]['CorrespondingK'], DataToCheck[key]['WallSecondChamber']) == False:
        print(f"\tError: Set of Mukai polarizations not corresponding")
        continue

    skip_key = False

    for mp in thisData['ListDelta']:
        entry = thisData['ListDelta'][mp]

        if entry['DataType'] == 'Two curves':
            if len(entry['Data']) != 2:
                print(f"\tError: Number of curves different than two")
                skip_key = True
                break

            if checkIntersection(mp, entry['Data'], thisData['RatsUsed'], thisData['AutUsed']) == False:
                print(f"\tError: Intersection non zero between curves and Mukai polarization")
                skip_key = True
                break

        elif entry['DataType'] == 'Sequence':
            # Build the saturated length-8 sequence from stored descriptions (elliptic data + automorphisms),
            # then append the stored 9th vector.
            Sequence = [
                getVectorFromDescription(itemDescription, thisData['EllipticConfUsed'], thisData['AutUsed'])
                for itemDescription in entry['Data']['SaturatedEightSequence']
            ]
            Sequence = Sequence + [entry['Data']['9th-vector']]

            if CheckIsotropic(Sequence) == False:
                print(f"\tError: Sequence non isotropic")
                skip_key = True
                break

            # Mukai polarization associated with a length-9 isotropic sequence: (sum v_i) / 2
            ComputedMukai = sum([vector(QQ, t) for t in Sequence]) / 2
            if all(t in ZZ for t in ComputedMukai) == False:
                print(f"\tError: Mukai polarization not divisible by two")
                skip_key = True
                break

            if tuple(ComputedMukai) != mp:
                print(f"\tError: Computed Mukai polarization not corresponding to key")
                skip_key = True
                break

            NinthCurveVector = getVectorFromDescription(entry['Data']['9th-RatCurve'], thisData['RatsUsed'], thisData['AutUsed'])
            if CheckNinthCurve(Sequence, NinthCurveVector, entry['Data']['9th-vector']) == False:
                print(f"\tError: Ninth vector not corresponding to Ninth curve")
                skip_key = True
                break

        else:
            print(f"\tError: Data type not recognized")
            skip_key = True
            break


    if skip_key == False:
        print(f"\tConcluded successfully without errors.\n")
