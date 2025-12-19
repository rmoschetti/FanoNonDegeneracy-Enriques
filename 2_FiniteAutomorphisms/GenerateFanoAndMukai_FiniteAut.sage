###############################################################################################
## STEP 0: Retrieve input data and auxiliary libraries
###############################################################################################

# Parameters controlling the run:
# - NumberTTB: which (tau, tau-bar)-generic entry to use as source of rational curves
# - NumberK:   finite-type label ("I", ..., "VII") used to select automorphism generators
# - ResultType: "Fano" or "Mukai" determines the target length and the divisibility test
NumberK    = "I"
NumberTTB  = "172"
ResultType = "Mukai"     # "Fano" or "Mukai"

# This is a possible correspondence between the labels of Enriques surfaces with finite automorphisms and the TTB invariants
# - KI   : TTB 172
# - KII  : TTB 184
# - KIII : TTB 171
# - KIV  : TTB 181
# - KV   : TTB 176
# - KVI  : TTB 179
# - KVII : TTB 182

# Data source: L10L26compdata.txt from [BS22b], read via libGAP.
libgap.Read("L10L26compdata.txt")


KData = None
for ikd in [0..17]:
    cand = libgap.Embs[ikd]
    if cand['NK'] == NumberK:
        KData = cand
        break
if KData is None:
    raise ValueError(f"Finite type K{NumberK} not found in L10L26compdata.txt")


# In this script, RatsTempList is the list of smooth rational curves, encoded as vectors in E10.
RatsTempList = [[int(y) for y in x] for x in KData['walls']]
print("Number of curves in ratstemp: ", len(RatsTempList))

# Requires CndFinder (cnd_finder.py) from [MRS22b]; see README file.
import logging
load('cnd_finder.py')
logging.getLogger().setLevel(logging.ERROR)


ListGenerators = []
for M in KData['generatorsOL10D']:
    ListGenerators.append(matrix(M))
   
import json
# Input data: expected to be a JSON file containing at least keys 'Walls' and 'no'
with open('BSData/E' + NumberTTB + '.txt') as json_file:
    data = json.load(json_file)

RatsTempListComparison = data['WallsV0']
if RatsTempListComparison == RatsTempList:
    print("The coefficients of the -2-curves are the same in the DB of TTB-generic enriques surfaces, and in the DB of irecs")
else:
    raise ValueError("Error: the coefficients of the curves are not the same")

###############################################################################################
## STEP 1: Use CndFinder to obtain saturated isotropic sequences
###############################################################################################

def GetBaseNum2(RatsTempList):
    """
    Extract a basis of Num(S) using smooth rational curves.
    Output basis is formatted to match cnd_finder.py conventions: vectors of length = #curves.
    """
    NumberOfRats = len(RatsTempList)

    # Find 10 curves whose matrix (as rows) has rank 10 (so they span Num over QQ).
    for S in Subsets(range(NumberOfRats), 10):
        ListOf10Rats = [RatsTempList[t] for t in S]
        if rank(matrix(ListOf10Rats)) == 10:
            break

    # Express the E10 basis in terms of the chosen 10 curves.
    E10BasisAsCombinationOfRats = matrix(ListOf10Rats)^(-1)
    E10BasisAsCombinationOfRatsAsList = [list(x) for x in list(E10BasisAsCombinationOfRats)]

    # Expand to vectors of length NumberOfRats, respecting cnd_finder.py indexing.
    FinalBasis = []
    for i in E10BasisAsCombinationOfRatsAsList:
        newBasisVector, CurveCounter = [], 0
        for j in range(NumberOfRats):
            if j in S:
                newBasisVector.append(i[CurveCounter])
                CurveCounter += 1
            else:
                newBasisVector.append(0)
        FinalBasis.append(newBasisVector)
    return FinalBasis

# Gram matrix for L10 in the chosen E10 basis
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
    [0, 0, 0, 0, 0, 0, 0, 0, 1, -2]
])

# Intersection matrix for the curves in RatsTempList (via GramL10)
IntersectionMatrixRatStemp = matrix(RatsTempList) * GramL10 * transpose(matrix(RatsTempList))
BasisNum = GetBaseNum2(RatsTempList)

# Run CndFinder; expected to return dict with key 'SaturatedSequences'
CndFinderResult = CndFinder(IntersectionMatrixRatStemp, BasisNum)

###############################################################################################
## STEP 2: Enumerate all subsets of saturated isotropic sequences (deduplicated as sets)
###############################################################################################

import itertools
ListWithAllPossibleSequences = [tuple([])]

def areSeqEqual(A, B):
    """Set-like equality for sequences (order ignored)."""
    if len(A) != len(B):
        return False
    for x in A:
        if x not in B:
            return False
    for x in B:
        if x not in A:
            return False
    return True

def IsInList(S):
    """Avoid duplicates among ListWithAllPossibleSequences (up to permutation)."""
    for x in ListWithAllPossibleSequences:
        if areSeqEqual(x, S):
            return True
    return False

for S in CndFinderResult['SaturatedSequences']:
    print("Check saturated sequences with lenght", S['Found cnd'])

    # Explore all sub-subsets of each saturated sequence returned by CndFinder.
    for TestSeq in S['ListSequencesVector']:
        for LL in [1..len(TestSeq)]:
            for NS in itertools.combinations(TestSeq, LL):
                if IsInList(list(NS)) == False:
                    ListWithAllPossibleSequences.append(NS)

print(f"Total number of isotropic sequences: {len(ListWithAllPossibleSequences)}")

###############################################################################################
## STEP 3: Extend each sequence using smooth rational curves, then filter by divisibility
###############################################################################################

def CheckNewVector(S, V, M):
    """
    Return True if V has intersection 1 with every vector in the existing isotropic sequence.
    Here M is the intersection matrix in the curve-basis conventions used by cnd_finder.py.
    """
    for g in S:
        for i in g:
            if intProduct(i, V, M)[0,0] != 1:
                return False
    return True

def AddVectors(A, B):
    return [A[i] + B[i] for i in range(len(A))]

def PrepareSearch(S, M, NumberOfNodes):
    # Currently unused (placeholder).
    return

def ExtendSaturatedSequence(Sequences, M, NumberOfNodes):
    """
    Recursively extend all candidate sequences by adding curve-classes (encoded as standard
    basis vectors in the curve-index space), until reaching the target length:
      - Fano: length 10
      - Mukai: length 9
    """
    if len(Sequences) == 0:
        return Sequences

    LengthCurrentSequence = len([item for sublist in Sequences[0]['list'] for item in sublist])

    if (ResultType == "Fano" and LengthCurrentSequence == 10):
        return Sequences
    if (ResultType == "Mukai" and LengthCurrentSequence == 9):
        return Sequences

    ExtendedSequences = []
    for S in Sequences:
        for g in range(len(S['list'])):
            for i in S['unused']:
                CurveAdd = i

                # Candidate extension: add the curve i (as a standard basis vector) to the last
                # vector in the g-th block.
                ToAdd = AddVectors(S['list'][g][-1], [1 if t == i else 0 for t in range(NumberOfNodes)])

                if CheckNewVector(S['list'], ToAdd, M) == True:
                    NewUnusedList = [x for x in S['unused'] if x != i]
                    NewSequence = [x[:] for x in S['list']]
                    NewSequence[g].append(ToAdd)

                    # Avoid duplicates (exact equality on the nested list structure).
                    newFound = True
                    for j in ExtendedSequences:
                        if NewSequence == j['list']:
                            newFound = False
                            break
                    if newFound == True:
                        ExtendedSequences.append({'list': NewSequence, 'unused': NewUnusedList, 'used': S['used'] + [CurveAdd]})

    return ExtendSaturatedSequence(ExtendedSequences, M, NumberOfNodes)

FinalResult = []

for iX, X in enumerate(ListWithAllPossibleSequences):
    ES = ExtendSaturatedSequence(
        [{'list': [[v] for v in X], 'unused': list(range(len(RatsTempList))), 'used': []}],
        IntersectionMatrixRatStemp,
        len(RatsTempList)
    )

    for SS in ES:
        Ext = [item for sublist in SS['list'] for item in sublist]
        A = matrix(QQ, Ext) * matrix(RatsTempList)  # convert from curve-index space back to E10 coordinates

        # Divisibility test selecting the valid polarizations:
        if (ResultType == "Fano"):
            DD = (vector(QQ, [1,1,1,1,1,1,1,1,1,1]) * A) / 3
        else:
            DD = (vector(QQ, [1,1,1,1,1,1,1,1,1]) * A) / 2

        if all(x in ZZ for x in DD.list()):
            FinalResult.append({'Original': X, 'Extended': Ext, 'Delta': DD, 'RatUsed': SS['used']})

DeltaList = [list(x['Delta']) for x in FinalResult]

###############################################################################################
## STEP 4: Save raw results
###############################################################################################

try:
    import pickle
except:
    from pickle5 import pickle

# Output files are written in the current working directory.
with open('K' + NumberK + ResultType + '.txt', 'w') as file:
    file.write(str(DeltaList))

with open('K' + NumberK + ResultType + '.pickle', 'wb') as file:
    pickle.dump(DeltaList, file, pickle.HIGHEST_PROTOCOL)

###############################################################################################
## STEP 5: Reduce modulo automorphisms (equivalence classes)
###############################################################################################

def LookForVector(x, ListOfVectors):
    """Return index of x in ListOfVectors, or -1 if not found."""
    for it, t in enumerate(ListOfVectors):
        if x == t:
            return it
    return -1

def ApplyAutToVector(A, x):
    """Apply automorphism A to vector x (row-vector convention)."""
    return (matrix(x) * matrix(A)).coefficients()

def ApplyAutToList(A, ListOfVectors, Equivalence):
    """
    Update Equivalence classes using automorphism A: whenever AV is present in ListOfVectors,
    merge the classes of v and AV.
    """
    for iv, v in enumerate(ListOfVectors):
        AV = ApplyAutToVector(A, v)
        iAV = LookForVector(AV, ListOfVectors)

        if iAV != -1:
            PosEquivalenceOld = iv
            PosEquivalenceNew = iAV
            if (Equivalence[PosEquivalenceOld] != Equivalence[PosEquivalenceNew]):
                print(f"List updated: substituted class at position {PosEquivalenceOld} (currently {Equivalence[PosEquivalenceOld]}) with the one at position {PosEquivalenceNew} (currently {Equivalence[PosEquivalenceNew]})")

            # Merge classes by relabeling all entries of the old class.
            for ix, _ in enumerate(Equivalence):
                if Equivalence[ix] == Equivalence[PosEquivalenceOld]:
                    Equivalence[ix] = Equivalence[PosEquivalenceNew]

ListVectors = [fd['Delta'].coefficients() for fd in FinalResult]
Equivalence = list(range(len(ListVectors)))

for A in ListGenerators:
    ApplyAutToList(A, ListVectors, Equivalence)

###############################################################################################
## STEP 6: Save reduced results (one representative per equivalence class)
###############################################################################################

DeltaAut = []
for eqcl in list(set(Equivalence)):
    DeltaAut.append(ListVectors[eqcl])

with open('K' + NumberK + ResultType + 'Aut.pickle', 'wb') as file:
    pickle.dump(DeltaAut, file, pickle.HIGHEST_PROTOCOL)

with open('K' + NumberK + ResultType + 'Aut.txt', 'w') as file:
    file.write(str(DeltaAut))
