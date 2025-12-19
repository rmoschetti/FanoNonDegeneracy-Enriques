###############################################################################################
## STEP 0: Load data and specify an explicit length-9 isotropic sequence
###############################################################################################

import json
with open('158.txt') as json_file:
    data = json.load(json_file)

# Data conventions (as in the rest of the repository):
# - data['Ratstemp']: smooth rational curves as vectors in the E10 basis
# - data['Autrec']: automorphisms as matrices (row-vector convention v -> v * A)
# - data['Ells']: elliptic/isotropic vectors used for effectivity tests
Automorphisms = data['Autrec']

print("Working on Enriques TTB No:", data['no'])  # (typo fix: was "TBT")
print("Number of curves in Ratstemp: ", len(data['Ratstemp']))
print("Number of automorphisms: ", len(data['Autrec']))

R = [vector(rr) for rr in data['Ratstemp']]

# A chosen non-degenerate isotropic sequence of length 9 giving a Mukai polarization.
# Some vectors are half-sums of rational curves, as is standard for half-fibers.
Sequence = [
    1/2*(R[0] + R[2]),
    1/2*(R[2] + R[16]),
    1/2*(R[2] + R[2] * matrix(Automorphisms[0])^(-1)),
    1/2*(R[3] + R[4]),
    1/2*(R[3] + R[12]),
    1/2*(R[2] + R[14]),
    (R[2] + R[8]),
    1/2*(R[1] + R[6] + 2*R[7] + R[8] + R[15]),
    1/2*(R[5] + R[6] + R[8] + R[9] + 2*R[11])
]

# Associated Mukai polarization is (sum v_i)/2 for a length-9 isotropic sequence.
AssociatedMukai = 1/2 * sum(Sequence)
print("Associated Mukai polarization:", AssociatedMukai)

# Gram matrix of L10 in the chosen E10 basis (intersection pairing on Num(Y)).
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

###############################################################################################
## STEP 1: Compute the (-2)-classes orthogonal to the Mukai polarization
###############################################################################################

# We solve for integral vectors Ort such that:
#   Ort · v_i = 0 for all v_i in the sequence, and Ort^2 = -2.
var("a0,a1,a2,a3,a4,a5,a6,a7,a8,a9")
Ort = vector([a0,a1,a2,a3,a4,a5,a6,a7,a8,a9])

Eqn = []
for s in Sequence:
    Eqn.append(Ort * GramL10 * s == 0)
Eqn.append(Ort * GramL10 * Ort == -2)

Coefficients = solve(Eqn, [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9])

# In this example, the system has two integral solutions (opposite vectors).
Ort1 = Ort.subs(Coefficients[0])
Ort2 = Ort.subs(Coefficients[1])

print(f"The two integral (-2)-vectors orthogonal to the Mukai polarization are {Ort1} and {Ort2}")

###############################################################################################
## STEP 2: Show the orthogonal (-2)-classes are not effective (hence the Mukai is ample)
###############################################################################################

# If a (-2)-class D is effective, then it has nonnegative intersection
# with every nef isotropic class (in particular, with the elliptic half-fibers in Ells).
# Here we exhibit an elliptic vector with negative intersection, hence D is not effective.

for ee in data['Ells']:
    if vector(ee) * GramL10 * Ort1 < 0:
        print(f"Orthogonal vector {Ort1} not effective")
        break

for ee in data['Ells']:
    if vector(ee) * GramL10 * Ort2 < 0:
        print(f"Orthogonal vector {Ort2} not effective")
        break
