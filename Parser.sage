"""
Parser.sage

Reads the GAP dataset `Enrs.txt` from [SB20] (Zenodo: 10.5281/zenodo.4327019) via libGAP
and exports, for each Enriques surface in `Enrs`, a JSON file containing selected fields.

Each output file is named `E<no>.txt`, where `<no>` is the surface number used in [SB20],
and is written into the folder `BSData/`.

Main exported keys:
  - no, irec
  - Ratstemp (rational curves as vectors)
  - Ells (isotropic vectors as in the source list Ellstemp)
  - Autrec (selected automorphisms)
  - WallsV0, VolumeIndex
  - ChamberList (adjacency and chamber automorphisms extracted from Autrec.V0)

Note: the JSON serialization uses `default=convert` to coerce GAP integer/list types to
plain Python types.
"""

import os
import json


# Output folder for parsed data
OUT_DIR = "BSData"
os.makedirs(OUT_DIR, exist_ok=True)


def convert(o):
    """
    Convert libGAP objects to JSON-serializable Python objects.

    Currently supports:
      - GapElement_Integer -> int
      - GapElement_List of integers -> list[int]
      - GapElement_List of lists of integers -> list[list[int]]
    """
    if isinstance(o, GapElement_Integer):
        return int(o)

    if isinstance(o, GapElement_List):
        if len(o) == 0:
            return []
        if isinstance(o[0], GapElement_Integer):
            return [int(x) for x in o]
        if isinstance(o[0], GapElement_List):
            return [[int(y) for y in x] for x in o]

    try:
        return str(o)
    except Exception:
        raise TypeError(f"Object of type {type(o)} is not JSON-convertible by convert().")


# Read the file "Enrs.txt" produced by [SB20].
# If GAP runs out of memory, you may need to increase the GAP memory pool; see:
# https://doc.sagemath.org/html/en/reference/interfaces/sage/interfaces/gap.html#sage.interfaces.gap.set_gap_memory_pool_size
libgap.Read("Enrs.txt")
Enrs = libgap.Enrs
LengthListEnrs = len(Enrs)
print("Importing", LengthListEnrs, "Enriques")


for i in range(LengthListEnrs):
    NumberInTheList = Enrs[i]["no"]
    print("\tImporting the", i, "-th element, corresponding to number", NumberInTheList, "in the list.")

    # --- Rational curves: Rats.Ratstemp[*].rat.ratY
    ListRatsTemp = []
    RatsListForThisEnriques = Enrs[i]["Rats"]["Ratstemp"]
    for j in range(len(RatsListForThisEnriques)):
        ListRatsTemp.append(RatsListForThisEnriques[j]["rat"]["ratY"])

    # --- Isotropic vectors: Ells.Ellstemp[*].ell
    ListEllsTemp = []
    EllsListForThisEnriques = Enrs[i]["Ells"]["Ellstemp"]
    for j in range(len(EllsListForThisEnriques)):
        ListEllsTemp.append(EllsListForThisEnriques[j]["ell"])

    # --- Automorphisms: Autrec.HHH[*].gY
    ListAutRec = []
    AutrecListForThisEnriques = Enrs[i]["Autrec"]["HHH"]
    for j in range(len(AutrecListForThisEnriques)):
        ListAutRec.append(AutrecListForThisEnriques[j]["gY"])

    # --- Additional invariants from SXrec and global name irecname
    ListWallsV0 = Enrs[i]["SXrec"]["walls"]
    VolumeIndex = Enrs[i]["SXrec"]["volumeindex"]
    IrecName = Enrs[i]["irecname"]

    # --- Chamber data extracted from Autrec.V0
    ChamberList = []
    AutRecV0 = Enrs[i]["Autrec"]["V0"]
    for V in AutRecV0:
        NewChamber = {
            "pos": V["pos"],
            "from": V["from"],
            "adjpos": V["adjpos"],
            "adjrecs": [],
            "autchamb": [],
        }

        for w in V["adjrecs"]:
            if (w["israt"] == False and w["isnew"] == True):
                NewChamber["adjrecs"].append({
                    "wall": w["wall"],
                    "indexTo": w["positionInV0"],
                })

        for A in V["autcham"]:
            NewChamber["autchamb"].append(A["gY"])

        ChamberList.append(NewChamber)

    NewEnriques = {
        "no": NumberInTheList,
        "irec": IrecName,
        "Ratstemp": ListRatsTemp,
        "Ells": ListEllsTemp,
        "Autrec": ListAutRec,
        "WallsV0": ListWallsV0,
        "VolumeIndex": VolumeIndex,
        "ChamberList": ChamberList,
    }

    out_path = os.path.join(OUT_DIR, "E" + str(NumberInTheList) + ".txt")
    with open(out_path, "w") as fp:
        json.dump(NewEnriques, fp, default=convert)
