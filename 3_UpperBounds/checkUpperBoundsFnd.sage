import os
TTB_Fano_folder_path = 'TTB_Fano'
FiniteAut_folder_path = '../2_FiniteAutomorphisms/FiniteAut_Fano'

DataToCheck={
    '84Fano.pickle' : {'CorrespondingK':'3', 'WallSecondChamber': [2, 1, 2, 3, 2, 1, 0, 0, 0, 0], 'ExpectedFnd':9},
    '85Fano.pickle' : {'CorrespondingK':'5', 'WallSecondChamber':None,                            'ExpectedFnd':7},
    '121Fano.pickle': {'CorrespondingK':'3', 'WallSecondChamber':None,                            'ExpectedFnd':9},
    '122Fano.pickle': {'CorrespondingK':'5', 'WallSecondChamber':None,                            'ExpectedFnd':7},
    '123Fano.pickle': {'CorrespondingK':'5', 'WallSecondChamber':None,                            'ExpectedFnd':7},
    '143Fano.pickle': {'CorrespondingK':'2', 'WallSecondChamber': [2, 1, 2, 3, 3, 3, 2, 2, 1, 0], 'ExpectedFnd':8},
    '144Fano.pickle': {'CorrespondingK':'3', 'WallSecondChamber':None,                            'ExpectedFnd':8},
    '158Fano.pickle': {'CorrespondingK':'3', 'WallSecondChamber':None,                            'ExpectedFnd':9},
    '159Fano.pickle': {'CorrespondingK':'5', 'WallSecondChamber':None,                            'ExpectedFnd':7},
    '171Fano.pickle': {'CorrespondingK':'3', 'WallSecondChamber':None,                            'ExpectedFnd':8},
    '176Fano.pickle': {'CorrespondingK':'5', 'WallSecondChamber':None,                            'ExpectedFnd':7},
}


try:
    import pickle
except:
    from pickle5 import pickle


def getVectorFromDescription(itemDescription,curvesList,automorphismList):
    StartingVector=vector(curvesList[itemDescription[0]])
    AutomorphismList=[matrix(automorphismList[itemDescription[t]]) for t in range(1,len(itemDescription)) ]
    for mAut in AutomorphismList:
        StartingVector = StartingVector * mAut
    return list(StartingVector)

GramL10=matrix([ [ -2, 0, 0, 1, 0, 0, 0, 0, 0, 0 ], [ 0, -2, 1, 0, 0, 0, 0, 0, 0, 0 ], [ 0, 1, -2, 1, 0, 0, 0, 0, 0, 0 ], 
  [ 1, 0, 1, -2, 1, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 1, -2, 1, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 1, -2, 1, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 1, -2, 1, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 1, -2, 1, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 1, -2, 1 ], [ 0, 0, 0, 0, 0, 0, 0, 0, 1, -2 ] ])
    

    
def checkIntersection(thisData):
    for fp in thisData['ListDelta']:
        for rat in thisData['ListDelta'][fp]:
            VectCurve=getVectorFromDescription(rat,thisData['RatsUsed'],thisData['AutUsed'])
            if (vector(fp) * GramL10 * vector(VectCurve) != 0):
                return False
    return True

def checkExpectedCnd(thisData, expectedCnd):
    """
    For each polarization fp, the stored list thisData['ListDelta'][fp] contains curves
    with intersection 0 with fp. For an upper bound ExpectedFnd, we require at least
    10 - ExpectedFnd such curves.
    """
    required_min = 10 - expectedCnd
    for fp in thisData['ListDelta']:
        if len(thisData['ListDelta'][fp]) < required_min:
            return False
    return True

def myintProduct(x, y):
    return vector(x)*GramL10*vector(y)

def Reflection(V,W):
    C=-2*myintProduct(V,W)/myintProduct(W,W)
    return [V[0]+C*W[0],V[1]+C*W[1],V[2]+C*W[2],V[3]+C*W[3],V[4]+C*W[4],V[5]+C*W[5],V[6]+C*W[6],V[7]+C*W[7],V[8]+C*W[8],V[9]+C*W[9]]



def checkFanoPolarization(thisData, CorrespondingK, WallSecondChamber):
    SetCurrentFp=thisData['ListDelta'].keys()
    with open(FiniteAut_folder_path + '/K' + CorrespondingK + '_Fano.pickle', 'rb') as f:
        ListFpImportedFromKondo = pickle.load(f)
    
    SetFpImportedFromKondo=set([tuple(v) for v in ListFpImportedFromKondo])
    
    if WallSecondChamber!=None:
        for ifp in ListFpImportedFromKondo:
            SetFpImportedFromKondo.add(tuple(Reflection(ifp,WallSecondChamber)))
        
    return SetCurrentFp==SetFpImportedFromKondo
    


for key in DataToCheck:
    print(f"Checking TTB {key}")
    with open(TTB_Fano_folder_path + '/' + key, 'rb') as f:
        thisData = pickle.load(f)
        
    if checkIntersection(thisData)==False:
        print(f"\tError: Intersection non zero between curves and Fano polarization")
        continue

    if checkExpectedCnd(thisData, DataToCheck[key]['ExpectedFnd'] )==False:
        print(f"\tError: Upper bound not realized")
        continue
        
    if checkFanoPolarization(thisData, DataToCheck[key]['CorrespondingK'], DataToCheck[key]['WallSecondChamber'])==False:
        print(f"\tError: Set of Fano polarization not corresponding")
        continue
            
            
    print(f"\tCondluded successfully without errors.\n")
              
    
    
