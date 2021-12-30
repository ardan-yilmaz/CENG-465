import math
import sys



class Point:
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)


class ResiduePairs:
    def __init__(self, r1, r2, groupID):
        self.r1 = r1
        self.r2 = r2
        self.groupID = groupID


class Residue:
    def __init__(self, atom, resName, chainID, resID, order, x, y, z):
        self.atom = atom
        self.resName = resName
        self.chainID = chainID
        self.resID = resID
        self.order = order
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)


def computeDistance(p1, p2):
    d = (p2.x - p1.x) ** 2 + (p2.y - p1.y) ** 2 + (p2.z - p1.z) ** 2
    return math.sqrt(d)


# finds pairs of interacting residues between the two chains
def findInteractingPairs(residuesOfChainA, residuesOfChainB):
    # pairs is the list to hold the interacting residues
    pairs = []
    lenA = len(residuesOfChainA)
    lenB = len(residuesOfChainB)
    for i in range(0, lenA):
        p1 = Point(residuesOfChainA[i].x, residuesOfChainA[i].y, residuesOfChainA[i].z)
        for j in range(0, lenB):
            p2 = Point(residuesOfChainB[j].x, residuesOfChainB[j].y, residuesOfChainB[j].z)
            distance = computeDistance(p1, p2)
            if (distance < 8):
                pair = ResiduePairs(residuesOfChainA[i], residuesOfChainB[j], -1)
                pairs.append(pair)
            else:
                continue
    findHotSpots(pairs)


def findHotSpots(pairs):
    length = len(pairs)
    pairs[0].groupID = 1
    for i in range(1, length):
        flag = 0
        ifGrouped = 0
        pair1 = pairs[i]
        resA1 = pair1.r1
        resB1 = pair1.r2
        PA1 = Point(resA1.x, resA1.y, resA1.z)
        PB1 = Point(resB1.x, resB1.y, resB1.z)
        for j in range(i - 1, -1, -1):
            pair2 = pairs[j]
            resA2 = pair2.r1
            resB2 = pair2.r2
            PA2 = Point(resA2.x, resA2.y, resA2.z)
            PB2 = Point(resB2.x, resB2.y, resB2.z)
            if (ifHotSpot(PA1, PB1, PA2, PB2)):
                ifGrouped = 1
                if (pairs[j].groupID == 1):
                    pairs[i].groupID = 1
                    helper(pairs, i)
                    break
                else:
                    if (pairs[i].groupID == -1):
                        flag = 1
                        pairs[i].groupID = pairs[j].groupID
                        # helper(pairs, i)
                    else:
                        if (pairs[i].groupID > pairs[j].groupID):
                            flag = 1
                            pairs[i].groupID = pairs[j].groupID
                            # helper(pairs, i)
            else:
                continue


        if(ifGrouped == 0):
            maxID = 1
            for k in range(0, i):
                if(pairs[k].groupID > maxID):
                    maxID = pairs[k].groupID
            pairs[i].groupID = maxID+1

        elif(ifGrouped == 1 and flag == 1): helper(pairs, i)

    groupID_output = 1
    for i in range(0, length):
        if(pairs[i].groupID > groupID_output):
            groupID_output = pairs[i].groupID

    strH = "There are " + str(length) + " interacting pairs" + "\n"
    file1 = open("output.txt", "w")
    file1.writelines(strH)
    for i in range(0, len(pairs)):
        group = "Group " + str(pairs[i].groupID) + ": "
        res1 = str(pairs[i].r1.resName) + "(" + str(pairs[i].r1.order) + ")-"
        res2 = str(pairs[i].r2.resName) + "(" + str(pairs[i].r2.order) + ")"
        str1 = group + res1 + res2 + "\n"
        file1.writelines(str1)
    strF = "Number of groups = " + str(groupID_output)
    file1.writelines(strF)



def helper(pairs, index):
    pair1 = pairs[index]
    resA1 = pair1.r1
    resB1 = pair1.r2
    PA1 = Point(resA1.x, resA1.y, resA1.z)
    PB1 = Point(resB1.x, resB1.y, resB1.z)
    for i in range(0, index):
        if (pairs[i].groupID > pairs[index].groupID):
            pair2 = pairs[i]
            resA2 = pair2.r1
            resB2 = pair2.r2
            PA2 = Point(resA2.x, resA2.y, resA2.z)
            PB2 = Point(resB2.x, resB2.y, resB2.z)
            if (ifHotSpot(PA1, PB1, PA2, PB2)):
                pairs[i].groupID = pairs[index].groupID
                helper(pairs, i)





def ifHotSpot(pA1, pB1, pA2, pB2):
    d1 = computeDistance(pA1, pA2)
    d2 = computeDistance(pA1, pB2)
    d3 = computeDistance(pB1, pA2)
    d4 = computeDistance(pB1, pB2)
    if (min(d1, d2, d3, d4) < 8):
        return 1
    else:
        return 0


# main func
def driverFunc():
    file_name = sys.argv[1]
    inputt = './' + file_name

    file1 = open(inputt, 'r')
    line = file1.readline()

    residuesOfChainA = []
    residuesOfChainB = []
    cnt = 0
    while (line):
        if (line[0:4] == "ATOM"):
            tmp = line.split()
            if (tmp[3] != "GLY"):
                if (tmp[2] == "CB"):
                    if (tmp[4] == 'A'):
                        cnt += 1
                        res = Residue(tmp[2], tmp[3], tmp[4], tmp[5], cnt, tmp[6], tmp[7], tmp[8])
                        residuesOfChainA.append(res)
                    elif (tmp[4] == 'B'):
                        cnt += 1
                        res = Residue(tmp[2], tmp[3], tmp[4], tmp[5], cnt, tmp[6], tmp[7], tmp[8])
                        residuesOfChainB.append(res)
            else:
                if (tmp[2] == "CA"):
                    if (tmp[4] == 'A'):
                        cnt += 1
                        res = Residue(tmp[2], tmp[3], tmp[4], tmp[5], cnt, tmp[6], tmp[7], tmp[8])
                        residuesOfChainA.append(res)
                    elif (tmp[4] == 'B'):
                        cnt += 1
                        res = Residue(tmp[2], tmp[3], tmp[4], tmp[5], cnt, tmp[6], tmp[7], tmp[8])
                        residuesOfChainB.append(res)
        line = file1.readline()

    findInteractingPairs(residuesOfChainA, residuesOfChainB)





driverFunc()

