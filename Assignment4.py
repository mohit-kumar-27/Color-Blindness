import numpy as np
import functools as func
import operator as operate
import time
"""
Please refer to lecture slides.
Please refer to README file.
All the functions that you define must be able to handle corner cases/exceptions
"""

"""
Starting and ending locations (indices) of red and green exons in the reference sequence - Begins

1. Red Exon Locations
"""
RedExonPos = np.array([
    [149249757, 149249868], # R1
    [149256127, 149256423], # R2
    [149258412, 149258580], # R3
    [149260048, 149260213], # R4
    [149261768, 149262007], # R5
    [149264290, 149264400]  # R6
    ])
"""
2. Green Exon Locations
"""
GreenExonPos = np.array([
    [149288166, 149288277], # G1
    [149293258, 149293554], # G2
    [149295542, 149295710], # G3
    [149297178, 149297343], # G4
    [149298898, 149299137], # G5
    [149301420, 149301530]  # G6
    ])
"""
Starting and ending locations (indices) of red and green exons in the reference sequence - Ends
"""    

def loadLastCol(filename):
    """
    Input: Path of the file corresponding the last column (BWT).
    Output: The last column (BWT) in string format.
    """
    # function body - Begins
    LastCol = ''.join(np.loadtxt(filename, dtype=str))
    # function body - Ends
    return LastCol #string data type

def loadRefSeq(filename):
    """
    Input: Path of the file containing the reference sequence.
    Output: The reference sequence in string format.
    """
    # function body - Begins
    RefSeq = ''.join(np.loadtxt(filename, dtype=str)[1:])
    # function body - Ends
    return RefSeq # string data type

def loadReads(filename):
    """
    Input: Path of the file containing all the reads.
    Output: A list containing all the reads.
    """
    # function body - Begins
    Reads = np.loadtxt(filename, dtype=str)
    # function body - Ends
    return Reads # list of strings

def loadMapToRefSeq(filename):
    """
    Input: Path of the file containing mapping of the first column to the reference sequence.
    Output: numpy integer array containing the mapping.
    """
    # function body - Begins
    MapToRefSeq = np.loadtxt(filename, dtype=int)
    # function body - Ends
    return MapToRefSeq # numpy integer array

# This function finds the indices in the reference sequence where the read exactly matches the reference sequence
def get_exact_match_idcs(read):
    read.replace("N", "A") # replace N with A
    batch = FirstCol
    try:
        batch_start, batch_end = batch.index(read[-1]), batch.rindex(read[-1]) # finds the start index and end index of the batch of circular permutations starting with the last character of the read.
    except:
        return []
    for i in range(2,len(read)+1):
        last = LastCol[batch_start : batch_end+1] # batch containing the rightmost part of the string as the leftmost part in the BWT table.
        try:
            sub_batch_start, sub_batch_end = Rank[batch_start + last.index(read[-i])], Rank[batch_start + last.rindex(read[-i])] # finding the next left character from the right in the string in the last column of the batch obtained from above and finding its rank in the first column
        except:
            return []
        batch_start, batch_end = sub_batch_start, sub_batch_end 
        batch = FirstCol[batch_start: batch_end+1]
    return [Map[i] for i in np.arange(batch_start,batch_end+1)]

# This function converts a string to np array of its characters
def string_to_nparray(string):    
    return np.fromiter(string, (str,1))

def MatchReadToLoc(read):
    """
    Input: a read (string)
    Output: list of potential locations at which the read may match the reference sequence. 
    Refer to example in the README file.
    IMPORTANT: This function must also handle the following:
        1. cases where the read does not match with the reference sequence
        2. any other special case that you may encounter
    """
    # function body - Begins
    # allow max 2 mismatches
    l = len(read)

    # divide the read into 3 parts and apply get_exact_match_idcs
    first_part = np.array(get_exact_match_idcs(read[:l//3]))
    second_part = np.array(get_exact_match_idcs(read[l//3:2*(l//3)]))
    second_part = second_part - l//3        # calculate the initial index of the complete read
    third_part = np.array(get_exact_match_idcs(read[2*(l//3):]))
    third_part = third_part- 2*(l//3)       # calculate the initial index of the complete read

    # take union of indices obtained by the three parts
    positions = np.unique(np.concatenate((first_part,second_part,third_part)).astype('int')).tolist()

    for i in positions:
        if i+l > len(RefSeq) or np.sum((string_to_nparray(RefSeq[i: i+l]) == string_to_nparray(read))) < l-2: # check if mismatches are greater than 2, then discard
            positions.remove(i)
    # function body - Ends
    return positions # list of potential locations at which the read may match the reference sequence.

def WhichExon(positions):
    """
    Input: list of potential locations at which the read may match the reference sequence.
    Output: Update(increment) to the counts of the 12 exons
    IMPORTANT: This function must also handle the following:
        1. cases where the read does not match with the reference sequence
        2. cases where there are more than one matches (match at two exons)
        3. any other special case that you may encounter
    """
    r1,r2,r3,r4,r5,r6,g1,g2,g3,g4,g5,g6 = 0,0,0,0,0,0,0,0,0,0,0,0
    # function body - Begins
    if any([(RedExonPos[0][0] <= j <= RedExonPos[0][1]) for j in positions]) and any([(GreenExonPos[0][0] <= j <= GreenExonPos[0][1]) for j in positions]):
        r1 += 0.5
        g1 += 0.5
    elif any([(RedExonPos[0][0] <= j <= RedExonPos[0][1]) for j in positions]):
        r1 += 1
    elif any([(GreenExonPos[0][0] <= j <= GreenExonPos[0][1]) for j in positions]):
        g1 += 1
    
    elif any([(RedExonPos[1][0] <= j <= RedExonPos[1][1]) for j in positions]) and any([(GreenExonPos[1][0] <= j <= GreenExonPos[1][1]) for j in positions]):
        r2 += 0.5
        g2 += 0.5
    elif any([(RedExonPos[1][0] <= j <= RedExonPos[1][1]) for j in positions]):
        r2 += 1
    elif any([(GreenExonPos[1][0] <= j <= GreenExonPos[1][1]) for j in positions]):
        g2 += 1
    
    elif any([(RedExonPos[2][0] <= j <= RedExonPos[2][1]) for j in positions]) and any([(GreenExonPos[2][0] <= j <= GreenExonPos[2][1]) for j in positions]):
        r3 += 0.5
        g3 += 0.5
    elif any([(RedExonPos[2][0] <= j <= RedExonPos[2][1]) for j in positions]):
        r3 += 1
    elif any([(GreenExonPos[2][0] <= j <= GreenExonPos[2][1]) for j in positions]):
        g3 += 1

    elif any([(RedExonPos[3][0] <= j <= RedExonPos[3][1]) for j in positions]) and any([(GreenExonPos[3][0] <= j <= GreenExonPos[3][1]) for j in positions]):
        r4 += 0.5
        g4 += 0.5
    elif any([(RedExonPos[3][0] <= j <= RedExonPos[3][1]) for j in positions]):
        r4 += 1
    elif any([(GreenExonPos[3][0] <= j <= GreenExonPos[3][1]) for j in positions]):
        g4 += 1
    
    elif any([(RedExonPos[4][0] <= j <= RedExonPos[4][1]) for j in positions]) and any([(GreenExonPos[4][0] <= j <= GreenExonPos[4][1]) for j in positions]):
        r5 += 0.5
        g5 += 0.5
    elif any([(RedExonPos[4][0] <= j <= RedExonPos[4][1]) for j in positions]):
        r5 += 1
    elif any([(GreenExonPos[4][0] <= j <= GreenExonPos[4][1]) for j in positions]):
        g5 += 1
    
    elif any([(RedExonPos[5][0] <= j <= RedExonPos[5][1]) for j in positions]) and any([(GreenExonPos[5][0] <= j <= GreenExonPos[5][1]) for j in positions]):
        r6 += 0.5
        g6 += 0.5
    elif any([(RedExonPos[5][0] <= j <= RedExonPos[5][1]) for j in positions]):
        r6 += 1
    elif any([(GreenExonPos[5][0] <= j <= GreenExonPos[5][1]) for j in positions]):
        g6 += 1
    # function body - Ends    
    return np.array([r1,r2,r3,r4,r5,r6,g1,g2,g3,g4,g5,g6])

def ncr(n, r):
    r = min(r, n-r)
    num = func.reduce(operate.mul, range(n, n-r, -1), 1)
    den = func.reduce(operate.mul, range(1, r+1), 1)
    return num/den

def ComputeProb(ExonMatchCounts):
    """
    Input: The counts for each exon
    Output: Probabilities of each of the four configurations (a list of four real numbers)
    """
    # function body - Begins

    red_exon_cnt = [int(c) for c in ExonMatchCounts[:6]]
    green_exon_cnt = [int(c) for c in ExonMatchCounts[6:]]
    total_cnt = []
    
    #Truncating Matches to the integral values
    for i in range(6):
        cnt = red_exon_cnt[i] + green_exon_cnt[i]
        total_cnt.append(cnt)

    config1red = [1/3, 1/3, 1/3, 1/3]  
    config1green = [2/3,2/3,2/3,2/3]

    config2red = [1/2, 1/2, 0, 0]      
    config2green = [1/2, 1/2, 1.0, 1.0]

    config3red = [1/4, 1/4, 1/2,1/2]   
    config3green = [3/4, 3/4, 1/2, 1/2]

    config4red = [1/4,1/4,1/4, 1/2]    
    config4green = [3/4, 3/4, 3/4,1/2]

    P0 = 1
    for i in range(len(config1red)):
        P0 *= ncr(total_cnt[i+1], red_exon_cnt[i+1])*pow(config1red[i],red_exon_cnt[i+1])*pow(config1green[i],green_exon_cnt[i+1])
    P1 = 1
    for i in range(len(config2red)):
        P1 *= ncr(total_cnt[i+1], red_exon_cnt[i+1])*pow(config2red[i],red_exon_cnt[i+1])*pow(config2green[i],green_exon_cnt[i+1])  
    P2 = 1
    for i in range(len(config3red)):
        P2 *= ncr(total_cnt[i+1], red_exon_cnt[i+1])*pow(config3red[i],red_exon_cnt[i+1])*pow(config3green[i],green_exon_cnt[i+1])
    P3 = 1
    for i in range(len(config4red)):
        P3 *= ncr(total_cnt[i+1], red_exon_cnt[i+1])*pow (config4red[i],red_exon_cnt[i+1])*pow(config4green[i],green_exon_cnt[i+1])
    # function body - ends
    return [P0, P1, P2, P3]

def BestMatch(ListProb):
    """
    Input: Probabilities of each of the four configurations (a list of four real numbers)
    Output: Most likely configuration (an integer). Refer to lecture slides
    """
    # function body - Begins
    MostLikelyConfiguration = np.asarray(ListProb).argmax()
    # function body - ends
    return MostLikelyConfiguration # it holds 0, 1, 2, or 3

if __name__ == "__main__":
    start = time.time()
    # load all the data files
    LastCol = loadLastCol("../data/chrX_last_col.txt") # loads the last column
    RefSeq = loadRefSeq("../data/chrX.fa") # loads the reference sequence
    Reads = loadReads("../data/reads") # loads the reads
    Map = loadMapToRefSeq("../data/chrX_map.txt") # loads the mapping to the reference sequence
    a,c,g,t = LastCol.count('A'), LastCol.count('C'), LastCol.count('G'), LastCol.count('T')
    FirstCol = 'A'*a + 'C'*c + 'G'*g + 'T'*t +'$' # create the first column as a string

    # find rank in the first column for each element in the last column and store it as a list(Rank) 
    A, C, G, T = 0, a, a+c, a+c+g
    Rank = []
    for i in LastCol:
        if i == 'A':
            Rank.append(A)
            A += 1
        elif i == 'C':
            Rank.append(C)
            C += 1
        elif i == 'G':
            Rank.append(G)
            G += 1
        elif i == 'T':
            Rank.append(T)
            T += 1
        elif i == '$':
            Rank.append(len(FirstCol)-1)


    # run the functions
    ExonMatchCounts = np.zeros(12) # initialize the counts for exons
    for read in Reads[2936000:2948000]: # update the counts for exons
        positions = MatchReadToLoc(read) # get the list of potential match locations
        ExonMatchCounts += WhichExon(positions) # update the counts of exons, if applicable
    print("ExonMatchCounts : ", ExonMatchCounts)
    ListProb = ComputeProb(ExonMatchCounts) # compute probabilities of each of the four configurations
    print("ListProb : ", ListProb)
    MostLikely = BestMatch(ListProb) # find the most likely configuration
    print("Configuration %d is the best match"%MostLikely)
    end = time.time()
    print(f"Total time elapsed = {end-start:0.2f} s.\n")