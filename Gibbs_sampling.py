import random

''' This script implements the Gibbs sampling algorithm for motif finding'''
def findMotifs(dna, k):
    """ This function Randomly chooses starting positions
     and form the set of k-mers associated with these starting positions."""
    motifs = []
    for seq in dna:
        start = random.randrange(len(seq)-k+1)
        subseq = seq[start:start + k]  # The chosen k-mer.
        motifs.append(subseq)
    return motifs

def counting(motifs):
    """ This function generates a position specific counting matrix. """
    k = len(motifs[0])  # The given k-mer length.
    nucleotides = ['A', 'T', 'C', 'G']
    counts = []
    for i in range(k):
        currCounts = {nuc: 0 for nuc in nucleotides}  # Initializes the matrix.
        for motif in motifs:
            currCounts[motif[i]] += 1  # Adds the occurrence of the nucleotide at the current position.
        counts.append([currCounts[nuc] for nuc in nucleotides])
    return counts


def Score(motifs):
    """ This function counts the number of mismatches between nucleotides
     in the list "motifs" and returns the mismatch count as the score
      - the higher the value the worse the score"""
    counts = counting(motifs)
    score = 0
    for count in counts:
        score += sum(count) - max(count)  # counts all the comparisons and removes the match.
    return score


def findProfile(motifs):
    """ This function generates a 4xk matrix containing the probabilities that
    the nth base in length k motif will be a specific base """
    k = float(len(motifs[0]))  # The length of the motif being searched for.
    counts = counting(motifs)  # The counting matrix.
    # Calculates the probabilty by dividing each cell in the matrix in k (adds pseudo 1 in order to avoid zeros).
    profile = [[(element+1)/k for element in count] for count in counts]
    return profile


def profileProb(kmer, profile):
    """ This function computes the probability of all k-mers in the removed sequence """
    prob = 1
    indices = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    for i in range(len(kmer)):
        currNucIndex = indices[kmer[i]]  # Extracts the current base.
        prob *= profile[i][currNucIndex]
    return prob


def randomPose(probs):
    """ This function selects the start positions of the motifs
    according to the computed probabilities """
    total = sum(probs)
    rand = random.random()  # A random number between 0 to 1.
    partialSum = 0.0
    for i in range(len(probs)):
        partialSum += probs[i]
        if partialSum/total >= rand:  # The highest probabilities.
            return i
    return -1

def gibbsSampler(dna, k, N):
    """ This function implements the Gibbs sampling algorithm -
    A motif finding algorithm that finds one common motif in a list of "dna" reads and returns a list bestMotifs that
    contains the closes motif match from each string in "dna "
    The function gets 3 arguments:
    "dna" -- A list of DNA reads that are the same length.
    "k" -- An integer indicating the motif length being searched for.
    "N" -- The number of iterations for the algorithm. """
    t = len(dna)
    motifs = findMotifs(dna, k)
    bestMotifs = motifs
    bestScore = Score(bestMotifs)

    for n in range(N):
        i = random.randrange(t-1)  # Randomly choose one of the given sequences.
        profile = findProfile(motifs[:i] + motifs[i+1:])  # Create a profile from the other sequences.
        # For each position  in the removed sequence calculate the probability that the k-mer at
        #  the current position was generated by the profile.
        probability = [profileProb(dna[i][s:s + k], profile) for s in range(len(dna[i]) - k + 1)]
        # Choose a new starting position for the removed sequence at random, based on the probability.
        pose = randomPose(probability)
        motifs[i] = dna[i][pose:pose + k]
        score = Score(motifs)
        if score < bestScore:  # The lower the score, the better the results.
            bestMotifs = motifs
            bestScore = score
    return bestMotifs, bestScore

def repeatGibbsSampler(dna, k, N, repeats):
    """ This function repeats the whole process in order to reach convergence. """
    bestMotifs = findMotifs(dna, k)
    bestScore = Score(bestMotifs)
    for i in range(repeats):
        (motifs, score) = gibbsSampler(dna, k, N)
        if score < bestScore:
            bestMotifs = motifs
            bestScore = score
    return bestMotifs

def main():
    pathToFile = raw_input("Enter the path for the file where the reads are: ")
    k = raw_input("Enter the length of the motif you want to search for: ")
    N = raw_input("Enter how many iterations do you want: ")
    R = raw_input("Enter how many times you want the whole process to be performed: ")
    dna = file(pathToFile).read()

    # Example for parameters:
    # pathToFile = "DNA.txt"
    # k = 8
    # N = 100
    # R = 20

    dna = file(pathToFile).read()
    reads = dna.split("\n")
    bestMotifs = repeatGibbsSampler(reads, int(k), int(N), int(R))
    print("\nThe best motifs were found by the algorithm:")
    print '\n'.join(str(p) for p in bestMotifs)


main()

# ------------------------------ Output ------------------------------
# Enter the path for the file where the reads are: DNA.txt
# Enter the length of the motif you want to search for: 8
# Enter how many iterations do you want: 100
# Enter how many times you want the whole process to be performed: 20
#
# The best motifs were found by the algorithm:
# CCCTCTCG
# GCGAGGTA
# AAGAAGTA
# TTCAGGTG
# TCCACGTG