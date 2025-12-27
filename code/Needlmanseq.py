import time
import random
from utils import generate_random_DNA



def Needleman_Wunsch(ADN1, ADN2, Match_score, MisMatch_Penality, Gap_Penality):

    n = len(ADN1)
    m = len(ADN2)
    start = time.perf_counter()

    Matrix = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        Matrix[i][0] = i * Gap_Penality

    for j in range(1, m + 1):
        Matrix[0][j] = j * Gap_Penality

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if ADN1[i - 1] == ADN2[j - 1]:
                diag = Matrix[i - 1][j - 1] + Match_score
            else:
                diag = Matrix[i - 1][j - 1] + MisMatch_Penality

            up = Matrix[i - 1][j] + Gap_Penality
            left = Matrix[i][j - 1] + Gap_Penality

            Matrix[i][j] = max(diag, up, left)

    end = time.perf_counter()
    MatrixCreationTime = end - start
    start = time.perf_counter()

    NewADN1 = []
    NewADN2 = []

    i = n
    j = m

    while i > 0 or j > 0:

        if i > 0 and j > 0:
            if ADN1[i - 1] == ADN2[j - 1]:
                score_diag = Match_score
            else:
                score_diag = MisMatch_Penality

            if Matrix[i][j] == Matrix[i - 1][j - 1] + score_diag:
                NewADN1.append(ADN1[i - 1])
                NewADN2.append(ADN2[j - 1])
                i -= 1
                j -= 1
                continue

        if i > 0 and Matrix[i][j] == Matrix[i - 1][j] + Gap_Penality:
            NewADN1.append(ADN1[i - 1])
            NewADN2.append("-")
            i -= 1
        else:
            NewADN1.append("-")
            NewADN2.append(ADN2[j - 1])
            j -= 1

    NewADN1 = ''.join(reversed(NewADN1))
    NewADN2 = ''.join(reversed(NewADN2))

    end = time.perf_counter()
    TraceBackTime = end - start

    return NewADN1, NewADN2, Matrix[n][m], MatrixCreationTime, TraceBackTime

seq1 = generate_random_DNA(10000)
seq2 = generate_random_DNA(10000)
import time
start_time = time.time()
a1, a2, score, t_matrix, t_trace = Needleman_Wunsch(
    seq1,
    seq2,
    Match_score=1,
    MisMatch_Penality=-1,
    Gap_Penality=-2
)

print(a1)
print(a2)
print("Score:", score)
print("Matrix time:", t_matrix)
print("Traceback time:", t_trace)
print("Time:", time.time() - start_time)

import cProfile
import pstats
"""
def run():
    seq1 = generate_random_DNA(500)
    seq2 = generate_random_DNA(500)

    Needleman_Wunsch(
        seq1,
        seq2,
        Match_score=1,
        MisMatch_Penality=-1,
        Gap_Penality=-2
    )

cProfile.run("run()", "profiling_results")
stats = pstats.Stats("profiling_results")
stats.sort_stats("cumtime").print_stats(10)###
"""