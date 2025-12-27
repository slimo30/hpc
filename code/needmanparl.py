import numpy as np
import random
from multiprocessing import Pool, cpu_count
import time
from utils import generate_random_DNA


def compute_cell(args):
    i, j, seq1, seq2, F_prev_diag, F_prev_row, F_prev_col, match, mismatch, gap = args
    if seq1[i-1] == seq2[j-1]:
        score_diag = match
    else:
        score_diag = mismatch
    diag = F_prev_diag + score_diag
    up = F_prev_row + gap
    left = F_prev_col + gap
    return (i, j, max(diag, up, left))

def needleman_wunsch_parallel(seq1, seq2, match=1, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    F = np.zeros((n+1, m+1), dtype=int)
    
    F[1:,0] = np.arange(1, n+1) * gap
    F[0,1:] = np.arange(1, m+1) * gap
    
    pool = Pool(processes=cpu_count())
    
    start_matrix = time.perf_counter()
    for k in range(2, n+m+1):
        # Prepare tasks with the necessary values from F
        tasks = [(i, k-i, seq1, seq2, F[i-1, k-i-1], F[i-1, k-i], F[i, k-i-1], match, mismatch, gap) 
                 for i in range(1, n+1) if 1 <= k-i <= m]
        results = pool.map(compute_cell, tasks)
        for i, j, val in results:
            F[i,j] = val
    end_matrix = time.perf_counter()
    matrix_time = end_matrix - start_matrix
    
    # Traceback
    start_traceback = time.perf_counter()
    i, j = n, m
    align1, align2 = [], []
    while i>0 or j>0:
        if i>0 and j>0:
            score_diag = match if seq1[i-1]==seq2[j-1] else mismatch
            if F[i,j]==F[i-1,j-1]+score_diag:
                align1.append(seq1[i-1])
                align2.append(seq2[j-1])
                i-=1
                j-=1
                continue
        if i>0 and F[i,j]==F[i-1,j]+gap:
            align1.append(seq1[i-1])
            align2.append("-")
            i-=1
        else:
            align1.append("-")
            align2.append(seq2[j-1])
            j-=1
    end_traceback = time.perf_counter()
    traceback_time = end_traceback - start_traceback
    
    pool.close()
    pool.join()
    
    return "".join(reversed(align1)), "".join(reversed(align2)), F[n,m], matrix_time, traceback_time

if __name__ == "__main__":
    import time
    seq1 = generate_random_DNA(10000)
    seq2 = generate_random_DNA(10000)

    start_time = time.time()
    align1, align2, score, t_matrix, t_trace = needleman_wunsch_parallel(seq1, seq2)
    print("Alignment Score:", score)
    print("Time:", time.time() - start_time)
    print("mt",t_matrix,"mg",t_trace)