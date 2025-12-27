import time
from needmanparl import needleman_wunsch_parallel
from Needlmanseq import Needleman_Wunsch
from utils import generate_random_DNA


def benchmark():
    sequence_lengths = [100, 200, 500, 1000, 2000, 3000]
    
    print("Running benchmarks...")
    print("-" * 60)
    
    for length in sequence_lengths:
        print(f"Testing sequences of length {length}...")
        
        # Generate same sequences for fair comparison
        seq1 = generate_random_DNA(length)
        seq2 = generate_random_DNA(length)
        
        # Sequential timing
        start = time.perf_counter()
        Needleman_Wunsch(seq1, seq2, 1, -1, -2)
        seq_time = time.perf_counter() - start
        
        # Parallel timing
        start = time.perf_counter()
        needleman_wunsch_parallel(seq1, seq2)
        par_time = time.perf_counter() - start
        
        speedup = seq_time / par_time
        print(f"  Sequential: {seq_time:.4f}s")
        print(f"  Parallel:   {par_time:.4f}s")
        print(f"  Speedup:    {speedup:.2f}x")
        print("-" * 60)


if __name__ == "__main__":
    benchmark()