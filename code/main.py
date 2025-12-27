import time
import matplotlib.pyplot as plt
from needmanparl import needleman_wunsch_parallel
from Needlmanseq import Needleman_Wunsch
from utils import generate_random_DNA


def benchmark():
    sequence_lengths = [100, 200, 500, 1000, 2000, 3000]
    sequential_times = []
    parallel_times = []
    
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
        sequential_times.append(seq_time)
        
        # Parallel timing
        start = time.perf_counter()
        needleman_wunsch_parallel(seq1, seq2)
        par_time = time.perf_counter() - start
        parallel_times.append(par_time)
        
        speedup = seq_time / par_time
        print(f"  Sequential: {seq_time:.4f}s")
        print(f"  Parallel:   {par_time:.4f}s")
        print(f"  Speedup:    {speedup:.2f}x")
        print("-" * 60)
    
    return sequence_lengths, sequential_times, parallel_times


def plot_comparison(lengths, seq_times, par_times):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Execution Time Comparison
    ax1.plot(lengths, seq_times, 'o-', label='Sequential', linewidth=2, markersize=8, color='blue')
    ax1.plot(lengths, par_times, 's-', label='Parallel', linewidth=2, markersize=8, color='orange')
    ax1.set_xlabel('Sequence Length', fontsize=12)
    ax1.set_ylabel('Execution Time (seconds)', fontsize=12)
    ax1.set_title('Needleman-Wunsch: Sequential vs Parallel', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Speedup
    speedups = [s/p for s, p in zip(seq_times, par_times)]
    ax2.plot(lengths, speedups, 'o-', color='green', linewidth=2, markersize=8)
    ax2.axhline(y=1, color='red', linestyle='--', label='No speedup (1x)', linewidth=2)
    ax2.set_xlabel('Sequence Length', fontsize=12)
    ax2.set_ylabel('Speedup Factor', fontsize=12)
    ax2.set_title('Parallel Speedup Over Sequential', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('needleman_wunsch_comparison.png', dpi=300, bbox_inches='tight')
    print("\nGraph saved as 'needleman_wunsch_comparison.png'")
    plt.show()


if __name__ == "__main__":
    lengths, seq_times, par_times = benchmark()
    plot_comparison(lengths, seq_times, par_times)
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    speedups = [s/p for s, p in zip(seq_times, par_times)]
    print(f"Average Speedup: {sum(speedups)/len(speedups):.2f}x")
    print(f"Maximum Speedup: {max(speedups):.2f}x")
    print(f"Total Sequential Time: {sum(seq_times):.2f}s")
    print(f"Total Parallel Time: {sum(par_times):.2f}s")