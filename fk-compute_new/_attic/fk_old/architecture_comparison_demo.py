#!/usr/bin/env python3
"""
Architecture Comparison Demo
Demonstrates the new architecture-independent forecasting capabilities.
"""

from profiling_utils import parse_profiling_data, print_forecast_table
import argparse

def main():
    parser = argparse.ArgumentParser(description='Compare performance across different CPU architectures')
    parser.add_argument('--input', '-i', default='profiling_results.txt',
                       help='Input profiling file (default: profiling_results.txt)')
    parser.add_argument('--degree-min', type=int, default=25,
                       help='Minimum degree to show (default: 25)')
    parser.add_argument('--degree-max', type=int, default=35,
                       help='Maximum degree to show (default: 35)')
    parser.add_argument('--braid-length', '-b', type=int, default=8,
                       help='Braid length to analyze (default: 8)')
    
    args = parser.parse_args()
    
    # Load profiling data
    print(f"🔍 Loading profiling data from {args.input}...")
    data = parse_profiling_data(args.input)
    
    if not data:
        print("❌ No data found!")
        return
    
    if args.braid_length not in data:
        print(f"❌ No data found for braid length {args.braid_length}")
        print(f"Available braid lengths: {sorted(data.keys())}")
        return
    
    degrees, times = zip(*data[args.braid_length])
    degree_range = range(args.degree_min, args.degree_max + 1)
    
    print(f"\n🚀 ARCHITECTURE PERFORMANCE COMPARISON")
    print(f"{'='*70}")
    print(f"📊 Based on M1 2020 MacBook Air @ 3.2 GHz ground truth data")
    print(f"🎯 Braid Length: {args.braid_length} | Degrees: {args.degree_min}-{args.degree_max}")
    print(f"📈 Exponential extrapolation for unmeasured degrees")
    print()
    
    # Define CPU architectures to compare
    cpu_architectures = [
        ("🍎 Apple M1 2020", 3.2, "Current reference system"),
        ("🍎 Apple M1 Pro", 3.2, "Higher performance variant"),
        ("🍎 Apple M2", 3.49, "Next generation Apple Silicon"),
        ("🍎 Apple M3", 4.05, "Latest Apple Silicon"),
        ("⚡ Intel i7-12700K", 3.6, "High-end Intel desktop"),
        ("⚡ Intel i9-13900K", 5.8, "Top-tier Intel desktop"),
        ("🔥 AMD Ryzen 7 7700X", 4.5, "High-end AMD desktop"),  
        ("🔥 AMD Ryzen 9 7950X", 4.5, "Flagship AMD desktop"),
        ("🚀 Future CPU (8GHz)", 8.0, "Hypothetical future architecture"),
        ("🔬 Research CPU (12GHz)", 12.0, "Theoretical research system")
    ]
    
    # Show comparison for each architecture
    for cpu_name, cpu_freq, description in cpu_architectures:
        print(f"\n{cpu_name} @ {cpu_freq} GHz")
        print(f"   {description}")
        print_forecast_table(degrees, times, args.braid_length, degree_range, cpu_freq, 3.2)
        print()
    
    # Summary table showing key degree performance
    print(f"\n📋 PERFORMANCE SUMMARY TABLE")
    print(f"{'='*120}")
    print(f"{'Architecture':<25} {'Freq':<8} {'Degree 25':<12} {'Degree 30':<12} {'Degree 35':<12} {'Speedup vs M1':<15}")
    print(f"{'-'*120}")
    
    # Import generate_forecast_table for summary
    from profiling_utils import generate_forecast_table
    
    m1_results, _ = generate_forecast_table(degrees, times, [25, 30, 35], 3.2, 3.2)
    m1_times = {r['degree']: r['predicted_time_seconds'] for r in m1_results}
    
    for cpu_name, cpu_freq, description in cpu_architectures:
        results, _ = generate_forecast_table(degrees, times, [25, 30, 35], cpu_freq, 3.2)
        
        deg25_time = next(r['predicted_formatted_time'] for r in results if r['degree'] == 25)
        deg30_time = next(r['predicted_formatted_time'] for r in results if r['degree'] == 30)
        deg35_time = next(r['predicted_formatted_time'] for r in results if r['degree'] == 35)
        
        # Calculate speedup vs M1 for degree 30
        result_30 = next(r for r in results if r['degree'] == 30)
        speedup = m1_times[30] / result_30['predicted_time_seconds']
        
        clean_name = cpu_name.replace('🍎 ', '').replace('⚡ ', '').replace('🔥 ', '').replace('🚀 ', '').replace('🔬 ', '')
        
        print(f"{clean_name:<25} {cpu_freq:<8.1f} {deg25_time:<12} {deg30_time:<12} {deg35_time:<12} {speedup:<15.2f}x")
    
    print(f"\n💡 KEY INSIGHTS:")
    print(f"   • Ground truth cycles remain constant across architectures")
    print(f"   • Higher frequency CPUs provide linear speedup benefits")
    print(f"   • M3 @ 4.05GHz: ~27% faster than M1 @ 3.2GHz")
    print(f"   • Intel i9-13900K @ 5.8GHz: ~81% faster than M1")
    print(f"   • Future 8GHz CPU: ~2.5x faster than current M1")
    print(f"\n⚠️  Note: Real-world performance may vary due to:")
    print(f"   • Architectural differences (IPC, cache, memory bandwidth)")
    print(f"   • Thermal throttling under sustained load")
    print(f"   • Compiler optimizations for specific architectures")

if __name__ == "__main__":
    main()