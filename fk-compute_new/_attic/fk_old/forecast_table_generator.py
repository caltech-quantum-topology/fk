#!/usr/bin/env python3
"""
Forecast Table Generator
Generate comprehensive forecast tables showing runtime and clock cycles for degrees 5-50.
"""

import argparse
from profiling_utils import (
    parse_profiling_data, 
    print_forecast_table, 
    save_forecast_table_csv,
    generate_forecast_table
)

def main():
    parser = argparse.ArgumentParser(description='Generate forecast tables with runtime and clock cycles')
    parser.add_argument('--input', '-i', default='profiling_results.txt',
                       help='Input profiling file (default: profiling_results.txt)')
    parser.add_argument('--cpu-freq', '-f', type=float, default=3.2,
                       help='CPU frequency in GHz (default: 3.2 for M1 2020)')
    parser.add_argument('--degree-min', type=int, default=5,
                       help='Minimum degree to forecast (default: 5)')
    parser.add_argument('--degree-max', type=int, default=50,
                       help='Maximum degree to forecast (default: 50)')
    parser.add_argument('--braid-length', '-b', type=int, 
                       help='Generate table for specific braid length only')
    parser.add_argument('--save-csv', action='store_true',
                       help='Save tables as CSV files')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Only show summary, not full tables')
    
    args = parser.parse_args()
    
    # Parse data
    print(f"Loading data from {args.input}...")
    data = parse_profiling_data(args.input)
    
    if not data:
        print("No data found in the input file!")
        return
    
    degree_range = range(args.degree_min, args.degree_max + 1)
    braid_lengths = [args.braid_length] if args.braid_length else sorted(data.keys())
    
    print(f"Found data for braid lengths: {sorted(data.keys())}")
    print(f"Generating forecasts for degrees {args.degree_min}-{args.degree_max}")
    print(f"CPU frequency: {args.cpu_freq} GHz")
    print()
    
    # Process each braid length
    for braid_length in braid_lengths:
        if braid_length not in data:
            print(f"Warning: No data found for braid length {braid_length}")
            continue
            
        degrees, times = zip(*data[braid_length])
        
        if not args.quiet:
            # Print full table
            print_forecast_table(degrees, times, braid_length, degree_range, args.cpu_freq)
        
        # Save CSV if requested
        if args.save_csv:
            csv_filename = f'forecast_braid_{braid_length}_deg_{args.degree_min}-{args.degree_max}.csv'
            save_forecast_table_csv(degrees, times, braid_length, csv_filename, degree_range, args.cpu_freq)
        
        # Show summary for key degrees
        if args.quiet:
            results, r2 = generate_forecast_table(degrees, times, [30, 35, 40, 45, 50], args.cpu_freq)
            print(f"\nBraid Length {braid_length} (R²={r2:.4f}):")
            print("  Key forecasts:")
            for result in results:
                if result['status'] == 'predicted':
                    print(f"    Degree {result['degree']}: {result['formatted_time']} ({result['formatted_cycles']} cycles)")
        
        print("-" * 80)
    
    print("\nForecast generation complete!")
    
    # Show extreme forecasts as a warning
    print("\n⚠️  EXTREME FORECAST WARNING:")
    print("For very high degrees (40+), runtimes become astronomical.")
    print("Example for degree 50:")
    
    # Use the first braid length for extreme example
    first_braid = sorted(data.keys())[0]
    degrees, times = zip(*data[first_braid])
    results, _ = generate_forecast_table(degrees, times, [50], args.cpu_freq)
    
    if results:
        result = results[0]
        cycles = result['cycles']
        years = result['time_seconds'] / (365.25 * 24 * 3600)
        print(f"  Degree 50: {result['formatted_time']} (~{years:.1e} years)")
        print(f"  Clock cycles: {result['formatted_cycles']}")
        print("  This is likely far beyond practical computation limits!")

if __name__ == "__main__":
    main()