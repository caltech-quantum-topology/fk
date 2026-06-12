#!/usr/bin/env python3
"""
Batch Export Forecast Tool
Automatically generates and saves forecast figures and tables for all braid lengths,
organized by crossing number from CSV profiling data.
"""

import os
import csv
import matplotlib.pyplot as plt
import numpy as np
from profiling_utils_suffix import (
    parse_profiling_data_csv,
    generate_forecast_table,
    exponential_func,
    fit_exponential_to_data
)

def get_crossing_number_mapping(filename):
    """
    Extract the mapping between braid length and crossing number from CSV data.
    
    Args:
        filename (str): Path to CSV file
        
    Returns:
        dict: Mapping from braid_length to crossing_number
    """
    mapping = {}
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                braid_length = int(row['Braid_Length'])
                crossing_number = int(row['Knot_Index'])
                if braid_length not in mapping:
                    mapping[braid_length] = crossing_number
            except (ValueError, KeyError):
                continue
    return mapping

def create_forecast_plot(degrees, times, braid_length, crossing_number, degree_range, target_cpu_freq=3.2):
    """
    Create and return a forecast plot for given data.
    
    Args:
        degrees (list): List of degrees
        times (list): List of computation times
        braid_length (int): Braid length
        crossing_number (int): Crossing number
        degree_range (range): Range of degrees to forecast
        target_cpu_freq (float): Target CPU frequency
        
    Returns:
        tuple: (fig, ax, r_squared)
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    degrees = np.array(degrees)
    times = np.array(times)
    
    # Plot measured data
    ax.scatter(degrees, times, color='blue', alpha=0.7, s=60,
               label='Measured Data', zorder=3)
    
    # Fit exponential curve
    params, r_squared = fit_exponential_to_data(degrees, times)
    
    if params is not None:
        # Plot fitted curve
        x_range = np.linspace(min(degrees), max(degree_range), 200)
        y_fit = exponential_func(x_range, *params)
        
        ax.plot(x_range, y_fit, 'r-', linewidth=2,
                label=f'Exponential Fit (R²={r_squared:.4f})', zorder=2)
        
        # Plot forecast points
        forecast_degrees = [d for d in degree_range if d > max(degrees)]
        if forecast_degrees:
            forecast_times = [exponential_func(d, *params) for d in forecast_degrees]
            ax.scatter(forecast_degrees, forecast_times,
                      color='red', s=80, marker='^',
                      label='Forecasted Points', alpha=0.8, zorder=3)
        
        # Add equation
        a, b, c = params
        equation = f'f(x) = {a:.3f}·exp({b:.3f}·x) + {c:.3f}'
        ax.text(0.05, 0.95, equation, transform=ax.transAxes,
                fontsize=12, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Formatting
    ax.set_xlabel('Degree', fontsize=12)
    ax.set_ylabel('Computation Time (seconds)', fontsize=12)
    ax.set_title(f'Runtime Forecast - Crossing Number {crossing_number} (Braid Length {braid_length})\n'
                f'Target CPU: {target_cpu_freq} GHz', fontsize=14, fontweight='bold')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11)
    
    # Add statistics text box
    measured_count = len(degrees)
    predicted_count = len([d for d in degree_range if d > max(degrees)])
    total_count = len(degree_range)
    
    stats_text = (f'Crossing Number: {crossing_number}\n'
                 f'Braid Length: {braid_length}\n'
                 f'Measured Points: {measured_count}\n'
                 f'Predicted Points: {predicted_count}\n'
                 f'Total Degrees: {total_count}\n'
                 f'R² Fit Quality: {r_squared:.6f}')
    
    ax.text(0.02, 0.02, stats_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    return fig, ax, r_squared

def save_forecast_table_csv(degrees, times, braid_length, crossing_number, filename, 
                           degree_range, target_cpu_freq=3.2, reference_cpu_freq=3.2):
    """
    Save forecast table to CSV file with crossing number information.
    """
    results, r_squared = generate_forecast_table(degrees, times, degree_range, 
                                                target_cpu_freq, reference_cpu_freq)
    
    if not results:
        print(f"Failed to generate forecast table for crossing number {crossing_number}")
        return False
    
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header with metadata
        writer.writerow([f'# Forecast Table - Crossing Number {crossing_number} (Braid Length {braid_length})'])
        writer.writerow([f'# Reference CPU: {reference_cpu_freq} GHz (M1 2020)'])
        writer.writerow([f'# Target CPU: {target_cpu_freq} GHz'])
        writer.writerow([f'# Fit Quality (R²): {r_squared:.6f}'])
        writer.writerow([])
        
        # Column headers
        writer.writerow(['Degree', 'Ground_Truth_Time_Seconds', 'Ground_Truth_Formatted_Time',
                        'Ground_Truth_Cycles', 'Ground_Truth_Formatted_Cycles',
                        'Predicted_Time_Seconds', 'Predicted_Formatted_Time', 'Status'])
        
        # Data rows
        for result in results:
            writer.writerow([
                result['degree'],
                f"{result['ground_truth_time_seconds']:.6f}",
                result['ground_truth_formatted_time'],
                result['ground_truth_cycles'],
                result['ground_truth_formatted_cycles'],
                f"{result['predicted_time_seconds']:.6f}",
                result['predicted_formatted_time'],
                result['status']
            ])
    
    return True

def batch_export_forecasts(input_file='profiling_results2_suffix.txt', output_dir='forecast_exports',
                          degree_range=None, target_cpu_freqs=None):
    """
    Batch export all forecast figures and tables organized by crossing number.
    
    Args:
        input_file (str): Path to input CSV file
        output_dir (str): Output directory for exports
        degree_range (range): Range of degrees to forecast (default: 5-50)
        target_cpu_freqs (list): List of target CPU frequencies (default: [3.2])
    """
    if degree_range is None:
        degree_range = range(5, 51)
    
    if target_cpu_freqs is None:
        target_cpu_freqs = [3.2]
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse data
    print("Loading profiling data...")
    data = parse_profiling_data_csv(input_file)
    crossing_mapping = get_crossing_number_mapping(input_file)
    
    if not data:
        print("No data found in input file!")
        return
    
    print(f"Found data for {len(data)} braid lengths")
    print("Crossing number mapping:")
    for braid_length in sorted(crossing_mapping.keys()):
        print(f"  Braid Length {braid_length} → Crossing Number {crossing_mapping[braid_length]}")
    
    # Process each braid length
    export_summary = []
    
    for braid_length in sorted(data.keys()):
        crossing_number = crossing_mapping.get(braid_length, braid_length)
        degrees, times = zip(*data[braid_length])
        
        print(f"\nProcessing Crossing Number {crossing_number} (Braid Length {braid_length})...")
        print(f"  Measured degrees: {min(degrees)}-{max(degrees)} ({len(degrees)} points)")
        
        # Export for each target CPU frequency
        for target_cpu_freq in target_cpu_freqs:
            freq_str = f"{target_cpu_freq:.1f}GHz".replace('.', 'p')
            
            # Create and save plot
            fig, ax, r_squared = create_forecast_plot(degrees, times, braid_length, 
                                                    crossing_number, degree_range, target_cpu_freq)
            
            plot_filename = f"braid_length_{braid_length:02d}_forecast_{freq_str}.png"
            plot_path = os.path.join(output_dir, plot_filename)
            fig.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            
            # Save CSV table
            csv_filename = f"braid_length_{braid_length:02d}_table_{freq_str}.csv"
            csv_path = os.path.join(output_dir, csv_filename)
            success = save_forecast_table_csv(degrees, times, braid_length, crossing_number,
                                            csv_path, degree_range, target_cpu_freq, 3.2)
            
            if success:
                export_summary.append({
                    'crossing_number': crossing_number,
                    'braid_length': braid_length,
                    'target_cpu_freq': target_cpu_freq,
                    'r_squared': r_squared,
                    'measured_points': len(degrees),
                    'plot_file': plot_filename,
                    'csv_file': csv_filename
                })
                print(f"  ✓ Exported {freq_str}: R²={r_squared:.6f}")
            else:
                print(f"  ✗ Failed to export {freq_str}")
    
    # Create summary report
    summary_path = os.path.join(output_dir, 'export_summary.csv')
    with open(summary_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Braid_Length', 'Crossing_Number', 'Target_CPU_GHz', 'R_Squared',
                        'Measured_Points', 'Plot_File', 'CSV_File'])
        
        for item in export_summary:
            writer.writerow([
                item['braid_length'], item['crossing_number'], item['target_cpu_freq'],
                f"{item['r_squared']:.6f}", item['measured_points'],
                item['plot_file'], item['csv_file']
            ])
    
    print(f"\n🎉 Batch export complete!")
    print(f"📁 Output directory: {output_dir}")
    print(f"📊 Exported {len(export_summary)} forecast sets")
    print(f"📋 Summary saved to: {summary_path}")

def main():
    """Main function with customizable parameters."""
    print("🚀 Batch Forecast Export Tool")
    print("=" * 50)
    
    # Configuration
    input_file = 'profiling_results2_suffix.txt'
    output_dir = 'forecast_exports'
    degree_range = range(5, 51)  # Degrees 5-50
    target_cpu_freqs = [3.2]  # Default CPU frequency
    
    print(f"📁 Input file: {input_file}")
    print(f"📂 Output directory: {output_dir}")
    print(f"🎯 Degree range: {min(degree_range)}-{max(degree_range)}")
    print(f"💻 Target CPU frequencies: {target_cpu_freqs} GHz")
    print()
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"❌ Error: Input file '{input_file}' not found!")
        return
    
    # Run batch export
    batch_export_forecasts(input_file, output_dir, degree_range, target_cpu_freqs)

if __name__ == "__main__":
    main()