import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from collections import defaultdict
import argparse

def parse_profiling_file(filename):
    """Parse the profiling results file and extract data."""
    data = defaultdict(list)
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            # Match lines with braid length and degree data
            match = re.search(r'Braid length: (\d+), Degree: (\d+), Time taken: ([\d.]+) seconds', line)
            if match:
                braid_length = int(match.group(1))
                degree = int(match.group(2))
                time = float(match.group(3))
                data[braid_length].append((degree, time))
    
    # Sort data by degree for each braid length
    for braid_length in data:
        data[braid_length].sort(key=lambda x: x[0])
    
    return dict(data)

def exponential_func(x, a, b, c):
    """Exponential function: f(x) = a * exp(b * x) + c"""
    return a * np.exp(b * x) + c

def fit_exponential_curve(degrees, times):
    """Fit exponential curve to the data."""
    try:
        # Initial guess for parameters
        p0 = [1, 0.1, 0]
        
        # Fit the curve
        popt, pcov = curve_fit(exponential_func, degrees, times, p0=p0, maxfev=5000)
        
        # Calculate R-squared
        y_pred = exponential_func(np.array(degrees), *popt)
        ss_res = np.sum((np.array(times) - y_pred) ** 2)
        ss_tot = np.sum((np.array(times) - np.mean(times)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        
        return popt, r_squared
    except Exception as e:
        print(f"Failed to fit exponential curve: {e}")
        return None, 0

def plot_data_with_fits(data, output_file='profiling_analysis.png'):
    """Plot the data and fitted curves for all braid lengths."""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Profiling Results: Time vs Degree with Exponential Fits', fontsize=16)
    
    braid_lengths = sorted(data.keys())
    colors = ['blue', 'red', 'green', 'orange']
    
    for i, braid_length in enumerate(braid_lengths):
        if i >= 4:  # Only plot first 4 braid lengths
            break
            
        row = i // 2
        col = i % 2
        ax = axes[row, col]
        
        degrees, times = zip(*data[braid_length])
        degrees = np.array(degrees)
        times = np.array(times)
        
        # Plot original data
        ax.scatter(degrees, times, color=colors[i], alpha=0.7, s=50, 
                  label=f'Braid Length {braid_length} (Data)')
        
        # Fit exponential curve
        params, r_squared = fit_exponential_curve(degrees, times)
        
        if params is not None:
            # Generate smooth curve for plotting
            x_smooth = np.linspace(min(degrees), max(degrees), 100)
            y_smooth = exponential_func(x_smooth, *params)
            
            ax.plot(x_smooth, y_smooth, color=colors[i], linewidth=2, 
                   label=f'Exponential Fit (R²={r_squared:.3f})')
            
            # Add equation to plot
            a, b, c = params
            equation = f'f(x) = {a:.3f}*exp({b:.3f}*x) + {c:.3f}'
            ax.text(0.05, 0.95, equation, transform=ax.transAxes, 
                   fontsize=8, verticalalignment='top', 
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.set_xlabel('Degree')
        ax.set_ylabel('Time (seconds)')
        ax.set_title(f'Braid Length {braid_length}')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_yscale('log')  # Log scale for better visualization of exponential growth
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Plot saved as {output_file}")

def forecast_runtime(braid_length_data, target_degrees):
    """Forecast runtime for target degrees based on fitted exponential curves."""
    print("\n=== Runtime Forecasting ===")
    
    for braid_length, data_points in braid_length_data.items():
        degrees, times = zip(*data_points)
        degrees = np.array(degrees)
        times = np.array(times)
        
        params, r_squared = fit_exponential_curve(degrees, times)
        
        if params is not None:
            print(f"\nBraid Length {braid_length}:")
            print(f"  Fit quality (R²): {r_squared:.4f}")
            a, b, c = params
            print(f"  Equation: f(x) = {a:.3f}*exp({b:.3f}*x) + {c:.3f}")
            
            print("  Forecasted runtimes:")
            for target_degree in target_degrees:
                if target_degree > max(degrees):
                    predicted_time = exponential_func(target_degree, *params)
                    if predicted_time < 3600:  # Less than 1 hour
                        print(f"    Degree {target_degree}: {predicted_time:.2f} seconds")
                    elif predicted_time < 86400:  # Less than 1 day
                        print(f"    Degree {target_degree}: {predicted_time/3600:.2f} hours")
                    else:
                        print(f"    Degree {target_degree}: {predicted_time/86400:.2f} days")
                else:
                    print(f"    Degree {target_degree}: Already measured")

def main():
    parser = argparse.ArgumentParser(description='Analyze profiling results and fit exponential curves')
    parser.add_argument('--input', '-i', default='profiling_results.txt', 
                       help='Input profiling file (default: profiling_results.txt)')
    parser.add_argument('--output', '-o', default='profiling_analysis.png',
                       help='Output plot file (default: profiling_analysis.png)')
    parser.add_argument('--forecast', '-f', nargs='+', type=int, default=[26, 27, 28, 30],
                       help='Degrees to forecast (default: 26 27 28 30)')
    
    args = parser.parse_args()
    
    # Parse data
    print(f"Parsing data from {args.input}...")
    data = parse_profiling_file(args.input)
    
    if not data:
        print("No data found in the input file!")
        return
    
    print(f"Found data for braid lengths: {sorted(data.keys())}")
    
    # Plot data and fits
    plot_data_with_fits(data, args.output)
    
    # Forecast runtimes
    forecast_runtime(data, args.forecast)

if __name__ == "__main__":
    main()