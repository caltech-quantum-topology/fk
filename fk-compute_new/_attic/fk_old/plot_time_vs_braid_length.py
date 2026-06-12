import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

def exponential_model(x, a, b, c):
    """Exponential model: y = a * exp(b * x) + c"""
    return a * np.exp(b * x) + c

def load_timing_data(filename, braid_length):
    """Load timing data from CSV file"""
    data = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Skip header lines that start with #
    start_line = 0
    for i, line in enumerate(lines):
        if not line.strip().startswith('#') and 'Degree,Time_Seconds' in line:
            start_line = i + 1
            break
    
    for line in lines[start_line:]:
        if line.strip() and 'measured' in line:
            parts = line.strip().split(',')
            degree = int(parts[0])
            time_seconds = float(parts[1])
            data.append({'braid_length': braid_length, 'degree': degree, 'time_seconds': time_seconds})
    
    return data

def plot_time_vs_braid_length(all_data, selected_degrees=None):
    """Plot time vs braid length for each degree (one curve per degree)"""
    
    # Get all degrees and braid lengths
    all_degrees = sorted(set(d['degree'] for d in all_data))
    all_braid_lengths = sorted(set(d['braid_length'] for d in all_data))
    
    # Select degrees to plot (default: every 5th degree plus some key ones)
    if selected_degrees is None:
        selected_degrees = [10, 15, 20, 25]  # Start with a few key degrees
        # Add more if available
        additional_degrees = [d for d in all_degrees if d % 5 == 0 and d not in selected_degrees]
        selected_degrees.extend(additional_degrees[:6])  # Add up to 6 more
        selected_degrees = sorted(selected_degrees)
    
    plt.figure(figsize=(12, 8))
    
    # Colors for different degrees
    colors = plt.cm.tab20(np.linspace(0, 1, len(selected_degrees)))
    
    fit_results = []
    
    for i, degree in enumerate(selected_degrees):
        # Get data for this degree
        degree_data = [d for d in all_data if d['degree'] == degree]
        
        if len(degree_data) < 3:  # Need at least 3 points
            continue
            
        braid_lengths = np.array([d['braid_length'] for d in degree_data])
        times = np.array([d['time_seconds'] for d in degree_data])
        
        # Sort by braid length
        sort_idx = np.argsort(braid_lengths)
        braid_lengths = braid_lengths[sort_idx]
        times = times[sort_idx]
        
        color = colors[i]
        
        # Plot measured points
        plt.scatter(braid_lengths, times, color=color, s=80, alpha=0.8, 
                   label=f'Degree {degree} (measured)', zorder=3)
        
        # Try to fit exponential
        try:
            # Fit exponential model
            popt, pcov = curve_fit(exponential_model, braid_lengths, times, 
                                 p0=[1, 0.1, 0], maxfev=10000)
            
            a, b, c = popt
            
            # Calculate R²
            predicted = exponential_model(braid_lengths, a, b, c)
            r2 = 1 - np.sum((times - predicted)**2) / np.sum((times - np.mean(times))**2)
            
            # Generate smooth curve for plotting
            x_smooth = np.linspace(braid_lengths.min() - 0.5, 
                                 braid_lengths.max() + 2, 100)
            y_smooth = exponential_model(x_smooth, a, b, c)
            
            # Plot fitted curve
            plt.plot(x_smooth, y_smooth, color=color, linestyle='-', 
                    linewidth=2, alpha=0.7, 
                    label=f'Degree {degree} fit (R²={r2:.3f})')
            
            fit_results.append({
                'degree': degree,
                'a': a, 'b': b, 'c': c,
                'r_squared': r2,
                'formula': f't = {a:.2e}×e^({b:.3f}×n) + {c:.2f}'
            })
            
        except Exception as e:
            print(f"Failed to fit degree {degree}: {e}")
            # Just connect the points with lines
            plt.plot(braid_lengths, times, color=color, linestyle='--', 
                    linewidth=1, alpha=0.5, label=f'Degree {degree} (linear)')
    
    plt.xlabel('Braid Length', fontsize=12)
    plt.ylabel('Computation Time (seconds)', fontsize=12)
    plt.title('Sign Diagram Computation Time vs Braid Length\n(Exponential Fits)', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.yscale('log')  # Log scale for time since it grows exponentially
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig('time_vs_braid_length.png', dpi=300, bbox_inches='tight')
    print("Saved: time_vs_braid_length.png")
    plt.show()
    
    # Print fit results
    print("\n# Exponential Fit Results")
    print("# Model: time = a × exp(b × braid_length) + c")
    print()
    print("| Degree | Coefficient a | Exp Rate b | Offset c | R² | Formula |")
    print("|--------|---------------|------------|----------|-----|---------|")
    
    for result in fit_results:
        print(f"| {result['degree']:6d} | {result['a']:13.3e} | {result['b']:10.6f} | {result['c']:8.2f} | {result['r_squared']:4.3f} | {result['formula']} |")

def main():
    # Load all timing data
    files = [
        ('forecast_table_braid_8.csv', 8),
        ('forecast_table_braid_10.csv', 10),
        ('forecast_table_braid_12.csv', 12),
        ('forecast_table_braid_14.csv', 14)
    ]
    
    all_data = []
    for filename, braid_length in files:
        try:
            data = load_timing_data(filename, braid_length)
            all_data.extend(data)
            print(f"Loaded {len(data)} measurements from {filename}")
        except FileNotFoundError:
            print(f"File not found: {filename}")
    
    if not all_data:
        print("No data loaded. Please check file paths.")
        return
    
    print(f"Total measurements loaded: {len(all_data)}")
    
    # Create the plot
    plot_time_vs_braid_length(all_data)

if __name__ == "__main__":
    main()