"""
Utility functions for profiling analysis and forecasting.
"""

import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from collections import defaultdict

def exponential_func(x, a, b, c):
    """Exponential function: f(x) = a * exp(b * x) + c"""
    return a * np.exp(b * x) + c

def parse_profiling_data(filename):
    """
    Parse profiling results file and return structured data.
    
    Args:
        filename (str): Path to profiling results file
        
    Returns:
        dict: Dictionary with braid_length as keys, list of (degree, time) tuples as values
    """
    data = defaultdict(list)
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            match = re.search(r'Braid length: (\d+), Degree: (\d+), Time taken: ([\d.]+) seconds', line)
            if match:
                braid_length = int(match.group(1))
                degree = int(match.group(2))
                time = float(match.group(3))
                data[braid_length].append((degree, time))
    
    for braid_length in data:
        data[braid_length].sort(key=lambda x: x[0])
    
    return dict(data)

def fit_exponential_to_data(degrees, times):
    """
    Fit exponential curve to degree vs time data.
    
    Args:
        degrees (list): List of degrees
        times (list): List of corresponding times
        
    Returns:
        tuple: (parameters, r_squared) or (None, 0) if fitting fails
    """
    try:
        p0 = [1, 0.1, 0]
        popt, pcov = curve_fit(exponential_func, degrees, times, p0=p0, maxfev=5000)
        
        y_pred = exponential_func(np.array(degrees), *popt)
        ss_res = np.sum((np.array(times) - y_pred) ** 2)
        ss_tot = np.sum((np.array(times) - np.mean(times)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        
        return popt, r_squared
    except Exception as e:
        print(f"Failed to fit exponential curve: {e}")
        return None, 0

def forecast_time(degrees, times, target_degree):
    """
    Forecast time for a target degree based on existing data.
    
    Args:
        degrees (list): List of known degrees
        times (list): List of corresponding times
        target_degree (int): Degree to forecast
        
    Returns:
        tuple: (predicted_time, r_squared, equation_params) or (None, 0, None) if fitting fails
    """
    params, r_squared = fit_exponential_to_data(degrees, times)
    
    if params is not None:
        predicted_time = exponential_func(target_degree, *params)
        return predicted_time, r_squared, params
    
    return None, 0, None

def plot_with_forecast(degrees, times, braid_length, target_degrees=None, save_path=None):
    """
    Plot data with exponential fit and optional forecasts.
    
    Args:
        degrees (list): List of degrees
        times (list): List of times
        braid_length (int): Braid length for title
        target_degrees (list): Optional list of degrees to forecast
        save_path (str): Optional path to save plot
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    degrees = np.array(degrees)
    times = np.array(times)
    
    # Plot original data
    ax.scatter(degrees, times, color='blue', alpha=0.7, s=50, label='Measured Data')
    
    # Fit and plot exponential curve
    params, r_squared = fit_exponential_to_data(degrees, times)
    
    if params is not None:
        x_range = np.linspace(min(degrees), max(degrees), 100)
        if target_degrees:
            x_range = np.linspace(min(degrees), max(max(degrees), max(target_degrees)), 100)
        
        y_fit = exponential_func(x_range, *params)
        ax.plot(x_range, y_fit, 'r-', linewidth=2, label=f'Exponential Fit (R²={r_squared:.4f})')
        
        # Add forecasts if requested
        if target_degrees:
            forecast_times = [exponential_func(d, *params) for d in target_degrees]
            ax.scatter(target_degrees, forecast_times, color='red', s=80, 
                      marker='^', label='Forecasted Points')
            
            for deg, time in zip(target_degrees, forecast_times):
                ax.annotate(f'Deg {deg}\n{time:.1f}s', 
                           xy=(deg, time), xytext=(5, 5), 
                           textcoords='offset points', fontsize=8)
        
        # Add equation
        a, b, c = params
        equation = f'f(x) = {a:.3f}·exp({b:.3f}·x) + {c:.3f}'
        ax.text(0.05, 0.95, equation, transform=ax.transAxes, 
               fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax.set_xlabel('Degree')
    ax.set_ylabel('Time (seconds)')
    ax.set_title(f'Runtime Analysis - Braid Length {braid_length}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
    
    plt.show()

def calculate_clock_cycles(time_seconds, cpu_frequency_ghz=3.2):
    """
    Calculate approximate clock cycles for given time and CPU frequency.
    
    Args:
        time_seconds (float): Time in seconds
        cpu_frequency_ghz (float): CPU frequency in GHz (default: 3.2 for M1 2020)
        
    Returns:
        int: Approximate number of clock cycles
    """
    cycles_per_second = cpu_frequency_ghz * 1e9
    return int(time_seconds * cycles_per_second)

def calculate_ground_truth_cycles(time_seconds, reference_cpu_frequency_ghz=3.2):
    """
    Calculate ground truth clock cycles based on reference CPU (M1 2020 @ 3.2 GHz).
    These cycles represent the computational work and remain constant across architectures.
    
    Args:
        time_seconds (float): Ground truth time in seconds
        reference_cpu_frequency_ghz (float): Reference CPU frequency (default: 3.2 GHz for M1 2020)
        
    Returns:
        int: Architecture-independent computational cycles
    """
    return calculate_clock_cycles(time_seconds, reference_cpu_frequency_ghz)

def predict_time_for_cpu(ground_truth_cycles, target_cpu_frequency_ghz):
    """
    Predict execution time on a different CPU architecture.
    
    Args:
        ground_truth_cycles (int): Architecture-independent computational cycles
        target_cpu_frequency_ghz (float): Target CPU frequency in GHz
        
    Returns:
        float: Predicted execution time in seconds
    """
    cycles_per_second = target_cpu_frequency_ghz * 1e9
    return ground_truth_cycles / cycles_per_second

def format_time(seconds):
    """
    Format time in a human-readable way.
    
    Args:
        seconds (float): Time in seconds
        
    Returns:
        str: Formatted time string
    """
    if seconds < 60:
        return f"{seconds:.2f}s"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.2f}m"
    elif seconds < 86400:
        hours = seconds / 3600
        return f"{hours:.2f}h"
    else:
        days = seconds / 86400
        return f"{days:.2f}d"

def format_cycles(cycles):
    """
    Format clock cycles in a human-readable way.
    
    Args:
        cycles (int): Number of clock cycles
        
    Returns:
        str: Formatted cycles string
    """
    if cycles < 1e3:
        return f"{cycles}"
    elif cycles < 1e6:
        return f"{cycles/1e3:.2f}K"
    elif cycles < 1e9:
        return f"{cycles/1e6:.2f}M"
    elif cycles < 1e12:
        return f"{cycles/1e9:.2f}G"
    elif cycles < 1e15:
        return f"{cycles/1e12:.2f}T"
    else:
        return f"{cycles/1e15:.2f}P"

def generate_forecast_table(degrees, times, degree_range=None, target_cpu_frequency_ghz=3.2, reference_cpu_frequency_ghz=3.2):
    """
    Generate a comprehensive forecast table for a range of degrees with ground truth and predicted times.
    
    Args:
        degrees (list): List of known degrees
        times (list): List of corresponding ground truth times
        degree_range (range or list): Range of degrees to forecast (default: 5-50)
        target_cpu_frequency_ghz (float): Target CPU frequency for predictions (default: 3.2 GHz)
        reference_cpu_frequency_ghz (float): Reference CPU frequency for ground truth (default: 3.2 GHz for M1 2020)
        
    Returns:
        tuple: (results_list, r_squared)
    """
    if degree_range is None:
        degree_range = range(5, 51)
    
    params, r_squared = fit_exponential_to_data(degrees, times)
    
    if params is None:
        return [], 0
    
    results = []
    known_degrees = set(degrees)
    
    for degree in degree_range:
        if degree in known_degrees:
            # Use actual measured time as ground truth
            idx = degrees.index(degree)
            ground_truth_time = times[idx]
            status = "measured"
        else:
            # Use predicted time as ground truth (based on exponential fit)
            ground_truth_time = exponential_func(degree, *params)
            status = "predicted"
        
        # Calculate ground truth cycles (architecture-independent)
        ground_truth_cycles = calculate_ground_truth_cycles(ground_truth_time, reference_cpu_frequency_ghz)
        
        # Calculate predicted time for target CPU
        predicted_time = predict_time_for_cpu(ground_truth_cycles, target_cpu_frequency_ghz)
        
        results.append({
            'degree': degree,
            'ground_truth_time_seconds': ground_truth_time,
            'ground_truth_formatted_time': format_time(ground_truth_time),
            'ground_truth_cycles': ground_truth_cycles,
            'ground_truth_formatted_cycles': format_cycles(ground_truth_cycles),
            'predicted_time_seconds': predicted_time,
            'predicted_formatted_time': format_time(predicted_time),
            'status': status
        })
    
    return results, r_squared

def print_forecast_table(degrees, times, braid_length, degree_range=None, target_cpu_frequency_ghz=3.2, reference_cpu_frequency_ghz=3.2):
    """
    Print a formatted forecast table for a range of degrees with ground truth and predictions.
    
    Args:
        degrees (list): List of known degrees
        times (list): List of corresponding times
        braid_length (int): Braid length for the title
        degree_range (range or list): Range of degrees to forecast (default: 5-50)
        target_cpu_frequency_ghz (float): Target CPU frequency for predictions (default: 3.2 GHz)
        reference_cpu_frequency_ghz (float): Reference CPU frequency (default: 3.2 GHz for M1 2020)
    """
    results, r_squared = generate_forecast_table(degrees, times, degree_range, target_cpu_frequency_ghz, reference_cpu_frequency_ghz)
    
    if not results:
        print("Failed to generate forecast table")
        return
    
    print(f"\n{'='*110}")
    print(f"FORECAST TABLE - BRAID LENGTH {braid_length}")
    print(f"Reference CPU: {reference_cpu_frequency_ghz} GHz (M1 2020 MacBook Air) | Target CPU: {target_cpu_frequency_ghz} GHz")
    print(f"Fit Quality (R²): {r_squared:.6f}")
    print(f"{'='*110}")
    print(f"{'Degree':<8} {'Ground Truth':<25} {'Arch-Independent':<18} {'Target CPU Time':<15} {'Status':<10}")
    print(f"{'':8} {'Time (M1 @ 3.2GHz)':<25} {'Cycles':<18} {'(@ ' + f'{target_cpu_frequency_ghz}' + 'GHz)':<15}")
    print(f"{'-'*110}")
    
    for result in results:
        print(f"{result['degree']:<8} {result['ground_truth_formatted_time']:<25} "
              f"{result['ground_truth_formatted_cycles']:<18} "
              f"{result['predicted_formatted_time']:<15} {result['status']:<10}")

def save_forecast_table_csv(degrees, times, braid_length, filename, degree_range=None, target_cpu_frequency_ghz=3.2, reference_cpu_frequency_ghz=3.2):
    """
    Save forecast table to CSV file with ground truth and predictions.
    
    Args:
        degrees (list): List of known degrees
        times (list): List of corresponding times
        braid_length (int): Braid length
        filename (str): Output CSV filename
        degree_range (range or list): Range of degrees to forecast (default: 5-50)
        target_cpu_frequency_ghz (float): Target CPU frequency (default: 3.2 GHz)
        reference_cpu_frequency_ghz (float): Reference CPU frequency (default: 3.2 GHz for M1 2020)
    """
    results, r_squared = generate_forecast_table(degrees, times, degree_range, target_cpu_frequency_ghz, reference_cpu_frequency_ghz)
    
    if not results:
        print("Failed to generate forecast table")
        return
    
    with open(filename, 'w') as f:
        f.write(f"# Forecast Table - Braid Length {braid_length}\n")
        f.write(f"# Reference CPU: {reference_cpu_frequency_ghz} GHz (M1 2020)\n")
        f.write(f"# Target CPU: {target_cpu_frequency_ghz} GHz\n")
        f.write(f"# Fit Quality (R²): {r_squared:.6f}\n")
        f.write("Degree,Ground_Truth_Time_Seconds,Ground_Truth_Formatted_Time,Ground_Truth_Cycles,Ground_Truth_Formatted_Cycles,Predicted_Time_Seconds,Predicted_Formatted_Time,Status\n")
        
        for result in results:
            f.write(f"{result['degree']},{result['ground_truth_time_seconds']:.6f},"
                   f"{result['ground_truth_formatted_time']},{result['ground_truth_cycles']},"
                   f"{result['ground_truth_formatted_cycles']},{result['predicted_time_seconds']:.6f},"
                   f"{result['predicted_formatted_time']},{result['status']}\n")
    
    print(f"Forecast table saved to {filename}")

# Example usage
if __name__ == "__main__":
    # Example of how to use these functions
    data = parse_profiling_data('profiling_results.txt')
    
    print("=== DEMO: Architecture-Independent Forecasting ===")
    print("Ground truth measured on M1 2020 @ 3.2 GHz")
    print("Predictions for different CPU architectures\n")
    
    # Test with different CPU frequencies
    cpu_configs = [
        ("M1 2020 (Same)", 3.2),
        ("Intel i7-10700K", 3.8),
        ("AMD Ryzen 9 5950X", 4.9),
        ("Intel i9-13900K", 5.8),
        ("Hypothetical Future CPU", 8.0)
    ]
    
    # Show example for first braid length
    first_braid = sorted(data.keys())[0]
    degrees, times = zip(*data[first_braid])
    
    for cpu_name, cpu_freq in cpu_configs:
        print(f"\n🖥️  {cpu_name} @ {cpu_freq} GHz")
        print_forecast_table(degrees, times, first_braid, range(25, 31), cpu_freq, 3.2)
        
    # Generate CSV for comparison
    for cpu_name, cpu_freq in cpu_configs:
        safe_name = cpu_name.replace(" ", "_").replace("(", "").replace(")", "")
        csv_filename = f'forecast_braid_{first_braid}_{safe_name}_{cpu_freq}GHz.csv'
        save_forecast_table_csv(degrees, times, first_braid, csv_filename, range(25, 31), cpu_freq, 3.2)