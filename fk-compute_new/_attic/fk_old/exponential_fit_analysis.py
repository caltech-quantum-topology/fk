import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
import csv

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

def fit_exponential_by_degree(all_data):
    """Fit exponential models for each degree across braid lengths"""
    degrees = sorted(set(d['degree'] for d in all_data))
    
    results = []
    
    for degree in degrees:
        degree_data = [d for d in all_data if d['degree'] == degree]
        braid_lengths = np.array([d['braid_length'] for d in degree_data])
        times = np.array([d['time_seconds'] for d in degree_data])
        
        if len(braid_lengths) >= 3:  # Need at least 3 points for fitting
            try:
                # Fit exponential model
                popt, pcov = curve_fit(exponential_model, braid_lengths, times, 
                                     p0=[0.1, 0.5, 0], maxfev=10000)
                
                a, b, c = popt
                
                # Calculate R²
                predicted = exponential_model(braid_lengths, a, b, c)
                r2 = 1 - np.sum((times - predicted)**2) / np.sum((times - np.mean(times))**2)
                
                # Calculate correlation coefficient
                corr_coef = pearsonr(braid_lengths, times)[0]
                
                results.append({
                    'degree': degree,
                    'a': a,
                    'b': b, 
                    'c': c,
                    'r_squared': r2,
                    'correlation': corr_coef,
                    'n_points': len(braid_lengths),
                    'braid_lengths': list(braid_lengths),
                    'measured_times': list(times),
                    'predicted_times': list(predicted)
                })
            except Exception as e:
                print(f"Failed to fit degree {degree}: {e}")
    
    return results

def create_summary_table(fit_results):
    """Create a summary table of exponential fits"""
    table_data = []
    
    print("# Exponential Fit Analysis: Sign Diagram Computation Time vs Braid Length")
    print("# Model: time = a * exp(b * braid_length) + c")
    print()
    print("| Degree | a (coeff) | b (exp rate) | c (offset) | R² | Corr | Points | Formula |")
    print("|--------|-----------|--------------|------------|----|----- |--------|---------|")
    
    for result in fit_results:
        degree = result['degree']
        a = result['a']
        b = result['b']
        c = result['c']
        r2 = result['r_squared']
        corr = result['correlation']
        n_points = result['n_points']
        
        formula = f"t = {a:.3e}×e^({b:.3f}×n) + {c:.3f}"
        
        print(f"| {degree:6d} | {a:9.3e} | {b:12.6f} | {c:10.3f} | {r2:4.3f} | {corr:5.3f} | {n_points:6d} | {formula} |")
        
        table_data.append({
            'degree': degree,
            'coefficient_a': a,
            'exponential_rate_b': b,
            'offset_c': c,
            'r_squared': r2,
            'correlation': corr,
            'n_points': n_points,
            'formula': formula
        })
    
    return table_data

def create_detailed_predictions_table(fit_results, max_braid_length=20):
    """Create detailed predictions table"""
    print("\n# Detailed Time Predictions by Degree and Braid Length")
    print()
    
    # Get all unique degrees
    degrees = sorted([r['degree'] for r in fit_results])
    braid_lengths = list(range(8, max_braid_length + 1, 2))  # 8, 10, 12, ..., 20
    
    # Header
    header = "| Braid Length |"
    for degree in degrees:
        header += f" Deg {degree:2d} (s) |"
    print(header)
    
    # Separator
    separator = "|--------------|"
    for _ in degrees:
        separator += "-----------|"
    print(separator)
    
    # Data rows
    for braid_length in braid_lengths:
        row = f"| {braid_length:12d} |"
        
        for degree in degrees:
            # Find fit result for this degree
            result = next((r for r in fit_results if r['degree'] == degree), None)
            if result:
                a, b, c = result['a'], result['b'], result['c']
                predicted_time = exponential_model(braid_length, a, b, c)
                
                if predicted_time < 60:
                    row += f" {predicted_time:9.3f} |"
                elif predicted_time < 3600:
                    row += f" {predicted_time/60:7.1f}m |"
                elif predicted_time < 86400:
                    row += f" {predicted_time/3600:7.1f}h |"
                else:
                    row += f" {predicted_time/86400:7.1f}d |"
            else:
                row += "       N/A |"
        
        print(row)

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
        data = load_timing_data(filename, braid_length)
        all_data.extend(data)
    
    print(f"Loaded {len(all_data)} timing measurements")
    
    # Fit exponential models
    fit_results = fit_exponential_by_degree(all_data)
    
    # Create and display tables
    summary_data = create_summary_table(fit_results)
    create_detailed_predictions_table(fit_results)
    
    # Save results to CSV
    with open('exponential_fit_summary.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=summary_data[0].keys())
        writer.writeheader()
        writer.writerows(summary_data)
    
    print(f"\n# Summary saved to: exponential_fit_summary.csv")
    print(f"# Analyzed {len(fit_results)} degree levels")

if __name__ == "__main__":
    main()