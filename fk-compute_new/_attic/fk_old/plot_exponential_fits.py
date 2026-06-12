import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def exponential_func(x, a, b):
    return a * np.exp(b * x)

data = []
with open('profiling_results.txt', 'r') as f:
    for line in f:
        if 'Time taken for sign assignment at length' in line:
            match = re.search(r'length (\d+): ([\d.]+) seconds', line)
            if match:
                length = int(match.group(1))
                time = float(match.group(2))
                data.append((length, time))

if not data:
    print("No sign assignment data found")
    exit()

data.sort(key=lambda x: x[0])
lengths = np.array([x[0] for x in data])
times = np.array([x[1] for x in data])

print(f"Found {len(data)} sign assignment timing points:")
for length, time in data:
    print(f"Length {length}: {time:.2f} seconds")

try:
    popt, pcov = curve_fit(exponential_func, lengths, times, p0=[1, 0.1])
    a_fit, b_fit = popt
    
    perr = np.sqrt(np.diag(pcov))
    a_err, b_err = perr
    
    print(f"\nExponential fit: y = {a_fit:.4f} * exp({b_fit:.4f} * x)")
    print(f"Parameter errors: a ± {a_err:.4f}, b ± {b_err:.4f}")
    
    r_squared = 1 - np.sum((times - exponential_func(lengths, *popt))**2) / np.sum((times - np.mean(times))**2)
    print(f"R-squared: {r_squared:.4f}")
    
    x_fit = np.linspace(min(lengths), max(lengths), 100)
    y_fit = exponential_func(x_fit, *popt)
    
    plt.figure(figsize=(10, 6))
    plt.scatter(lengths, times, color='red', s=100, label='Data points', zorder=3)
    plt.plot(x_fit, y_fit, 'b-', linewidth=2, 
             label=f'Exponential fit: y = {a_fit:.3f}*exp({b_fit:.3f}*x)\nR² = {r_squared:.4f}')
    
    plt.xlabel('Braid Length')
    plt.ylabel('Time (seconds)')
    plt.title('Exponential Fit to Sign Assignment Timing Data')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    
    plt.tight_layout()
    plt.savefig('exponential_fit_plot.png', dpi=300, bbox_inches='tight')
    plt.show()
    
except Exception as e:
    print(f"Error fitting exponential: {e}")
    
    plt.figure(figsize=(10, 6))
    plt.scatter(lengths, times, color='red', s=100, label='Data points')
    plt.xlabel('Braid Length')
    plt.ylabel('Time (seconds)')
    plt.title('Sign Assignment Timing Data (No fit available)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig('timing_data_plot.png', dpi=300, bbox_inches='tight')
    plt.show()
