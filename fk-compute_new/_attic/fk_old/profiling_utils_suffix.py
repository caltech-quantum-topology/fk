"""
Utility functions for profiling analysis and forecasting - specialized for CSV suffix format.
Based on profiling_utils.py but adapted for the CSV format in profiling_results2_suffix.txt
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from collections import defaultdict

def exponential_func(x, a, b, c):
    """Exponential function: f(x) = a * exp(b * x) + c"""
    return a * np.exp(b * x) + c

def parse_profiling_data_csv(filename):
    """
    Parse CSV profiling results file and return structured data.
    
    Expected CSV format:
    Braid_Length,Degree,Sign_Assignment_Time,Sign_Assignment_Memory,Computation_Time,Computation_Memory,Knot_Index
    
    Args:
        filename (str): Path to CSV profiling results file
        
    Returns:
        dict: Dictionary with braid_length as keys, list of (degree, computation_time) tuples as values
    """
    data = defaultdict(list)
    
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                braid_length = int(row['Braid_Length'])
                degree = int(row['Degree'])
                computation_time = float(row['Computation_Time'])
                data[braid_length].append((degree, computation_time))
            except (ValueError, KeyError) as e:
                # Skip invalid rows
                continue
    
    # Sort each braid length's data by degree
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