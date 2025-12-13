#!/usr/bin/env python3
"""
Visualize benchmark results from res.csv
- Histograms for each degree showing time and memory distributions
- Growth plots for average and max times/memory vs degree
"""

import csv
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend for headless environments
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd

# Read the data
data = pd.read_csv("benchmark_results.csv")

print(f"Total rows: {len(data)}")
print(f"Successful runs: {len(data[data['success']])}")
print(f"Degrees in dataset: {data['degree'].nunique()}")

# Within each degree, remove the top outlier
data_cleaned = pd.DataFrame()
for degree in data["degree"].unique():
    subset = data[data["degree"] == degree]
    if len(subset) > 1:
        max_time = subset["time_seconds"].max()
        subset = subset[subset["time_seconds"] < max_time]
    data_cleaned = pd.concat([data_cleaned, subset], ignore_index=True)
data = data_cleaned

# ===== PART 1: Histograms for each degree =====
print("\nCreating histograms for each degree...")

degrees = sorted(data["degree"].unique())
# Calculate number of rows needed for subplots
n_degrees = len(degrees)
n_cols = 4
n_rows = (n_degrees + n_cols - 1) // n_cols

# Create figure for time histograms
fig1, axes1_raw = plt.subplots(n_rows, n_cols, figsize=(20, 5 * n_rows))
fig1.suptitle("Time Distribution by Degree", fontsize=16, fontweight="bold", y=0.995)
if n_rows * n_cols > 1:
    axes1 = axes1_raw.flatten()
else:
    axes1 = [axes1_raw]

# Create figure for memory histograms
fig2, axes2_raw = plt.subplots(n_rows, n_cols, figsize=(20, 5 * n_rows))
fig2.suptitle(
    "Peak Memory Distribution by Degree", fontsize=16, fontweight="bold", y=0.995
)
if n_rows * n_cols > 1:
    axes2 = axes2_raw.flatten()
else:
    axes2 = [axes2_raw]

for idx, degree in enumerate(degrees):
    times = data[data["degree"] == degree]["time_seconds"].values
    memories = data[data["degree"] == degree]["peak_memory_mb"].values

    # Time histogram
    ax1 = axes1[idx]
    ax1.hist(times, bins=20, edgecolor="black", alpha=0.7, color="steelblue")
    ax1.set_xlabel("Time (seconds)", fontsize=10)
    ax1.set_ylabel("Frequency", fontsize=10)
    ax1.set_title(
        f"Degree {degree}\n(n={len(times)}, mean={times.mean():.4f}s, max={times.max():.4f}s)",
        fontsize=11,
    )
    ax1.grid(True, alpha=0.3)

    # Memory histogram
    ax2 = axes2[idx]
    ax2.hist(memories, bins=20, edgecolor="black", alpha=0.7, color="coral")
    ax2.set_xlabel("Peak Memory (MB)", fontsize=10)
    ax2.set_ylabel("Frequency", fontsize=10)
    ax2.set_title(
        f"Degree {degree}\n(n={len(memories)}, mean={memories.mean():.2f}MB, max={memories.max():.2f}MB)",
        fontsize=11,
    )
    ax2.grid(True, alpha=0.3)

# Hide unused subplots
for idx in range(n_degrees, len(axes1)):
    axes1[idx].set_visible(False)
    axes2[idx].set_visible(False)

plt.tight_layout()
fig1.savefig("histograms_time.png", dpi=300, bbox_inches="tight")
fig2.savefig("histograms_memory.png", dpi=300, bbox_inches="tight")
print("Saved: histograms_time.png, histograms_memory.png")

# ===== PART 2: Growth Analysis by Degree =====
print("\nCreating growth analysis plots...")

# Calculate statistics by degree
degrees_array = np.array(degrees)
time_means = np.array(
    [np.mean(data[data["degree"] == d]["time_seconds"]) for d in degrees]
)
time_totals = np.array(
    [np.sum(data[data["degree"] == d]["time_seconds"]) for d in degrees]
)
time_maxs = np.array(
    [np.max(data[data["degree"] == d]["time_seconds"]) for d in degrees]
)
time_mins = np.array(
    [np.min(data[data["degree"] == d]["time_seconds"]) for d in degrees]
)
time_medians = np.array(
    [np.median(data[data["degree"] == d]["time_seconds"]) for d in degrees]
)
memory_means = np.array(
    [np.mean(data[data["degree"] == d]["peak_memory_mb"]) for d in degrees]
)
memory_maxs = np.array(
    [np.max(data[data["degree"] == d]["peak_memory_mb"]) for d in degrees]
)
memory_mins = np.array(
    [np.min(data[data["degree"] == d]["peak_memory_mb"]) for d in degrees]
)
memory_medians = np.array(
    [np.median(data[data["degree"] == d]["peak_memory_mb"]) for d in degrees]
)
counts = np.array([len(data[data["degree"] == d]) for d in degrees])

print("\nStatistics by degree:")
print(
    f"{'Degree':<10} {'Count':<8} {'Time (mean)':<15} {'Time (max)':<15} {'Time (total)':<15} {'Memory (mean)':<15} {'Memory (max)':<15}"
)
print("-" * 110)
for i, deg in enumerate(degrees):
    count = counts[i]
    print(
        f"{deg:<10} {count:<8} {time_means[i]:<15.4f} {time_maxs[i]:<15.4f} {time_totals[i]:<15.2f} {memory_means[i]:<15.2f} {memory_maxs[i]:<15.2f}"
    )

# Create growth analysis plots - now with 3 columns to include total time
fig3, axes = plt.subplots(2, 3, figsize=(24, 12))
ax_time_mean = axes[0, 0]
ax_time_max = axes[0, 1]
ax_time_total = axes[0, 2]
ax_mem_mean = axes[1, 0]
ax_mem_max = axes[1, 1]
ax_time_total_poly = axes[1, 2]
fig3.suptitle("Computational Growth Analysis by Degree", fontsize=16, fontweight="bold")


# Helper function to try multiple fits and choose the best
def fit_and_predict(x, y, x_pred, fit_types=["exponential", "polynomial"]):
    """Try different fits and return the best one based on R²"""
    best_fit = None
    best_r2 = -np.inf
    best_pred = None
    best_name = None

    # Exponential fit: y = a * e^(b*x)
    if "exponential" in fit_types and np.all(y > 0):
        try:
            log_y = np.log(y)
            coeffs = np.polyfit(x, log_y, 1)
            y_fit = np.exp(coeffs[1]) * np.exp(coeffs[0] * x)
            residuals = y - y_fit
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r2 = 1 - (ss_res / ss_tot)

            if r2 > best_r2:
                best_r2 = r2
                best_pred = np.exp(coeffs[1]) * np.exp(coeffs[0] * x_pred)
                best_name = f"Exponential (R²={r2:.3f})"
                best_fit = lambda x: np.exp(coeffs[1]) * np.exp(coeffs[0] * x)
        except:
            pass

    # Polynomial fit (degree 2)
    if "polynomial" in fit_types:
        try:
            coeffs = np.polyfit(x, y, 2)
            y_fit = np.polyval(coeffs, x)
            residuals = y - y_fit
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r2 = 1 - (ss_res / ss_tot)

            if r2 > best_r2:
                best_r2 = r2
                best_pred = np.polyval(coeffs, x_pred)
                best_name = f"Polynomial deg=2 (R²={r2:.3f})"
                best_fit = lambda x: np.polyval(coeffs, x)
        except:
            pass

    # Power law fit: y = a * x^b
    if "power" in fit_types and np.all(y > 0) and np.all(x > 0):
        try:
            log_x = np.log(x)
            log_y = np.log(y)
            coeffs = np.polyfit(log_x, log_y, 1)
            y_fit = np.exp(coeffs[1]) * (x ** coeffs[0])
            residuals = y - y_fit
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r2 = 1 - (ss_res / ss_tot)

            if r2 > best_r2:
                best_r2 = r2
                best_pred = np.exp(coeffs[1]) * (x_pred ** coeffs[0])
                best_name = f"Power Law (R²={r2:.3f})"
                best_fit = lambda x: np.exp(coeffs[1]) * (x ** coeffs[0])
        except:
            pass

    return best_fit, best_pred, best_name


# Prediction range
max_degree = degrees_array.max()
pred_degrees = np.linspace(degrees_array.min(), 40, 100)

# ---- Time Mean Plot ----
ax_time_mean.scatter(
    degrees_array,
    time_means,
    label="Mean Time",
    s=100,
    alpha=0.7,
    color="blue",
    zorder=3,
)
fit_func, pred_vals, fit_name = fit_and_predict(
    degrees_array, time_means, pred_degrees, ["exponential", "polynomial", "power"]
)
if fit_func is not None:
    ax_time_mean.plot(
        pred_degrees,
        pred_vals,
        "--",
        label=fit_name,
        linewidth=2,
        color="darkblue",
        alpha=0.8,
    )
ax_time_mean.axvline(max_degree, color="gray", linestyle=":", alpha=0.5, linewidth=1.5)
ax_time_mean.set_xlabel("Degree", fontsize=12)
ax_time_mean.set_ylabel("Mean Time (seconds)", fontsize=12)
ax_time_mean.set_title("Mean Computation Time Growth", fontsize=13, fontweight="bold")
ax_time_mean.legend()
ax_time_mean.grid(True, alpha=0.3)


# ---- Time Mean Plot ----
ax_time_mean.scatter(
    degrees_array,
    time_means,
    label="Mean Time",
    s=100,
    alpha=0.7,
    color="blue",
    zorder=3,
)
fit_func, pred_vals, fit_name = fit_and_predict(
    degrees_array, time_means, pred_degrees, ["exponential", "polynomial", "power"]
)
if fit_func is not None:
    ax_time_mean.plot(
        pred_degrees,
        pred_vals,
        "--",
        label=fit_name,
        linewidth=2,
        color="darkblue",
        alpha=0.8,
    )
ax_time_mean.axvline(max_degree, color="gray", linestyle=":", alpha=0.5, linewidth=1.5)
ax_time_mean.set_xlabel("Degree", fontsize=12)
ax_time_mean.set_ylabel("Mean Time (seconds)", fontsize=12)
ax_time_mean.set_title("Mean Computation Time Growth", fontsize=13, fontweight="bold")
ax_time_mean.legend()
ax_time_mean.grid(True, alpha=0.3)

# ---- Time Max Plot ----
ax_time_max.scatter(
    degrees_array, time_maxs, label="Max Time", s=100, alpha=0.7, color="red", zorder=3
)
fit_func, pred_vals, fit_name = fit_and_predict(
    degrees_array, time_maxs, pred_degrees, ["exponential", "polynomial", "power"]
)
if fit_func is not None:
    ax_time_max.plot(
        pred_degrees,
        pred_vals,
        "--",
        label=fit_name,
        linewidth=2,
        color="darkred",
        alpha=0.8,
    )
ax_time_max.axvline(max_degree, color="gray", linestyle=":", alpha=0.5, linewidth=1.5)
ax_time_max.set_xlabel("Degree", fontsize=12)
ax_time_max.set_ylabel("Max Time (seconds)", fontsize=12)
ax_time_max.set_title("Maximum Computation Time Growth", fontsize=13, fontweight="bold")
ax_time_max.legend()
ax_time_max.grid(True, alpha=0.3)

# ---- Time Total Plot ----
ax_time_total.scatter(
    degrees_array,
    time_totals,
    label="Total Time",
    s=100,
    alpha=0.7,
    color="purple",
    zorder=3,
)
fit_func, pred_vals, fit_name = fit_and_predict(
    degrees_array, time_totals, pred_degrees, ["exponential"]
)
if fit_func is not None:
    ax_time_total.plot(
        pred_degrees,
        pred_vals,
        "--",
        label=fit_name,
        linewidth=2,
        color="indigo",
        alpha=0.8,
    )
ax_time_total.axvline(max_degree, color="gray", linestyle=":", alpha=0.5, linewidth=1.5)
ax_time_total.set_yscale("log")
ax_time_total.set_xlabel("Degree", fontsize=12)
ax_time_total.set_ylabel("Total Time (seconds)", fontsize=12)
ax_time_total.set_title("Total Computation Time Growth", fontsize=13, fontweight="bold")
ax_time_total.legend()
ax_time_total.grid(True, alpha=0.3)

# ---- Time Total Plot Polynomial ----
ax_time_total_poly.scatter(
    degrees_array,
    time_totals,
    label="Total Time",
    s=100,
    alpha=0.7,
    color="purple",
    zorder=3,
)
fit_func, pred_vals, fit_name = fit_and_predict(
    degrees_array, time_totals, pred_degrees, ["polynomial"]
)
if fit_func is not None:
    ax_time_total_poly.plot(
        pred_degrees,
        pred_vals,
        "--",
        label=fit_name,
        linewidth=2,
        color="indigo",
        alpha=0.8,
    )
ax_time_total_poly.axvline(max_degree, color="gray", linestyle=":", alpha=0.5, linewidth=1.5)
ax_time_total_poly.set_xlabel("Degree", fontsize=12)
ax_time_total_poly.set_ylabel("Total Time (seconds)", fontsize=12)
ax_time_total_poly.set_title("Total Computation Time Growth", fontsize=13, fontweight="bold")
ax_time_total_poly.legend()
ax_time_total_poly.grid(True, alpha=0.3)

# ---- Memory Mean Plot ----
ax_mem_mean.scatter(
    degrees_array,
    memory_means,
    label="Mean Memory",
    s=100,
    alpha=0.7,
    color="green",
    zorder=3,
)
fit_func, pred_vals, fit_name = fit_and_predict(
    degrees_array, memory_means, pred_degrees, ["exponential", "polynomial", "power"]
)
if fit_func is not None:
    ax_mem_mean.plot(
        pred_degrees,
        pred_vals,
        "--",
        label=fit_name,
        linewidth=2,
        color="darkgreen",
        alpha=0.8,
    )
ax_mem_mean.axvline(max_degree, color="gray", linestyle=":", alpha=0.5, linewidth=1.5)
ax_mem_mean.set_xlabel("Degree", fontsize=12)
ax_mem_mean.set_ylabel("Mean Memory (MB)", fontsize=12)
ax_mem_mean.set_title("Mean Memory Usage Growth", fontsize=13, fontweight="bold")
ax_mem_mean.legend()
ax_mem_mean.grid(True, alpha=0.3)

# ---- Memory Max Plot ----
ax_mem_max.scatter(
    degrees_array,
    memory_maxs,
    label="Max Memory",
    s=100,
    alpha=0.7,
    color="orange",
    zorder=3,
)
fit_func, pred_vals, fit_name = fit_and_predict(
    degrees_array, memory_maxs, pred_degrees, ["exponential", "polynomial", "power"]
)
if fit_func is not None:
    ax_mem_max.plot(
        pred_degrees,
        pred_vals,
        "--",
        label=fit_name,
        linewidth=2,
        color="darkorange",
        alpha=0.8,
    )
ax_mem_max.axvline(max_degree, color="gray", linestyle=":", alpha=0.5, linewidth=1.5)
ax_mem_max.set_xlabel("Degree", fontsize=12)
ax_mem_max.set_ylabel("Max Memory (MB)", fontsize=12)
ax_mem_max.set_title("Maximum Memory Usage Growth", fontsize=13, fontweight="bold")
ax_mem_max.legend()
ax_mem_max.grid(True, alpha=0.3)

plt.tight_layout()
fig3.savefig("growth_analysis.png", dpi=300, bbox_inches="tight")
print("Saved: growth_analysis.png")

# ===== Print predictions =====
print("\n" + "=" * 80)
print("PREDICTIONS FOR HIGHER DEGREES")
print("=" * 80)

# Get fitted functions for predictions
fit_time_mean, _, _ = fit_and_predict(
    degrees_array, time_means, pred_degrees, ["exponential", "polynomial", "power"]
)
fit_time_max, _, _ = fit_and_predict(
    degrees_array, time_maxs, pred_degrees, ["exponential", "polynomial", "power"]
)
fit_time_total, _, _ = fit_and_predict(
    degrees_array, time_totals, pred_degrees, ["exponential"]
)
fit_mem_mean, _, _ = fit_and_predict(
    degrees_array, memory_means, pred_degrees, ["exponential", "polynomial", "power"]
)
fit_mem_max, _, _ = fit_and_predict(
    degrees_array, memory_maxs, pred_degrees, ["exponential", "polynomial", "power"]
)

# Target degrees for prediction
target_degrees = [20, 25, 30, 35, 40]

print("\nNote: Predictions are based on fitted models and may have high uncertainty,")
print("especially when extrapolating far beyond the observed data range.")
print()

for deg in target_degrees:
    print(f"\n{'='*60}")
    print(f"Degree {deg} Predictions:")
    print(f"{'='*60}")

    if fit_time_mean:
        pred_t_mean = fit_time_mean(deg)
        # Estimate total time assuming same number of samples as current average
        avg_count = np.mean(counts)
        total_time_est = pred_t_mean * avg_count
        print(
            f"  Time per computation (mean): {pred_t_mean:.4f} seconds ({pred_t_mean/60:.2f} minutes)"
        )
        print(
            f"  Estimated total time ({int(avg_count)} samples): {total_time_est:.2f} seconds ({total_time_est/60:.2f} minutes / {total_time_est/3600:.2f} hours)"
        )

    if fit_time_max:
        pred_t_max = fit_time_max(deg)
        print(
            f"  Time per computation (max):  {pred_t_max:.4f} seconds ({pred_t_max/60:.2f} minutes)"
        )

    if fit_time_total:
        pred_t_total = fit_time_total(deg)
        print(
            f"  Total time (predicted):       {pred_t_total:.2f} seconds ({pred_t_total/60:.2f} minutes / {pred_t_total/3600:.2f} hours)"
        )

    if fit_mem_mean:
        pred_m_mean = fit_mem_mean(deg)
        print(f"  Memory (mean): {pred_m_mean:.2f} MB ({pred_m_mean/1024:.3f} GB)")

    if fit_mem_max:
        pred_m_max = fit_mem_max(deg)
        print(f"  Memory (max):  {pred_m_max:.2f} MB ({pred_m_max/1024:.3f} GB)")

print("\n" + "=" * 80)
print("Visualization complete! Created 3 plots:")
print("  1. histograms_time.png - Time distribution for each degree")
print("  2. histograms_memory.png - Memory distribution for each degree")
print("  3. growth_analysis.png - Growth trends with predictions for degrees 20-40")
print("=" * 80)
