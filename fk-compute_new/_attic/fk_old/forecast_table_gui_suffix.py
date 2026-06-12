#!/usr/bin/env python3
"""
Forecast Table GUI for CSV Suffix Format
Interactive GUI application for displaying runtime forecast tables with filtering and visualization.
Specialized for the CSV format in profiling_results2_suffix.txt files.
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import threading
from profiling_utils_suffix import (
    parse_profiling_data_csv,
    generate_forecast_table,
    exponential_func,
    fit_exponential_to_data
)

class ForecastTableGUISuffix:
    def __init__(self, root):
        self.root = root
        self.root.title("Profiling Forecast Table Viewer (CSV Suffix Format)")
        self.root.geometry("1200x800")

        # Data storage
        self.data = {}
        self.current_results = []
        self.current_braid_length = None

        # Create GUI components
        self.create_widgets()

        # Try to load default data
        self.load_default_data()

    def create_widgets(self):
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(2, weight=1)

        # Control panel
        self.create_control_panel(main_frame)

        # Notebook for tables and plots
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.grid(row=2, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(10, 0))

        # Table tab
        self.table_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.table_frame, text="Forecast Table")
        self.create_table_tab()

        # Plot tab
        self.plot_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.plot_frame, text="Visualization")
        self.create_plot_tab()

    def create_control_panel(self, parent):
        # File loading
        file_frame = ttk.LabelFrame(parent, text="Data File", padding="5")
        file_frame.grid(row=0, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 10))

        ttk.Button(file_frame, text="Load CSV Profiling File",
                  command=self.load_file).pack(side=tk.LEFT)

        self.file_label = ttk.Label(file_frame, text="No file loaded")
        self.file_label.pack(side=tk.LEFT, padx=(10, 0))

        # Settings frame
        settings_frame = ttk.LabelFrame(parent, text="Settings", padding="5")
        settings_frame.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 10))

        # Braid length selection
        ttk.Label(settings_frame, text="Braid Length:").grid(row=0, column=0, sticky=tk.W)
        self.braid_var = tk.StringVar()
        self.braid_combo = ttk.Combobox(settings_frame, textvariable=self.braid_var,
                                       state="readonly", width=10)
        self.braid_combo.grid(row=0, column=1, padx=(5, 10))
        self.braid_combo.bind('<<ComboboxSelected>>', self.on_braid_change)

        # Target CPU frequency
        ttk.Label(settings_frame, text="Target CPU (GHz):").grid(row=0, column=2, sticky=tk.W)
        self.target_cpu_freq_var = tk.StringVar(value="3.2")
        ttk.Entry(settings_frame, textvariable=self.target_cpu_freq_var, width=8).grid(row=0, column=3, padx=(5, 10))

        # Degree range
        ttk.Label(settings_frame, text="Degree Range:").grid(row=0, column=4, sticky=tk.W)
        self.deg_min_var = tk.StringVar(value="5")
        self.deg_max_var = tk.StringVar(value="50")
        ttk.Entry(settings_frame, textvariable=self.deg_min_var, width=5).grid(row=0, column=5, padx=(5, 2))
        ttk.Label(settings_frame, text="to").grid(row=0, column=6)
        ttk.Entry(settings_frame, textvariable=self.deg_max_var, width=5).grid(row=0, column=7, padx=(2, 10))

        # Update button
        ttk.Button(settings_frame, text="Update Table",
                  command=self.update_table).grid(row=0, column=8, padx=(10, 0))

    def create_table_tab(self):
        # Table with scrollbars
        table_container = ttk.Frame(self.table_frame)
        table_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Treeview for table
        columns = ('degree', 'ground_truth_time', 'ground_truth_cycles', 'target_time', 'status')
        self.tree = ttk.Treeview(table_container, columns=columns, show='headings', height=20)

        # Configure column headings
        self.tree.heading('degree', text='Degree')
        self.tree.heading('ground_truth_time', text='Ground Truth Time (M1@3.2GHz)')
        self.tree.heading('ground_truth_cycles', text='Ground Truth Cycles')
        self.tree.heading('target_time', text='Target CPU Time')
        self.tree.heading('status', text='Status')

        # Configure column widths
        self.tree.column('degree', width=80, anchor=tk.CENTER)
        self.tree.column('ground_truth_time', width=200, anchor=tk.CENTER)
        self.tree.column('ground_truth_cycles', width=180, anchor=tk.CENTER)
        self.tree.column('target_time', width=150, anchor=tk.CENTER)
        self.tree.column('status', width=100, anchor=tk.CENTER)

        # Scrollbars
        v_scrollbar = ttk.Scrollbar(table_container, orient=tk.VERTICAL, command=self.tree.yview)
        h_scrollbar = ttk.Scrollbar(table_container, orient=tk.HORIZONTAL, command=self.tree.xview)
        self.tree.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)

        # Pack table and scrollbars
        self.tree.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        v_scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        h_scrollbar.grid(row=1, column=0, sticky=(tk.W, tk.E))

        # Configure grid weights
        table_container.columnconfigure(0, weight=1)
        table_container.rowconfigure(0, weight=1)

        # Status and statistics
        stats_frame = ttk.LabelFrame(self.table_frame, text="Statistics", padding="5")
        stats_frame.pack(fill=tk.X, padx=10, pady=(0, 10))

        self.stats_label = ttk.Label(stats_frame, text="No data loaded")
        self.stats_label.pack(anchor=tk.W)

    def create_plot_tab(self):
        # Matplotlib figure
        self.fig, self.ax = plt.subplots(figsize=(10, 6))
        self.ax.set_xlabel('Degree')
        self.ax.set_ylabel('Time (seconds)')
        self.ax.set_yscale('log')
        self.ax.grid(True, alpha=0.3)

        # Embed plot in tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, self.plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Plot controls
        plot_controls = ttk.Frame(self.plot_frame)
        plot_controls.pack(fill=tk.X, padx=10, pady=(0, 10))

        ttk.Button(plot_controls, text="Refresh Plot",
                  command=self.update_plot).pack(side=tk.LEFT)

        self.show_forecast_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(plot_controls, text="Show Forecasted Points",
                       variable=self.show_forecast_var).pack(side=tk.LEFT, padx=(10, 0))

    def load_default_data(self):
        """Try to load default profiling_results2_suffix.txt file"""
        try:
            self.data = parse_profiling_data_csv('profiling_results2_suffix.txt')
            if self.data:
                self.file_label.config(text="profiling_results2_suffix.txt")
                self.update_braid_combo()
                self.update_table()
        except:
            pass

    def load_file(self):
        """Load profiling data from CSV file"""
        filename = filedialog.askopenfilename(
            title="Select CSV Profiling Results File",
            filetypes=[("CSV files", "*.csv"), ("Text files", "*.txt"), ("All files", "*.*")]
        )

        if filename:
            try:
                self.data = parse_profiling_data_csv(filename)
                if self.data:
                    self.file_label.config(text=filename.split('/')[-1])
                    self.update_braid_combo()
                    self.update_table()
                else:
                    messagebox.showerror("Error", "No valid data found in the file")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file: {str(e)}")

    def update_braid_combo(self):
        """Update braid length combo box"""
        braid_lengths = sorted(self.data.keys())
        self.braid_combo['values'] = braid_lengths
        if braid_lengths:
            self.braid_var.set(str(braid_lengths[0]))
            self.current_braid_length = braid_lengths[0]

    def on_braid_change(self, event=None):
        """Handle braid length selection change"""
        try:
            self.current_braid_length = int(self.braid_var.get())
            self.update_table()
        except ValueError:
            pass

    def update_table(self):
        """Update the forecast table"""
        if not self.data or self.current_braid_length is None:
            return

        try:
            # Get parameters
            target_cpu_freq = float(self.target_cpu_freq_var.get())
            deg_min = int(self.deg_min_var.get())
            deg_max = int(self.deg_max_var.get())

            # Get data for current braid length
            degrees, times = zip(*self.data[self.current_braid_length])

            # Generate forecast table
            degree_range = range(deg_min, deg_max + 1)
            results, r_squared = generate_forecast_table(degrees, times, degree_range, target_cpu_freq, 3.2)

            self.current_results = results

            # Clear existing table
            for item in self.tree.get_children():
                self.tree.delete(item)

            # Populate table
            for result in results:
                # Color code measured vs predicted
                tags = ('measured',) if result['status'] == 'measured' else ('predicted',)

                self.tree.insert('', tk.END,
                               values=(result['degree'],
                                      result['ground_truth_formatted_time'],
                                      result['ground_truth_formatted_cycles'],
                                      result['predicted_formatted_time'],
                                      result['status']),
                               tags=tags)

            # Configure row colors
            self.tree.tag_configure('measured', background='#e8f5e8')
            self.tree.tag_configure('predicted', background='#f0f8ff')

            # Update statistics
            measured_count = sum(1 for r in results if r['status'] == 'measured')
            predicted_count = len(results) - measured_count

            stats_text = (f"Braid Length: {self.current_braid_length} | "
                         f"Target CPU: {target_cpu_freq} GHz | "
                         f"R² = {r_squared:.6f} | "
                         f"Measured: {measured_count} | "
                         f"Predicted: {predicted_count} | "
                         f"Total: {len(results)} degrees")

            self.stats_label.config(text=stats_text)

            # Update plot
            self.update_plot()

        except Exception as e:
            messagebox.showerror("Error", f"Failed to update table: {str(e)}")

    def update_plot(self):
        """Update the visualization plot"""
        if not self.data or self.current_braid_length is None:
            return

        try:
            self.ax.clear()

            # Get data
            degrees, times = zip(*self.data[self.current_braid_length])
            degrees = np.array(degrees)
            times = np.array(times)

            # Plot measured data
            self.ax.scatter(degrees, times, color='blue', alpha=0.7, s=50,
                           label='Measured Data', zorder=3)

            # Fit curve and plot
            params, r_squared = fit_exponential_to_data(degrees, times)

            if params is not None:
                # Plot fitted curve
                deg_min = int(self.deg_min_var.get())
                deg_max = int(self.deg_max_var.get())
                x_range = np.linspace(min(degrees), max(deg_max, max(degrees)), 200)
                y_fit = exponential_func(x_range, *params)

                self.ax.plot(x_range, y_fit, 'r-', linewidth=2,
                           label=f'Exponential Fit (R²={r_squared:.4f})', zorder=2)

                # Plot forecast points if enabled
                if self.show_forecast_var.get() and self.current_results:
                    forecast_degrees = []
                    forecast_times = []

                    for result in self.current_results:
                        if (result['status'] == 'predicted' and
                            result['degree'] > max(degrees)):
                            forecast_degrees.append(result['degree'])
                            forecast_times.append(result['ground_truth_time_seconds'])

                    if forecast_degrees:
                        self.ax.scatter(forecast_degrees, forecast_times,
                                      color='red', s=60, marker='^',
                                      label='Forecasted Points', alpha=0.8, zorder=3)

                # Add equation
                a, b, c = params
                equation = f'f(x) = {a:.3f}·exp({b:.3f}·x) + {c:.3f}'
                self.ax.text(0.05, 0.95, equation, transform=self.ax.transAxes,
                           fontsize=10, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

            # Formatting
            self.ax.set_xlabel('Degree')
            self.ax.set_ylabel('Time (seconds)')
            self.ax.set_title(f'Runtime Analysis - Braid Length {self.current_braid_length} (CSV Format)')
            self.ax.set_yscale('log')
            self.ax.grid(True, alpha=0.3)
            self.ax.legend()

            # Refresh canvas
            self.canvas.draw()

        except Exception as e:
            messagebox.showerror("Error", f"Failed to update plot: {str(e)}")

def main():
    root = tk.Tk()
    app = ForecastTableGUISuffix(root)
    root.mainloop()

if __name__ == "__main__":
    main()