# Profiling Analysis and Forecasting Tools

A comprehensive suite of tools for analyzing profiling results, fitting exponential curves, and forecasting computation times with clock cycle calculations.

## 🚀 Quick Start

### Option 1: Use the Launcher (Recommended)
```bash
# Command line interface
python forecast_launcher.py cli

# Desktop GUI application
python forecast_launcher.py gui

# Web browser interface
python forecast_launcher.py web
```

### Option 2: Direct Tool Usage
```bash
# Command line analysis
python profiling_analyzer.py

# Desktop GUI
python forecast_table_gui.py

# Web interface
python forecast_web_viewer.py
```

## 📁 Files Overview

### Core Analysis Tools
- **`profiling_analyzer.py`** - Main command-line analysis tool with plotting
- **`profiling_utils.py`** - Core utility functions and libraries
- **`forecast_table_generator.py`** - Comprehensive table generation

### GUI Interfaces
- **`forecast_table_gui.py`** - Desktop application with interactive tables and plots
- **`forecast_web_viewer.py`** - Modern web interface accessible via browser
- **`forecast_launcher.py`** - Easy launcher for all interfaces

### Generated Files
- **`profiling_analysis.png`** - Comprehensive analysis plots
- **`forecast_table_braid_*.csv`** - Detailed forecast tables in CSV format
- **`braid_length_*_forecast.png`** - Individual braid length plots

## 🔧 Features

### Exponential Curve Fitting
- Uses model: `f(x) = a * exp(b * x) + c`
- High accuracy fitting (R² > 0.999 for provided data)
- Automatic parameter estimation with error handling

### Clock Cycle Calculations
- Optimized for M1 2020 MacBook Air (3.2 GHz)
- Configurable CPU frequency
- Human-readable formatting (K/M/G/T/P scales)

### Multiple Interfaces
1. **Command Line**: Fast, scriptable, CSV export
2. **Desktop GUI**: Interactive tables, real-time plots, parameter adjustment
3. **Web Interface**: Modern responsive design, accessible anywhere

### Comprehensive Forecasting
- Degrees 5-50 (configurable range)
- Mixed measured/predicted data
- Time formatting (seconds/minutes/hours/days)
- Extreme forecast warnings for impractical computations

## 📊 Data Format

Input file format (like `profiling_results.txt`):
```
Braid length: 8, Degree: 5, Time taken: 1.25 seconds, Max memory: 101.71875 MiB
Braid length: 8, Degree: 6, Time taken: 0.90 seconds, Max memory: 101.734375 MiB
...
```

## 💻 Usage Examples

### Command Line Interface
```bash
# Basic analysis with default settings
python profiling_analyzer.py

# Custom input file and forecasts
python profiling_analyzer.py -i my_data.txt -f 26 27 28 30

# Generate tables for specific range
python forecast_table_generator.py --degree-min 20 --degree-max 35 --save-csv

# Quiet mode with specific braid length
python forecast_table_generator.py -q -b 8 --cpu-freq 3.5
```

### Desktop GUI
- Load profiling files via file dialog
- Adjust parameters in real-time
- Interactive plots with zoom/pan
- Color-coded measured vs predicted data
- Export capabilities

### Web Interface
- Responsive design works on all devices
- Real-time table updates
- Beautiful modern UI
- No installation required beyond Python
- Accessible at `http://localhost:8080`

## 📈 Sample Results

For M1 2020 MacBook Air (3.2 GHz):

| Degree | Time      | Clock Cycles | Status    |
|--------|-----------|--------------|-----------|
| 25     | 3.66m     | 702.85G      | measured  |
| 30     | 26.91m    | 5.17T        | predicted |
| 35     | 3.30h     | 38.03T       | predicted |
| 40     | 1.01d     | 279.87T      | predicted |
| 50     | 54.83d    | 15.16P       | predicted |

## ⚠️ Forecast Warnings

For degrees > 40, computation times become astronomical:
- **Degree 50**: ~35-55 days (~10-15 petacycles)
- Beyond practical computation limits
- Use forecasts for planning and resource estimation only

## 🛠 Dependencies

```bash
pip install scipy matplotlib numpy
```

**Note**: `tkinter` is usually included with Python. For the web interface, only built-in Python libraries are used.

## 📱 Interface Comparison

| Feature | CLI | GUI | Web |
|---------|-----|-----|-----|
| Speed | ⚡⚡⚡ | ⚡⚡ | ⚡⚡ |
| Interactivity | ⚡ | ⚡⚡⚡ | ⚡⚡⚡ |
| Plots | ⚡⚡ | ⚡⚡⚡ | ⚡ |
| Automation | ⚡⚡⚡ | ⚡ | ⚡ |
| Accessibility | ⚡⚡ | ⚡⚡ | ⚡⚡⚡ |
| Dependencies | ⚡⚡ | ⚡⚡ | ⚡⚡⚡ |

## 🎯 Use Cases

- **Research**: Analyze algorithm complexity and scaling behavior
- **Planning**: Estimate computation time for higher degrees
- **Optimization**: Identify performance bottlenecks
- **Presentation**: Generate publication-ready plots and tables
- **Automation**: Integrate into CI/CD pipelines with CLI tools

## 🤝 Contributing

The tools are designed to be modular and extensible:
- Add new fitting functions in `profiling_utils.py`
- Customize GUI layouts in `forecast_table_gui.py`
- Extend web API in `forecast_web_viewer.py`
- Add new output formats in `forecast_table_generator.py`