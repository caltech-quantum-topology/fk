#!/usr/bin/env python3
"""
Forecast Table Launcher
Easy launcher for different forecast table viewing interfaces.
"""

import argparse
import sys
import os

def main():
    print("🚀 Profiling Forecast Table Launcher")
    print("=====================================")
    
    parser = argparse.ArgumentParser(description='Launch forecast table viewers')
    parser.add_argument('interface', choices=['cli', 'gui', 'web'], 
                       help='Interface type: cli (command line), gui (desktop), web (browser)')
    parser.add_argument('--input', '-i', default='profiling_results.txt',
                       help='Input profiling file (default: profiling_results.txt)')
    parser.add_argument('--port', '-p', type=int, default=8080,
                       help='Web server port for web interface (default: 8080)')
    
    # Additional CLI options
    parser.add_argument('--cpu-freq', '-f', type=float, default=3.2,
                       help='CPU frequency in GHz (default: 3.2)')
    parser.add_argument('--degree-min', type=int, default=5,
                       help='Minimum degree (default: 5)')
    parser.add_argument('--degree-max', type=int, default=50,
                       help='Maximum degree (default: 50)')
    parser.add_argument('--braid-length', '-b', type=int,
                       help='Specific braid length for CLI mode')
    parser.add_argument('--save-csv', action='store_true',
                       help='Save CSV files (CLI mode)')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"❌ Error: Input file '{args.input}' not found!")
        print(f"   Please make sure the profiling results file exists.")
        sys.exit(1)
    
    print(f"📁 Using input file: {args.input}")
    
    if args.interface == 'cli':
        print("🖥️  Launching Command Line Interface...")
        
        cmd_args = [
            'python', 'forecast_table_generator.py',
            '--input', args.input,
            '--cpu-freq', str(args.cpu_freq),
            '--degree-min', str(args.degree_min),
            '--degree-max', str(args.degree_max)
        ]
        
        if args.braid_length:
            cmd_args.extend(['--braid-length', str(args.braid_length)])
        
        if args.save_csv:
            cmd_args.append('--save-csv')
        
        os.execvp('python', cmd_args)
    
    elif args.interface == 'gui':
        print("🖼️  Launching Desktop GUI...")
        print("   Note: Make sure you have tkinter installed (usually comes with Python)")
        
        os.execvp('python', ['python', 'forecast_table_gui.py'])
    
    elif args.interface == 'web':
        print("🌐 Launching Web Interface...")
        print(f"   Server will start on http://localhost:{args.port}")
        print("   Your browser should open automatically")
        print("   Press Ctrl+C to stop the server")
        
        os.execvp('python', ['python', 'forecast_web_viewer.py', 
                            '--input', args.input, 
                            '--port', str(args.port)])

def show_help():
    print("""
🚀 Profiling Forecast Table Launcher Help
==========================================

This launcher provides three different ways to view forecast tables:

1. CLI (Command Line Interface):
   python forecast_launcher.py cli
   - Text-based tables printed to terminal
   - Good for scripting and automation
   - Can save CSV files
   - Options: --braid-length, --save-csv, --cpu-freq, --degree-min, --degree-max

2. GUI (Desktop Application):
   python forecast_launcher.py gui
   - Desktop application with interactive tables and plots
   - Real-time parameter adjustment
   - Visual plots with exponential curve fitting
   - Requires tkinter (usually included with Python)

3. Web (Browser Interface):
   python forecast_launcher.py web
   - Modern web interface accessible via browser
   - Interactive table with filtering
   - Responsive design works on mobile
   - No additional dependencies required

Examples:
  python forecast_launcher.py cli --braid-length 8 --save-csv
  python forecast_launcher.py gui
  python forecast_launcher.py web --port 8080
  python forecast_launcher.py cli --degree-min 20 --degree-max 30

All interfaces support:
- Multiple braid lengths
- Configurable CPU frequency (default: 3.2 GHz for M1 2020)
- Degree range selection
- Clock cycle calculations
- Exponential curve fitting with R² statistics
""")

if __name__ == "__main__":
    if len(sys.argv) == 1 or sys.argv[1] in ['--help', '-h', 'help']:
        show_help()
    else:
        main()