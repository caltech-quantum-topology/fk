#!/usr/bin/env python3
"""
Web-based Forecast Table Viewer
Simple web interface for viewing forecast tables using built-in HTTP server.
"""

import http.server
import socketserver
import webbrowser
import json
import urllib.parse
import threading
import time
from profiling_utils import (
    parse_profiling_data, 
    generate_forecast_table
)

class ForecastWebHandler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, data=None, **kwargs):
        self.data = data or {}
        super().__init__(*args, **kwargs)
    
    def do_GET(self):
        if self.path == '/' or self.path == '/index.html':
            self.serve_main_page()
        elif self.path.startswith('/api/forecast'):
            self.serve_forecast_api()
        elif self.path.startswith('/api/braids'):
            self.serve_braids_api()
        else:
            super().do_GET()
    
    def serve_main_page(self):
        """Serve the main HTML page"""
        html = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Profiling Forecast Viewer</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 10px;
            text-align: center;
            margin-bottom: 20px;
        }
        .controls {
            background: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            margin-bottom: 20px;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            align-items: end;
        }
        .control-group {
            display: flex;
            flex-direction: column;
        }
        .control-group label {
            margin-bottom: 5px;
            font-weight: 600;
            color: #333;
        }
        .control-group select, .control-group input, .control-group button {
            padding: 8px 12px;
            border: 2px solid #ddd;
            border-radius: 5px;
            font-size: 14px;
        }
        .control-group button {
            background: #667eea;
            color: white;
            border: none;
            cursor: pointer;
            transition: background-color 0.3s;
        }
        .control-group button:hover {
            background: #5a67d8;
        }
        .stats {
            background: white;
            padding: 15px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            margin-bottom: 20px;
            text-align: center;
            font-weight: 600;
        }
        .table-container {
            background: white;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            overflow: hidden;
        }
        table {
            width: 100%;
            border-collapse: collapse;
        }
        th {
            background: #667eea;
            color: white;
            padding: 12px;
            text-align: center;
            font-weight: 600;
        }
        td {
            padding: 10px;
            text-align: center;
            border-bottom: 1px solid #eee;
        }
        tr.measured {
            background-color: #e8f5e8;
        }
        tr.predicted {
            background-color: #f0f8ff;
        }
        tr:hover {
            background-color: #f8f9fa;
        }
        .loading {
            text-align: center;
            padding: 40px;
            font-size: 18px;
            color: #666;
        }
        .error {
            color: #e53e3e;
            background: #fed7d7;
            padding: 15px;
            border-radius: 5px;
            margin: 10px 0;
        }
        .degree-cell {
            font-weight: bold;
            color: #2d3748;
        }
        .time-cell {
            font-family: 'Monaco', 'Menlo', monospace;
        }
        .cycles-cell {
            font-family: 'Monaco', 'Menlo', monospace;
            color: #4a5568;
        }
        .status-measured {
            color: #38a169;
            font-weight: bold;
        }
        .status-predicted {
            color: #3182ce;
            font-weight: bold;
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>🚀 Profiling Forecast Viewer</h1>
        <p>Interactive table for runtime and clock cycle forecasting</p>
    </div>
    
    <div class="controls">
        <div class="control-group">
            <label for="braidLength">Braid Length:</label>
            <select id="braidLength">
                <option value="">Loading...</option>
            </select>
        </div>
        
        <div class="control-group">
            <label for="targetCpuFreq">Target CPU Freq (GHz):</label>
            <input type="number" id="targetCpuFreq" value="3.2" step="0.1" min="0.1" max="10">
        </div>
        
        <div class="control-group">
            <label for="degreeMin">Min Degree:</label>
            <input type="number" id="degreeMin" value="5" min="1" max="100">
        </div>
        
        <div class="control-group">
            <label for="degreeMax">Max Degree:</label>
            <input type="number" id="degreeMax" value="40" min="1" max="100">
        </div>
        
        <div class="control-group">
            <label>&nbsp;</label>
            <button onclick="updateTable()">Update Table</button>
        </div>
    </div>
    
    <div class="stats" id="stats">
        Select parameters and click "Update Table" to view forecast data
    </div>
    
    <div class="table-container">
        <div id="loading" class="loading">Ready to load data...</div>
        <table id="forecastTable" style="display: none;">
            <thead>
                <tr>
                    <th rowspan="2">Degree</th>
                    <th colspan="2">Ground Truth (M1 @ 3.2GHz)</th>
                    <th rowspan="2">Target CPU Time</th>
                    <th rowspan="2">Status</th>
                </tr>
                <tr>
                    <th>Time</th>
                    <th>Arch-Independent Cycles</th>
                </tr>
            </thead>
            <tbody id="tableBody">
            </tbody>
        </table>
    </div>

    <script>
        // Load available braid lengths on page load
        fetch('/api/braids')
            .then(response => response.json())
            .then(braids => {
                const select = document.getElementById('braidLength');
                select.innerHTML = '';
                braids.forEach(braid => {
                    const option = document.createElement('option');
                    option.value = braid;
                    option.textContent = braid;
                    select.appendChild(option);
                });
                if (braids.length > 0) {
                    select.value = braids[0];
                    updateTable();
                }
            })
            .catch(error => {
                document.getElementById('loading').innerHTML = 
                    '<div class="error">Error loading braid lengths: ' + error.message + '</div>';
            });
        
        function updateTable() {
            const braidLength = document.getElementById('braidLength').value;
            const targetCpuFreq = document.getElementById('targetCpuFreq').value;
            const degreeMin = document.getElementById('degreeMin').value;
            const degreeMax = document.getElementById('degreeMax').value;
            
            if (!braidLength) return;
            
            document.getElementById('loading').style.display = 'block';
            document.getElementById('forecastTable').style.display = 'none';
            document.getElementById('loading').textContent = 'Loading forecast data...';
            
            const params = new URLSearchParams({
                braid_length: braidLength,
                target_cpu_freq: targetCpuFreq,
                degree_min: degreeMin,
                degree_max: degreeMax
            });
            
            fetch('/api/forecast?' + params)
                .then(response => response.json())
                .then(data => {
                    if (data.error) {
                        throw new Error(data.error);
                    }
                    
                    // Update statistics
                    const measuredCount = data.results.filter(r => r.status === 'measured').length;
                    const predictedCount = data.results.length - measuredCount;
                    
                    document.getElementById('stats').innerHTML = 
                        `<strong>Braid Length:</strong> ${braidLength} | 
                         <strong>Target CPU:</strong> ${targetCpuFreq} GHz | 
                         <strong>R²:</strong> ${data.r_squared.toFixed(6)} | 
                         <strong>Measured:</strong> ${measuredCount} | 
                         <strong>Predicted:</strong> ${predictedCount} | 
                         <strong>Total:</strong> ${data.results.length} degrees`;
                    
                    // Populate table
                    const tbody = document.getElementById('tableBody');
                    tbody.innerHTML = '';
                    
                    data.results.forEach(result => {
                        const row = document.createElement('tr');
                        row.className = result.status;
                        
                        row.innerHTML = `
                            <td class="degree-cell">${result.degree}</td>
                            <td class="time-cell">${result.ground_truth_formatted_time}</td>
                            <td class="cycles-cell">${result.ground_truth_formatted_cycles}</td>
                            <td class="time-cell">${result.predicted_formatted_time}</td>
                            <td class="status-${result.status}">${result.status}</td>
                        `;
                        
                        tbody.appendChild(row);
                    });
                    
                    document.getElementById('loading').style.display = 'none';
                    document.getElementById('forecastTable').style.display = 'table';
                })
                .catch(error => {
                    document.getElementById('loading').innerHTML = 
                        '<div class="error">Error loading forecast data: ' + error.message + '</div>';
                });
        }
    </script>
</body>
</html>"""
        
        self.send_response(200)
        self.send_header('Content-type', 'text/html')
        self.end_headers()
        self.wfile.write(html.encode())
    
    def serve_braids_api(self):
        """API endpoint to get available braid lengths"""
        braids = sorted(self.data.keys()) if self.data else []
        
        self.send_response(200)
        self.send_header('Content-type', 'application/json')
        self.send_header('Access-Control-Allow-Origin', '*')
        self.end_headers()
        
        response = json.dumps(braids)
        self.wfile.write(response.encode())
    
    def serve_forecast_api(self):
        """API endpoint to get forecast data"""
        try:
            # Parse query parameters
            parsed_url = urllib.parse.urlparse(self.path)
            params = urllib.parse.parse_qs(parsed_url.query)
            
            braid_length = int(params.get('braid_length', [8])[0])
            target_cpu_freq = float(params.get('target_cpu_freq', [3.2])[0])
            degree_min = int(params.get('degree_min', [5])[0])
            degree_max = int(params.get('degree_max', [50])[0])
            
            if braid_length not in self.data:
                raise ValueError(f"No data for braid length {braid_length}")
            
            # Generate forecast
            degrees, times = zip(*self.data[braid_length])
            degree_range = range(degree_min, degree_max + 1)
            results, r_squared = generate_forecast_table(degrees, times, degree_range, target_cpu_freq, 3.2)
            
            response_data = {
                'results': results,
                'r_squared': r_squared,
                'braid_length': braid_length
            }
            
            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.send_header('Access-Control-Allow-Origin', '*')
            self.end_headers()
            
            response = json.dumps(response_data)
            self.wfile.write(response.encode())
            
        except Exception as e:
            self.send_response(500)
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            
            error_response = json.dumps({'error': str(e)})
            self.wfile.write(error_response.encode())

def start_web_server(data, port=8080):
    """Start the web server with forecast data"""
    handler = lambda *args, **kwargs: ForecastWebHandler(*args, data=data, **kwargs)
    
    with socketserver.TCPServer(("", port), handler) as httpd:
        print(f"🌐 Web server running at http://localhost:{port}")
        print(f"📊 Open the URL above to view forecast tables")
        print("Press Ctrl+C to stop the server")
        
        # Try to open browser automatically
        def open_browser():
            time.sleep(1)  # Give server time to start
            webbrowser.open(f'http://localhost:{port}')
        
        browser_thread = threading.Thread(target=open_browser)
        browser_thread.daemon = True
        browser_thread.start()
        
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\n🛑 Server stopped")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Start web-based forecast table viewer')
    parser.add_argument('--input', '-i', default='profiling_results.txt',
                       help='Input profiling file (default: profiling_results.txt)')
    parser.add_argument('--port', '-p', type=int, default=8080,
                       help='Web server port (default: 8080)')
    
    args = parser.parse_args()
    
    # Load data
    print(f"📁 Loading data from {args.input}...")
    try:
        data = parse_profiling_data(args.input)
        if not data:
            print("❌ No valid data found in the file!")
            return
        
        print(f"✅ Found data for braid lengths: {sorted(data.keys())}")
        
        # Start web server
        start_web_server(data, args.port)
        
    except FileNotFoundError:
        print(f"❌ File not found: {args.input}")
    except Exception as e:
        print(f"❌ Error loading data: {e}")

if __name__ == "__main__":
    main()