"""
BioHabs Dashboard - Flask Backend
Main application server for running pipelines and serving results
"""

from flask import Flask, jsonify, request, send_file, send_from_directory
from flask_cors import CORS
import os
import json
from pathlib import Path
import subprocess
import threading
from datetime import datetime

app = Flask(__name__, static_folder='../static', static_url_path='/static')
CORS(app)  # Enable CORS for frontend

# Configuration
BASE_DIR = Path(__file__).parent.parent.parent
RESULTS_DIR = BASE_DIR / 'results'
DATA_DIR = BASE_DIR / 'data'
BIOHABS_DIR = BASE_DIR / 'BioHabs'

# Pipeline status tracking
pipeline_status = {
    'h2o': 'idle',
    'ssdm': 'idle',
    'issa': 'idle'
}

pipeline_jobs = {}

# ===== UTILITY FUNCTIONS =====

def load_csv_data(filepath):
    """Load CSV file and return as list of dicts"""
    import csv
    if not os.path.exists(filepath):
        return []
    
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        return list(reader)

def load_geotiff_as_geojson(filepath):
    """Convert GeoTIFF to GeoJSON for web display"""
    try:
        import rasterio
        from rasterio.features import shapes
        import numpy as np
        
        with rasterio.open(filepath) as src:
            image = src.read(1)
            mask = image != src.nodata
            
            # Sample for performance (every 10th pixel)
            image_sampled = image[::10, ::10]
            mask_sampled = mask[::10, ::10]
            
            features = []
            for geom, value in shapes(image_sampled, mask=mask_sampled, transform=src.transform):
                # Handle NaN values for JSON serialization
                val = float(value)
                if np.isnan(val):
                    val = None
                    
                features.append({
                    'type': 'Feature',
                    'geometry': geom,
                    'properties': {'value': val}
                })
            
            return {
                'type': 'FeatureCollection',
                'features': features
            }
    except Exception as e:
        print(f"Error loading GeoTIFF: {e}")
        return {'type': 'FeatureCollection', 'features': []}

# ===== API ENDPOINTS =====

@app.route('/api/elephants', methods=['GET'])
def get_elephants():
    """Get list of available elephants"""
    elephants = [
        {'id': 'E1', 'name': 'E1', 'type': 'Bull', 'population': 1},
        {'id': 'E2', 'name': 'E2', 'type': 'Bull', 'population': 1},
        {'id': 'E3', 'name': 'E3', 'type': 'Herd', 'population': 30},
        {'id': 'E4', 'name': 'E4', 'type': 'Herd', 'population': 20},
        {'id': 'E5', 'name': 'E5', 'type': 'Herd', 'population': 20},
        {'id': 'E6', 'name': 'E6', 'type': 'Bull', 'population': 1}
    ]
    return jsonify(elephants)

@app.route('/api/results/<elephant>/<run>', methods=['GET'])
def get_results(elephant, run):
    """Get results for specific elephant and run"""
    results_path = RESULTS_DIR / 'uncertainty' / elephant / 'combined'
    
    if not results_path.exists():
        return jsonify({'error': 'Results not found'}), 404
    
    files = {
        'hybrid': f'{elephant}_{run}_Hybrid_Ensemble.tif',
        'uncertainty': f'{elephant}_{run}_Uncertainty_Total.tif',
        'carbon': f'{elephant}_{run}_CSP.tif'
    }
    
    available_files = {k: str(results_path / v) for k, v in files.items() 
                      if (results_path / v).exists()}
    
    return jsonify({
        'elephant': elephant,
        'run': run,
        'files': available_files
    })

@app.route('/api/raster/<elephant>/<run>/<layer>', methods=['GET'])
def get_raster(elephant, run, layer):
    """Get raster data as GeoJSON"""
    # Map layer names to file patterns
    # Note: Comparison rasters use E3B format, while Uncertainty rasters use E3_B format
    layer_files = {
        'hybrid': f'{elephant}_{run}_Hybrid_Ensemble.tif',
        'h2o': f'{elephant}{run}_H2O_mean.tif',
        'ssdm': f'{elephant}{run}_SSDM_mean.tif',
        'issa': f'{elephant}{run}_SSF.tif',
        'uncertainty': f'{elephant}_{run}_Uncertainty_Total.tif',
        'carbon': f'{elephant}_{run}_CSP.tif'
    }
    
    # Try different result directories
    search_paths = [
        RESULTS_DIR / 'uncertainty' / elephant / 'combined',
        RESULTS_DIR / 'compare' / run / '01_between_methods' / 'rasters'
    ]
    
    filepath = None
    for search_path in search_paths:
        potential_file = search_path / layer_files.get(layer, '')
        if potential_file.exists():
            filepath = potential_file
            break
    
    if not filepath:
        return jsonify({'error': 'Raster file not found'}), 404
    
    geojson = load_geotiff_as_geojson(str(filepath))
    return jsonify(geojson)

@app.route('/api/summary/<elephant>/<run>', methods=['GET'])
def get_summary(elephant, run):
    """Get summary statistics"""
    # Load comparison data
    comparison_file = RESULTS_DIR / 'compare' / run / '02_tables' / 'per_species_summary_long.csv'
    
    if not comparison_file.exists():
        return jsonify({
            'avg_suitability': 0.75,
            'avg_uncertainty': 0.15,
            'total_carbon': 1000,
            'models_completed': 3
        })
    
    data = load_csv_data(str(comparison_file))
    elephant_data = next((d for d in data if d.get('dataset') == f'{elephant}{run}'), {})
    
    return jsonify({
        'avg_suitability': float(elephant_data.get('h2o_vs_ssdm_pearson_mean', 0.75)),
        'avg_uncertainty': float(elephant_data.get('h2o_vs_ssdm_rmse_mean', 0.15)),
        'total_carbon': 1000,  # TODO: Calculate from CSP raster
        'models_completed': 3,
        'change_suitability': 5.2  # TODO: Calculate A vs B change
    })

@app.route('/api/comparison/<elephant>', methods=['GET'])
def get_comparison(elephant):
    """Get comparison metrics between models"""
    # Try both runs
    for run in ['A', 'B']:
        comparison_file = RESULTS_DIR / 'compare' / run / '02_tables' / 'per_species_summary_long.csv'
        if comparison_file.exists():
            data = load_csv_data(str(comparison_file))
            elephant_data = next((d for d in data if d.get('dataset') == f'{elephant}{run}'), {})
            
            if elephant_data:
                return jsonify({
                    'h2o_vs_ssdm_pearson': float(elephant_data.get('h2o_vs_ssdm_pearson_mean', 0)),
                    'h2o_vs_ssdm_rmse': float(elephant_data.get('h2o_vs_ssdm_rmse_mean', 0)),
                    'h2o_vs_ssdm_mae': float(elephant_data.get('h2o_vs_ssdm_mae_mean', 0)),
                    'h2o_vs_ssdm_jaccard': float(elephant_data.get('h2o_vs_ssdm_jaccard_mean', 0)),
                    'h2o_vs_ssf_pearson': float(elephant_data.get('h2o_vs_ssf_pearson_mean', 0)),
                    'h2o_vs_ssf_rmse': float(elephant_data.get('h2o_vs_ssf_rmse_mean', 0)),
                    'ssdm_vs_ssf_pearson': float(elephant_data.get('ssdm_vs_ssf_pearson_mean', 0)),
                    'avg_uncertainty': 0.15,
                    'total_carbon': 1000
                })
    
    return jsonify({'error': 'No comparison data found'}), 404

@app.route('/api/pipeline/status', methods=['GET'])
def get_pipeline_status():
    """Get current pipeline status"""
    return jsonify(pipeline_status)

@app.route('/api/pipeline/run', methods=['POST'])
def run_pipeline():
    """Run BioHabs pipeline"""
    config = request.json
    elephant = config.get('elephant', 'E3')
    run = config.get('run', 'B')
    mode = config.get('mode', 'FAST')
    
    job_id = f"{elephant}_{run}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    
    def run_script():
        try:
            # Update status
            pipeline_status['h2o'] = 'running'
            
            # Run H2O
            subprocess.run([
                'Rscript', str(BIOHABS_DIR / '03a_h2o_replicates.R'),
                '--run', run, '--mode', mode, '--species', elephant
            ], check=True)
            
            pipeline_status['h2o'] = 'complete'
            pipeline_status['ssdm'] = 'running'
            
            # Run SSDM
            subprocess.run([
                'Rscript', str(BIOHABS_DIR / '04a_ssdm_replicates.R'),
                '--run', run, '--mode', mode, '--species', elephant
            ], check=True)
            
            pipeline_status['ssdm'] = 'complete'
            pipeline_status['issa'] = 'running'
            
            # Run iSSA
            subprocess.run([
                'Rscript', str(BIOHABS_DIR / '01_a_iSSA.R'),
                '--run', run, '--mode', mode, '--species', elephant
            ], check=True)
            
            pipeline_status['issa'] = 'complete'
            
        except Exception as e:
            print(f"Pipeline error: {e}")
            for key in pipeline_status:
                if pipeline_status[key] == 'running':
                    pipeline_status[key] = 'error'
    
    # Run in background thread
    thread = threading.Thread(target=run_script)
    thread.start()
    pipeline_jobs[job_id] = thread
    
    return jsonify({'job_id': job_id, 'status': 'started'})

@app.route('/api/export/<elephant>/<run>', methods=['GET'])
def export_data(elephant, run):
    """Export results as zip file"""
    import zipfile
    import io
    
    # Create zip file in memory
    zip_buffer = io.BytesIO()
    
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        results_path = RESULTS_DIR / 'uncertainty' / elephant / 'combined'
        
        if results_path.exists():
            for file in results_path.glob(f'{elephant}_{run}_*.tif'):
                zip_file.write(file, file.name)
    
    zip_buffer.seek(0)
    
    return send_file(
        zip_buffer,
        mimetype='application/zip',
        as_attachment=True,
        download_name=f'{elephant}_{run}_results.zip'
    )

# ===== STATIC FILE SERVING =====

# Get absolute path to webapp directory
WEBAPP_DIR = Path(__file__).parent.parent

@app.route('/')
def serve_landing():
    """Serve the landing page"""
    return send_file(WEBAPP_DIR / 'index.html')

@app.route('/frontend/<path:path>')
def serve_frontend(path):
    """Serve frontend files"""
    return send_from_directory(WEBAPP_DIR / 'frontend', path)

@app.route('/api/figures/<elephant>', methods=['GET'])
def get_figures(elephant):
    """Get list of available figures for an elephant"""
    pipeline_files = []
    analysis_files = []
    
    # 1. Uncertainty Figures
    uncertainty_dir = RESULTS_DIR / 'uncertainty' / elephant / 'figures'
    if uncertainty_dir.exists():
        for f in uncertainty_dir.glob('*.png'):
            rel_path = f.relative_to(RESULTS_DIR).as_posix()
            if 'Scatter' in f.name or 'Violin' in f.name or 'metrics' in f.name:
                analysis_files.append(rel_path)
            else:
                pipeline_files.append(rel_path)
                
    # 2. Comparison Figures (Method comparison within runs)
    # Structure: results/compare/<run>/01_between_methods/plots
    compare_dir = RESULTS_DIR / 'compare'
    if compare_dir.exists():
        for run_dir in compare_dir.iterdir():
            if not run_dir.is_dir():
                continue
            plots_dir = run_dir / '01_between_methods' / 'plots'
            if plots_dir.exists():
                # Filter for files related to this elephant (e.g. E3A, E3B)
                for f in plots_dir.glob(f'*{elephant}*.png'):
                    rel_path = f.relative_to(RESULTS_DIR).as_posix()
                    # "diff_mean" are maps (Pipeline), others might be plots
                    if 'diff_mean' in f.name:
                        pipeline_files.append(rel_path)
                    else:
                        analysis_files.append(rel_path)
            
    return jsonify({
        'pipeline': sorted(pipeline_files),
        'analysis': sorted(analysis_files)
    })

@app.route('/api/results/file/<path:filename>')
def serve_result_file(filename):
    """Serve files from the results directory"""
    return send_from_directory(RESULTS_DIR, filename)

# ===== MAIN =====

if __name__ == '__main__':
    print("=" * 50)
    print("Starting BioHabs Dashboard Server...")
    print(f"Base directory: {BASE_DIR}")
    print(f"Results directory: {RESULTS_DIR}")
    print("Server running at: http://localhost:5000")
    print("=" * 50)
    
    app.run(debug=True, host='0.0.0.0', port=5000)
