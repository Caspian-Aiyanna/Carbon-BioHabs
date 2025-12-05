# BioHabs Interactive Dashboard

A modern web application for running and visualizing the BioHabs pipeline results.

## 🎯 Features

- **Interactive Maps**: Visualize habitat suitability, uncertainty, and carbon sequestration
- **Pipeline Control**: Run H2O, SSDM, iSSA models from the browser
- **Real-time Monitoring**: Track pipeline progress with live updates
- **Comparative Analysis**: Side-by-side comparison of Run A vs Run B
- **Data Export**: Download results as GeoTIFF, CSV, or PNG

## 📁 Structure

```
webapp/
├── frontend/           # HTML pages
├── backend/            # Flask API
├── static/
│   ├── css/           # Stylesheets
│   ├── js/            # JavaScript modules
│   └── images/        # Assets
├── templates/         # Jinja2 templates
├── data/              # Cached results
└── docs/              # API documentation
```

## 🚀 Quick Start

### 1. Install Dependencies

```bash
pip install -r requirements.txt
```

### 2. Start the Server

```bash
python webapp/backend/app.py
```

### 3. Open Browser

Navigate to: `http://localhost:5000`

## 🛠️ Technology Stack

- **Frontend**: HTML5, CSS3, JavaScript (ES6+)
- **Maps**: Leaflet.js with GeoTIFF support
- **Charts**: Chart.js for metrics visualization
- **Backend**: Flask (Python)
- **Data**: GeoJSON, GeoTIFF, CSV

## 📊 Available Visualizations

1. **Habitat Suitability Maps** (H2O, SSDM, iSSA)
2. **Uncertainty Maps** (Total, Within-Model, Between-Model)
3. **Carbon Sequestration Potential**
4. **Method Comparison Plots**
5. **Time Series Analysis** (Run A vs Run B)

## 🔧 Configuration

Edit `webapp/backend/config.py` to customize:
- Pipeline paths
- Server port
- Data refresh intervals
- Map tile providers

## 📖 API Documentation

See `webapp/docs/API.md` for full API reference.
