# 🚀 BioHabs Dashboard - Quick Start Guide

## Installation

### 1. Install Python Dependencies

```bash
cd webapp
pip install -r requirements.txt
```

### 2. Start the Backend Server

```bash
python backend/app.py
```

The server will start at `http://localhost:5000`

### 3. Open the Dashboard

Open `frontend/index.html` in your browser, or use a local server:

```bash
# Option 1: Python HTTP Server
cd frontend
python -m http.server 8080

# Option 2: Node.js HTTP Server (if you have Node.js)
npx http-server frontend -p 8080
```

Then navigate to: `http://localhost:8080`

## Usage

### View Results

1. **Select Elephant**: Click on E1-E6 in the sidebar
2. **Select Run**: Choose Run A (Before) or Run B (After)
3. **View Map**: Interactive map shows habitat suitability
4. **View Charts**: Model comparison and metrics
5. **Export Data**: Download results as GeoTIFF/CSV

### Run Pipeline from Browser

1. Click "Pipeline" in navigation
2. Select elephant and run
3. Click "Start Pipeline"
4. Monitor progress in real-time

## Features

✅ **Interactive Maps** - Leaflet.js with GeoTIFF rendering  
✅ **Real-time Charts** - Chart.js visualizations  
✅ **Pipeline Control** - Run models from browser  
✅ **Data Export** - Download results  
✅ **Responsive Design** - Works on all devices  

## Troubleshooting

### Backend not starting?
- Check Python version (3.8+)
- Install dependencies: `pip install -r requirements.txt`
- Check port 5000 is not in use

### Frontend not loading data?
- Ensure backend is running
- Check browser console for errors
- Verify CORS is enabled

### Maps not displaying?
- Check GeoTIFF files exist in `results/` directory
- Verify rasterio is installed
- Check browser console for errors

## API Endpoints

- `GET /api/elephants` - List elephants
- `GET /api/results/{elephant}/{run}` - Get results
- `GET /api/raster/{elephant}/{run}/{layer}` - Get map data
- `GET /api/summary/{elephant}/{run}` - Get statistics
- `POST /api/pipeline/run` - Run pipeline
- `GET /api/export/{elephant}/{run}` - Export data

## Development

### File Structure
```
webapp/
├── frontend/
│   └── index.html          # Main dashboard
├── backend/
│   └── app.py              # Flask server
├── static/
│   ├── css/                # Styles
│   └── js/                 # JavaScript modules
└── requirements.txt        # Python deps
```

### Adding New Features

1. **New API Endpoint**: Add to `backend/app.py`
2. **New Visualization**: Add to `static/js/charts.js`
3. **New Map Layer**: Update `static/js/map.js`

## Production Deployment

### Using Gunicorn (Recommended)

```bash
pip install gunicorn
gunicorn -w 4 -b 0.0.0.0:5000 backend.app:app
```

### Using Nginx (Reverse Proxy)

```nginx
server {
    listen 80;
    server_name your-domain.com;

    location / {
        proxy_pass http://localhost:5000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }

    location /static {
        alias /path/to/webapp/static;
    }
}
```

## Support

For issues or questions:
1. Check the troubleshooting section
2. Review browser console logs
3. Check Flask server logs

---

**Enjoy your BioHabs Dashboard!** 🎉
