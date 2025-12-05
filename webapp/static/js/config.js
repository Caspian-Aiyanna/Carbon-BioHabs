// Configuration for BioHabs Dashboard
const CONFIG = {
    // API Configuration
    API_BASE_URL: 'http://localhost:5000/api',
    API_TIMEOUT: 30000,
    
    // Paths
    RESULTS_PATH: '../results',
    DATA_PATH: '../data',
    
    // Map Configuration
    MAP: {
        center: [-1.5, 37.5], // Kenya coordinates
        zoom: 8,
        minZoom: 6,
        maxZoom: 15,
        tileLayer: 'https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png',
        attribution: '© OpenStreetMap contributors'
    },
    
    // Color Schemes
    COLORS: {
        suitability: ['#440154', '#31688e', '#35b779', '#fde724'],
        uncertainty: ['#0d0887', '#7e03a8', '#cc4778', '#f89540', '#f0f921'],
        carbon: ['#081d58', '#253494', '#225ea8', '#1d91c0', '#41b6c4', '#7fcdbb', '#c7e9b4', '#edf8b1', '#ffffd9'],
        diverging: ['#0c97f3', '#ffffff', '#e41b1e']
    },
    
    // Elephants Configuration
    ELEPHANTS: [
        { id: 'E1', name: 'E1', type: 'Bull', population: 1 },
        { id: 'E2', name: 'E2', type: 'Bull', population: 1 },
        { id: 'E3', name: 'E3', type: 'Herd', population: 30 },
        { id: 'E4', name: 'E4', type: 'Herd', population: 20 },
        { id: 'E5', name: 'E5', type: 'Herd', population: 20 },
        { id: 'E6', name: 'E6', type: 'Bull', population: 1 }
    ],
    
    // Runs Configuration
    RUNS: ['A', 'B'],
    
    // Layer Types
    LAYERS: {
        hybrid: { name: 'Hybrid Ensemble', colorScheme: 'suitability' },
        h2o: { name: 'H2O Suitability', colorScheme: 'suitability' },
        ssdm: { name: 'SSDM Suitability', colorScheme: 'suitability' },
        issa: { name: 'iSSA Suitability', colorScheme: 'suitability' },
        uncertainty: { name: 'Uncertainty', colorScheme: 'uncertainty' },
        carbon: { name: 'Carbon Potential', colorScheme: 'carbon' }
    },
    
    // Refresh Intervals (ms)
    REFRESH_INTERVAL: 5000,
    PIPELINE_CHECK_INTERVAL: 2000,
    
    // Chart Configuration
    CHART: {
        defaultType: 'correlation',
        animationDuration: 750,
        responsive: true,
        maintainAspectRatio: false
    }
};

// Export for ES6 modules
if (typeof module !== 'undefined' && module.exports) {
    module.exports = CONFIG;
}
