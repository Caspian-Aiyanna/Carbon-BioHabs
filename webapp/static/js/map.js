// Map Module for BioHabs Dashboard
class MapManager {
    constructor(containerId) {
        this.containerId = containerId;
        this.map = null;
        this.layers = {};
        this.currentLayer = null;
        this.legend = null;
        this.init();
    }

    init() {
        // Initialize Leaflet map
        this.map = L.map(this.containerId, {
            center: CONFIG.MAP.center,
            zoom: CONFIG.MAP.zoom,
            minZoom: CONFIG.MAP.minZoom,
            maxZoom: CONFIG.MAP.maxZoom,
            zoomControl: true
        });

        // Add base tile layer
        L.tileLayer(CONFIG.MAP.tileLayer, {
            attribution: CONFIG.MAP.attribution,
            maxZoom: CONFIG.MAP.maxZoom
        }).addTo(this.map);

        // Add scale control
        L.control.scale({ position: 'bottomleft' }).addTo(this.map);
    }

    async loadRaster(elephant, run, layerType) {
        try {
            // Show loading state
            this.showLoading(true);

            // Fetch GeoJSON data from API
            const data = await api.getRaster(elephant, run, layerType);

            // Remove existing layer
            if (this.currentLayer) {
                this.map.removeLayer(this.currentLayer);
            }

            // Get color scheme
            const colorScheme = CONFIG.LAYERS[layerType].colorScheme;
            const colors = CONFIG.COLORS[colorScheme];

            // Create GeoJSON layer with styling
            this.currentLayer = L.geoJSON(data, {
                style: (feature) => this.getFeatureStyle(feature, colors),
                onEachFeature: (feature, layer) => {
                    if (feature.properties) {
                        layer.bindPopup(this.createPopup(feature.properties));
                    }
                }
            }).addTo(this.map);

            // Fit bounds to layer
            this.map.fitBounds(this.currentLayer.getBounds());

            // Update legend
            this.updateLegend(layerType, colors, data);

            this.showLoading(false);
        } catch (error) {
            console.error('Error loading raster:', error);
            this.showError('Failed to load map layer');
            this.showLoading(false);
        }
    }

    getFeatureStyle(feature, colors) {
        const value = feature.properties.value;
        const color = this.getColor(value, colors);

        return {
            fillColor: color,
            weight: 0,
            opacity: 1,
            fillOpacity: 0.7
        };
    }

    getColor(value, colors) {
        // Normalize value to 0-1 range
        const normalized = Math.max(0, Math.min(1, value));
        const index = Math.floor(normalized * (colors.length - 1));
        return colors[index];
    }

    createPopup(properties) {
        return `
            <div class="map-popup">
                <strong>Value:</strong> ${properties.value.toFixed(3)}<br>
                ${properties.uncertainty ? `<strong>Uncertainty:</strong> ${properties.uncertainty.toFixed(3)}<br>` : ''}
                ${properties.carbon ? `<strong>Carbon:</strong> ${properties.carbon.toFixed(1)} MgC<br>` : ''}
            </div>
        `;
    }

    updateLegend(layerType, colors, data) {
        const legendEl = document.getElementById('mapLegend');
        if (!legendEl) return;

        // Calculate value range
        const values = data.features.map(f => f.properties.value);
        const min = Math.min(...values);
        const max = Math.max(...values);

        // Create legend HTML
        const steps = 5;
        let html = `<div class="legend-title">${CONFIG.LAYERS[layerType].name}</div>`;
        html += '<div class="legend-scale">';

        for (let i = 0; i < steps; i++) {
            const value = min + (max - min) * (i / (steps - 1));
            const color = this.getColor((value - min) / (max - min), colors);
            html += `
                <div class="legend-item">
                    <div class="legend-color" style="background: ${color}"></div>
                    <span>${value.toFixed(2)}</span>
                </div>
            `;
        }

        html += '</div>';
        legendEl.innerHTML = html;
    }

    showLoading(show) {
        const overlay = document.getElementById('loadingOverlay');
        if (overlay) {
            overlay.classList.toggle('active', show);
        }
    }

    showError(message) {
        // TODO: Implement error notification
        console.error(message);
    }

    clearLayers() {
        if (this.currentLayer) {
            this.map.removeLayer(this.currentLayer);
            this.currentLayer = null;
        }
    }

    resize() {
        if (this.map) {
            this.map.invalidateSize();
        }
    }
}

// Create global map instance
let mapManager = null;

// Initialize map when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
    mapManager = new MapManager('map');
});
