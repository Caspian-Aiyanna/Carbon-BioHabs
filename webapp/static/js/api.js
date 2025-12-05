// API Client for BioHabs Dashboard
class API {
    constructor(baseURL = CONFIG.API_BASE_URL) {
        this.baseURL = baseURL;
        this.timeout = CONFIG.API_TIMEOUT;
    }

    async request(endpoint, options = {}) {
        const url = `${this.baseURL}${endpoint}`;
        const controller = new AbortController();
        const timeoutId = setTimeout(() => controller.abort(), this.timeout);

        try {
            const response = await fetch(url, {
                ...options,
                signal: controller.signal,
                headers: {
                    'Content-Type': 'application/json',
                    ...options.headers
                }
            });

            clearTimeout(timeoutId);

            if (!response.ok) {
                throw new Error(`HTTP ${response.status}: ${response.statusText}`);
            }

            return await response.json();
        } catch (error) {
            if (error.name === 'AbortError') {
                throw new Error('Request timeout');
            }
            throw error;
        }
    }

    // Get list of available elephants
    async getElephants() {
        return this.request('/elephants');
    }

    // Get results for specific elephant and run
    async getResults(elephant, run) {
        return this.request(`/results/${elephant}/${run}`);
    }

    // Get comparison data
    async getComparison(elephant) {
        return this.request(`/comparison/${elephant}`);
    }

    // Get raster data (GeoJSON)
    async getRaster(elephant, run, layer) {
        return this.request(`/raster/${elephant}/${run}/${layer}`);
    }

    // Get summary statistics
    async getSummary(elephant, run) {
        return this.request(`/summary/${elephant}/${run}`);
    }

    // Get pipeline status
    async getPipelineStatus() {
        return this.request('/pipeline/status');
    }

    // Run pipeline
    async runPipeline(config) {
        return this.request('/pipeline/run', {
            method: 'POST',
            body: JSON.stringify(config)
        });
    }

    // Stop pipeline
    async stopPipeline(jobId) {
        return this.request(`/pipeline/stop/${jobId}`, {
            method: 'POST'
        });
    }

    // Get comparison metrics
    async getMetrics(run) {
        return this.request(`/metrics/${run}`);
    }

    // Get available figures
    async getFigures(elephant) {
        return this.request(`/figures/${elephant}`);
    }

    // Export data
    async exportData(elephant, run, format = 'geotiff') {
        const response = await fetch(`${this.baseURL}/export/${elephant}/${run}?format=${format}`);
        if (!response.ok) {
            throw new Error(`Export failed: ${response.statusText}`);
        }
        return response.blob();
    }
}

// Create global API instance
const api = new API();
