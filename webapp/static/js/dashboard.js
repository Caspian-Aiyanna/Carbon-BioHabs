// Dashboard State Management
class Dashboard {
    constructor() {
        this.state = {
            selectedElephant: 'E3',
            selectedRun: 'A',
            selectedLayer: 'hybrid',
            compareMode: false
        };
        this.init();
    }

    init() {
        this.loadElephants();
        this.setupEventListeners();
        this.loadInitialData();
    }

    loadElephants() {
        const elephantList = document.getElementById('elephantList');
        if (!elephantList) return;

        elephantList.innerHTML = CONFIG.ELEPHANTS.map(elephant => `
            <div class="elephant-item ${elephant.id === this.state.selectedElephant ? 'active' : ''}" 
                 data-elephant="${elephant.id}">
                <span class="elephant-name">${elephant.name}</span>
                <span class="elephant-meta">${elephant.type} • Pop: ${elephant.population}</span>
            </div>
        `).join('');
    }

    setupEventListeners() {
        // Elephant selection
        document.getElementById('elephantList')?.addEventListener('click', (e) => {
            const item = e.target.closest('.elephant-item');
            if (item) {
                this.selectElephant(item.dataset.elephant);
            }
        });

        // Run selection
        document.querySelectorAll('.run-btn').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const run = e.target.dataset.run;
                if (run === 'compare') {
                    this.toggleCompareMode();
                } else {
                    this.selectRun(run);
                }
            });
        });

        // Layer selection
        document.getElementById('mapLayerSelect')?.addEventListener('change', (e) => {
            this.selectLayer(e.target.value);
        });

        // Chart tabs
        document.querySelectorAll('.chart-tab').forEach(tab => {
            tab.addEventListener('click', (e) => {
                const chartType = e.target.dataset.chart;
                this.updateChart(chartType);
            });
        });

        // Refresh button
        document.getElementById('refreshBtn')?.addEventListener('click', () => {
            this.refresh();
        });

        // Export button
        document.getElementById('exportBtn')?.addEventListener('click', () => {
            this.exportData();
        });

        // Fullscreen button
        document.getElementById('fullscreenBtn')?.addEventListener('click', () => {
            this.toggleFullscreen();
        });

        // Navigation
        document.querySelectorAll('.nav-link').forEach(link => {
            link.addEventListener('click', (e) => {
                e.preventDefault();
                const targetId = link.getAttribute('href').substring(1);
                this.navigateTo(targetId);
            });
        });
    }

    async loadInitialData() {
        await this.updateDashboard();
    }

    async selectElephant(elephantId) {
        this.state.selectedElephant = elephantId;

        // Update UI
        document.querySelectorAll('.elephant-item').forEach(item => {
            item.classList.toggle('active', item.dataset.elephant === elephantId);
        });

        await this.updateDashboard();
        this.loadFigures();
    }

    async selectRun(run) {
        this.state.selectedRun = run;
        this.state.compareMode = false;

        // Update UI
        document.querySelectorAll('.run-btn').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.run === run);
        });

        await this.updateDashboard();
        this.loadFigures();
    }

    async selectLayer(layer) {
        this.state.selectedLayer = layer;
        await this.updateMap();
    }

    toggleCompareMode() {
        this.state.compareMode = !this.state.compareMode;

        document.querySelectorAll('.run-btn').forEach(btn => {
            if (btn.dataset.run === 'compare') {
                btn.classList.toggle('active', this.state.compareMode);
            } else {
                btn.classList.toggle('active', false);
            }
        });

        this.updateDashboard();
        this.loadFigures();
    }

    async updateDashboard() {
        await Promise.all([
            this.updateStats(),
            this.updateMap(),
            this.updateChart(),
            this.updateTable()
        ]);
    }

    async updateStats() {
        try {
            const summary = await api.getSummary(
                this.state.selectedElephant,
                this.state.selectedRun
            );

            // Update stat cards
            document.getElementById('modelsCompleted').textContent = summary.models_completed || 3;
            document.getElementById('avgSuitability').textContent = (summary.avg_suitability || 0).toFixed(2);
            document.getElementById('avgUncertainty').textContent = (summary.avg_uncertainty || 0).toFixed(2);
            document.getElementById('carbonPotential').textContent = Math.round(summary.total_carbon || 0);

            // Update change indicators
            if (summary.change_suitability) {
                const changeEl = document.getElementById('suitabilityChange');
                changeEl.textContent = `${summary.change_suitability > 0 ? '+' : ''}${summary.change_suitability.toFixed(1)}%`;
                changeEl.classList.toggle('positive', summary.change_suitability > 0);
                changeEl.classList.toggle('negative', summary.change_suitability < 0);
            }
        } catch (error) {
            console.error('Error updating stats:', error);
        }
    }

    async updateMap() {
        if (!mapManager) return;

        await mapManager.loadRaster(
            this.state.selectedElephant,
            this.state.selectedRun,
            this.state.selectedLayer
        );
    }

    async updateChart(chartType = 'correlation') {
        if (!chartManager) return;

        // Update active tab
        document.querySelectorAll('.chart-tab').forEach(tab => {
            tab.classList.toggle('active', tab.dataset.chart === chartType);
        });

        await chartManager.renderChart(
            chartType,
            this.state.selectedElephant,
            'A',
            'B'
        );
    }

    async updateTable() {
        try {
            const comparison = await api.getComparison(this.state.selectedElephant);
            const tbody = document.getElementById('resultsTableBody');
            if (!tbody) return;

            tbody.innerHTML = `
                <tr>
                    <td>${this.state.selectedElephant}</td>
                    <td>${this.state.selectedRun}</td>
                    <td>${(comparison.h2o_vs_ssdm_pearson || 0).toFixed(3)}</td>
                    <td>${(comparison.h2o_vs_ssf_pearson || 0).toFixed(3)}</td>
                    <td>${(comparison.h2o_vs_ssdm_rmse || 0).toFixed(3)}</td>
                    <td>${(comparison.avg_uncertainty || 0).toFixed(3)}</td>
                    <td>${Math.round(comparison.total_carbon || 0)}</td>
                </tr>
            `;
        } catch (error) {
            console.error('Error updating table:', error);
        }
    }

    async refresh() {
        const btn = document.getElementById('refreshBtn');
        if (btn) {
            btn.style.animation = 'spin 1s linear';
            setTimeout(() => {
                btn.style.animation = '';
            }, 1000);
        }

        await this.updateDashboard();
    }

    async exportData() {
        try {
            const blob = await api.exportData(
                this.state.selectedElephant,
                this.state.selectedRun,
                'geotiff'
            );

            // Create download link
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = `${this.state.selectedElephant}_${this.state.selectedRun}_results.zip`;
            document.body.appendChild(a);
            a.click();
            window.URL.revokeObjectURL(url);
            document.body.removeChild(a);
        } catch (error) {
            console.error('Error exporting data:', error);
            alert('Export failed. Please try again.');
        }
    }

    toggleFullscreen() {
        const mapPanel = document.querySelector('.map-panel');
        if (!document.fullscreenElement) {
            mapPanel.requestFullscreen();
        } else {
            document.exitFullscreen();
        }
    }

    navigateTo(sectionId) {
        // Update nav links
        document.querySelectorAll('.nav-link').forEach(link => {
            link.classList.toggle('active', link.getAttribute('href') === `#${sectionId}`);
        });

        // Update sections
        document.querySelectorAll('section').forEach(section => {
            if (section.id === `${sectionId}-section`) {
                section.classList.remove('hidden-section');
                section.classList.add('active-section');
            } else {
                section.classList.add('hidden-section');
                section.classList.remove('active-section');
            }
        });

        // Load content if needed
        if (sectionId === 'pipeline' || sectionId === 'analysis') {
            this.loadFigures();
        } else if (sectionId === 'export') {
            this.updateExportSection();
        }
    }

    updateExportSection() {
        const elephantSpan = document.getElementById('exportElephant');
        const runSpan = document.getElementById('exportRun');

        if (elephantSpan) elephantSpan.textContent = this.state.selectedElephant;
        if (runSpan) runSpan.textContent = this.state.compareMode ? 'Compare' : `Run ${this.state.selectedRun}`;

        // Setup download buttons
        const downloadCurrentBtn = document.getElementById('downloadCurrentBtn');
        if (downloadCurrentBtn) {
            downloadCurrentBtn.onclick = () => this.exportData();
        }
    }

    async loadFigures() {
        const pipelineGallery = document.getElementById('pipelineGallery');
        const analysisGallery = document.getElementById('analysisGallery');

        if (!pipelineGallery || !analysisGallery) return;

        // Only show loading if empty
        if (!pipelineGallery.children.length) pipelineGallery.innerHTML = '<div class="loading-spinner"></div>';
        if (!analysisGallery.children.length) analysisGallery.innerHTML = '<div class="loading-spinner"></div>';

        try {
            const data = await api.getFigures(this.state.selectedElephant);

            this.renderGallery(pipelineGallery, data.pipeline);
            this.renderGallery(analysisGallery, data.analysis);
        } catch (error) {
            console.error('Error loading figures:', error);
            pipelineGallery.innerHTML = '<div class="empty-state">Error loading figures</div>';
            analysisGallery.innerHTML = '<div class="empty-state">Error loading figures</div>';
        }
    }

    renderGallery(container, files) {
        let filteredFiles = files || [];

        // Filter based on mode/run
        if (this.state.compareMode) {
            filteredFiles = filteredFiles.filter(f => f.includes('Compare'));
        } else {
            const run = this.state.selectedRun;
            // Match _A_ or _B_ or E3A/E3B
            filteredFiles = filteredFiles.filter(f =>
                f.includes(`_${run}_`) ||
                f.includes(`_${run}.`) ||
                f.includes(`${this.state.selectedElephant}${run}`) ||
                f.includes(`${this.state.selectedElephant}_${run}`)
            );
        }

        if (filteredFiles.length === 0) {
            container.innerHTML = '<div class="empty-state">No visualizations available for this selection</div>';
            return;
        }

        container.innerHTML = filteredFiles.map(file => {
            const fileName = file.split('/').pop();
            return `
            <div class="gallery-item">
                <div class="gallery-image-container">
                    <img src="/api/results/file/${file}" 
                         alt="${fileName}" 
                         class="gallery-image"
                         onclick="window.open(this.src, '_blank')">
                </div>
                <div class="gallery-caption">
                    <div class="gallery-title">${fileName.replace(/_/g, ' ').replace('.png', '')}</div>
                    <div class="gallery-meta">
                        <span>${this.state.selectedElephant} ${this.state.compareMode ? 'Compare' : this.state.selectedRun}</span>
                        <button class="btn-small" onclick="window.open('/api/results/file/${file}', '_blank')">
                            View Full
                        </button>
                    </div>
                </div>
            </div>
            `;
        }).join('');
    }
}

// Create global dashboard instance
let dashboard = null;

// Initialize dashboard when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
    dashboard = new Dashboard();
});
