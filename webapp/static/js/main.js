// Main Application Entry Point
(function () {
    'use strict';

    // Application initialization
    function initApp() {
        console.log('🚀 BioHabs Dashboard initializing...');

        // Check if all required libraries are loaded
        if (typeof L === 'undefined') {
            console.error('Leaflet not loaded');
            return;
        }

        if (typeof Chart === 'undefined') {
            console.error('Chart.js not loaded');
            return;
        }

        // Initialize pipeline status checker
        startPipelineMonitor();

        // Setup global error handler
        setupErrorHandler();

        console.log('✅ BioHabs Dashboard ready');
    }

    // Monitor pipeline status
    function startPipelineMonitor() {
        setInterval(async () => {
            try {
                const status = await api.getPipelineStatus();
                updatePipelineStatus(status);
            } catch (error) {
                // Silently fail - backend might not be running
            }
        }, CONFIG.PIPELINE_CHECK_INTERVAL);
    }

    // Update pipeline status indicators
    function updatePipelineStatus(status) {
        const statusItems = document.querySelectorAll('.status-item');

        ['h2o', 'ssdm', 'issa'].forEach((model, index) => {
            const badge = statusItems[index]?.querySelector('.status-badge');
            if (!badge) return;

            const modelStatus = status[model] || 'idle';
            badge.className = `status-badge ${modelStatus}`;
            badge.textContent = modelStatus.charAt(0).toUpperCase() + modelStatus.slice(1);
        });
    }

    // Global error handler
    function setupErrorHandler() {
        window.addEventListener('error', (event) => {
            console.error('Global error:', event.error);
        });

        window.addEventListener('unhandledrejection', (event) => {
            console.error('Unhandled promise rejection:', event.reason);
        });
    }

    // Utility: Show notification
    window.showNotification = function (message, type = 'info') {
        // TODO: Implement toast notifications
        console.log(`[${type.toUpperCase()}] ${message}`);
    };

    // Utility: Format number
    window.formatNumber = function (num, decimals = 2) {
        if (typeof num !== 'number') return 'N/A';
        return num.toFixed(decimals);
    };

    // Utility: Format large numbers
    window.formatLargeNumber = function (num) {
        if (typeof num !== 'number') return 'N/A';
        if (num >= 1000000) return (num / 1000000).toFixed(1) + 'M';
        if (num >= 1000) return (num / 1000).toFixed(1) + 'K';
        return num.toFixed(0);
    };

    // Initialize when DOM is ready
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', initApp);
    } else {
        initApp();
    }
})();
