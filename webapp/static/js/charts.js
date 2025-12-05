// Charts Module for BioHabs Dashboard
class ChartManager {
    constructor(canvasId) {
        this.canvasId = canvasId;
        this.chart = null;
        this.currentType = 'correlation';
    }

    async renderChart(type, elephant, runA = 'A', runB = 'B') {
        this.currentType = type;

        try {
            // Fetch data based on chart type
            let data;
            switch (type) {
                case 'correlation':
                    data = await this.getCorrelationData(elephant);
                    break;
                case 'metrics':
                    data = await this.getMetricsData(elephant);
                    break;
                case 'timeseries':
                    data = await this.getTimeSeriesData(elephant, runA, runB);
                    break;
            }

            // Destroy existing chart
            if (this.chart) {
                this.chart.destroy();
            }

            // Create new chart
            const ctx = document.getElementById(this.canvasId).getContext('2d');
            this.chart = new Chart(ctx, this.getChartConfig(type, data));
        } catch (error) {
            console.error('Error rendering chart:', error);
        }
    }

    async getCorrelationData(elephant) {
        const comparison = await api.getComparison(elephant);
        return {
            labels: ['H2O vs SSDM', 'H2O vs iSSA', 'SSDM vs iSSA'],
            datasets: [{
                label: 'Pearson Correlation',
                data: [
                    comparison.h2o_vs_ssdm_pearson,
                    comparison.h2o_vs_ssf_pearson,
                    comparison.ssdm_vs_ssf_pearson
                ],
                backgroundColor: [
                    'rgba(102, 126, 234, 0.8)',
                    'rgba(240, 147, 251, 0.8)',
                    'rgba(67, 233, 123, 0.8)'
                ],
                borderColor: [
                    'rgb(102, 126, 234)',
                    'rgb(240, 147, 251)',
                    'rgb(67, 233, 123)'
                ],
                borderWidth: 2
            }]
        };
    }

    async getMetricsData(elephant) {
        const comparison = await api.getComparison(elephant);
        return {
            labels: ['RMSE', 'MAE', 'Jaccard'],
            datasets: [
                {
                    label: 'H2O vs SSDM',
                    data: [
                        comparison.h2o_vs_ssdm_rmse,
                        comparison.h2o_vs_ssdm_mae,
                        comparison.h2o_vs_ssdm_jaccard
                    ],
                    backgroundColor: 'rgba(102, 126, 234, 0.6)',
                    borderColor: 'rgb(102, 126, 234)',
                    borderWidth: 2
                },
                {
                    label: 'H2O vs iSSA',
                    data: [
                        comparison.h2o_vs_ssf_rmse,
                        comparison.h2o_vs_ssf_mae,
                        comparison.h2o_vs_ssf_jaccard
                    ],
                    backgroundColor: 'rgba(240, 147, 251, 0.6)',
                    borderColor: 'rgb(240, 147, 251)',
                    borderWidth: 2
                }
            ]
        };
    }

    async getTimeSeriesData(elephant, runA, runB) {
        const dataA = await api.getSummary(elephant, runA);
        const dataB = await api.getSummary(elephant, runB);

        return {
            labels: ['Suitability', 'Uncertainty', 'Carbon (MgC)'],
            datasets: [
                {
                    label: `Run ${runA}`,
                    data: [dataA.avg_suitability, dataA.avg_uncertainty, dataA.total_carbon],
                    backgroundColor: 'rgba(230, 159, 0, 0.6)',
                    borderColor: 'rgb(230, 159, 0)',
                    borderWidth: 2
                },
                {
                    label: `Run ${runB}`,
                    data: [dataB.avg_suitability, dataB.avg_uncertainty, dataB.total_carbon],
                    backgroundColor: 'rgba(86, 180, 233, 0.6)',
                    borderColor: 'rgb(86, 180, 233)',
                    borderWidth: 2
                }
            ]
        };
    }

    getChartConfig(type, data) {
        const baseConfig = {
            data: data,
            options: {
                responsive: true,
                maintainAspectRatio: false,
                animation: {
                    duration: CONFIG.CHART.animationDuration
                },
                plugins: {
                    legend: {
                        display: true,
                        position: 'top',
                        labels: {
                            color: 'rgba(255, 255, 255, 0.7)',
                            font: {
                                family: 'Inter',
                                size: 12
                            }
                        }
                    },
                    tooltip: {
                        backgroundColor: 'rgba(15, 15, 35, 0.9)',
                        titleColor: 'rgba(255, 255, 255, 0.9)',
                        bodyColor: 'rgba(255, 255, 255, 0.7)',
                        borderColor: 'rgba(102, 126, 234, 0.5)',
                        borderWidth: 1,
                        padding: 12,
                        displayColors: true
                    }
                },
                scales: {
                    y: {
                        beginAtZero: true,
                        grid: {
                            color: 'rgba(255, 255, 255, 0.1)'
                        },
                        ticks: {
                            color: 'rgba(255, 255, 255, 0.7)'
                        }
                    },
                    x: {
                        grid: {
                            color: 'rgba(255, 255, 255, 0.1)'
                        },
                        ticks: {
                            color: 'rgba(255, 255, 255, 0.7)'
                        }
                    }
                }
            }
        };

        // Chart type specific configurations
        switch (type) {
            case 'correlation':
                return {
                    type: 'bar',
                    ...baseConfig,
                    options: {
                        ...baseConfig.options,
                        scales: {
                            ...baseConfig.options.scales,
                            y: {
                                ...baseConfig.options.scales.y,
                                max: 1,
                                title: {
                                    display: true,
                                    text: 'Correlation Coefficient',
                                    color: 'rgba(255, 255, 255, 0.7)'
                                }
                            }
                        }
                    }
                };
            case 'metrics':
                return {
                    type: 'radar',
                    ...baseConfig,
                    options: {
                        ...baseConfig.options,
                        scales: {
                            r: {
                                beginAtZero: true,
                                grid: {
                                    color: 'rgba(255, 255, 255, 0.1)'
                                },
                                ticks: {
                                    color: 'rgba(255, 255, 255, 0.7)',
                                    backdropColor: 'transparent'
                                },
                                pointLabels: {
                                    color: 'rgba(255, 255, 255, 0.7)'
                                }
                            }
                        }
                    }
                };
            case 'timeseries':
                return {
                    type: 'bar',
                    ...baseConfig,
                    options: {
                        ...baseConfig.options,
                        scales: baseConfig.options.scales
                    }
                };
            default:
                return { type: 'bar', ...baseConfig };
        }
    }

    destroy() {
        if (this.chart) {
            this.chart.destroy();
            this.chart = null;
        }
    }
}

// Create global chart instance
let chartManager = null;

// Initialize chart when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
    chartManager = new ChartManager('comparisonChart');
});
