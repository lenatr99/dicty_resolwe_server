import { ReactElement, useEffect, useState, useRef, useCallback } from 'react';
import { connect, ConnectedProps } from 'react-redux';
import {
    Box,
    Typography,
} from '@mui/material';
import { ChartHandle } from '../../common/chart/chart';
import UmapScatterPlot from './umapScatterPlot';
import {
    UmapVisualizationContainer,
    UmapVisualizationControls,
    UmapPlotContainer,
} from './umapVisualization.styles';
import { RootState } from 'redux/rootReducer';
import useSize from 'components/genexpress/common/useSize';
import { getStorage } from 'api';
import { UmapPoint, UmapStorageData, SingleCellStorageData, GeneExpressionsStorageData } from 'types/application';
import { getSelectedGenes } from 'redux/stores/genes';

const mapStateToProps = (state: RootState) => ({
    selectedGenes: getSelectedGenes(state.genes),
});

const connector = connect(mapStateToProps, {
    // Add any required actions if needed
});

type PropsFromRedux = ConnectedProps<typeof connector>;

interface UmapVisualizationProps extends PropsFromRedux {
    storageId?: number;
    geneExpressionsStorageId?: number;
    highlightedCellIds?: string[];
}

const UmapVisualization = ({ storageId, geneExpressionsStorageId, highlightedCellIds = [], selectedGenes }: UmapVisualizationProps): ReactElement => {
    const [storageData, setStorageData] = useState<UmapStorageData | SingleCellStorageData | null>(null);
    const [geneExpressionsData, setGeneExpressionsData] = useState<GeneExpressionsStorageData | null>(null);
    const [umapPoints, setUmapPoints] = useState<UmapPoint[]>([]);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);
    const [selectedCells, setSelectedCells] = useState<string[]>([]);
    const [isSingleCellData, setIsSingleCellData] = useState(false);

    const chartRef = useRef<ChartHandle>(null);
    const containerRef = useRef<HTMLDivElement>(null);
    const { width } = useSize(containerRef);

    // Fetch main storage data (UMAP coordinates and metadata)
    useEffect(() => {
        const fetchStorageData = async () => {
            if (!storageId) {
                setError('No storage ID provided');
                return;
            }

            setIsLoading(true);
            setError(null);

            try {
                const data = await getStorage(storageId);
                
                // Check if this is single-cell data with expressions or just UMAP data
                const hasSingleCellData = data.json?.genes && data.json?.metadata?.n_genes;
                setIsSingleCellData(hasSingleCellData);
                
                setStorageData(data as unknown as UmapStorageData | SingleCellStorageData);
            } catch (err) {
                setError('Failed to load UMAP data');
                console.error('Error fetching storage data:', err);
            } finally {
                setIsLoading(false);
            }
        };

        fetchStorageData();
    }, [storageId]);

    // Fetch gene expressions data when needed
    useEffect(() => {
        const fetchGeneExpressionsData = async () => {
            if (!geneExpressionsStorageId || !isSingleCellData || selectedGenes.length === 0) {
                setGeneExpressionsData(null);
                return;
            }

            try {
                const data = await getStorage(geneExpressionsStorageId);
                setGeneExpressionsData(data as unknown as GeneExpressionsStorageData);
            } catch (err) {
                console.error('Error fetching gene expressions data:', err);
                // Don't set error state as this is optional functionality
            }
        };

        fetchGeneExpressionsData();
    }, [geneExpressionsStorageId, isSingleCellData, selectedGenes]);

    // Update UMAP points with expression data
    useEffect(() => {
        if (!storageData) return;

        const cells = storageData.json?.cells || {};
        
        const points: UmapPoint[] = Object.entries(cells).map(([cellId, cellData]) => {
            const typedCellData = cellData as any;
            const basePoint = {
                cellId,
                umap_1: typedCellData.umap_1,
                umap_2: typedCellData.umap_2,
            };
            
            // Calculate mean expression for selected genes if we have expression data
            if (isSingleCellData && selectedGenes.length > 0 && geneExpressionsData) {
                const expressions: number[] = [];
                
                selectedGenes.forEach(gene => {
                    const geneData = geneExpressionsData.json[gene.feature_id];
                    if (geneData && geneData.cells[cellId]) {
                        expressions.push(geneData.cells[cellId]);
                    }
                });
                
                const meanExpression = expressions.length > 0 
                    ? expressions.reduce((sum, val) => sum + val, 0) / expressions.length
                    : 0;
                
                return {
                    ...basePoint,
                    expressionValue: meanExpression
                };
            }
            
            return basePoint;
        });
        
        setUmapPoints(points);
    }, [storageData, geneExpressionsData, selectedGenes, isSingleCellData]);

    const handlePlotOnSelect = useCallback((selectedCellIds: string[]) => {
        setSelectedCells(selectedCellIds);
    }, []);

    if (isLoading) {
        return (
            <UmapVisualizationContainer>
                <Box display="flex" justifyContent="center" alignItems="center" height="100%">
                    <Typography>Loading UMAP data...</Typography>
                </Box>
            </UmapVisualizationContainer>
        );
    }

    if (error) {
        return (
            <UmapVisualizationContainer>
                <Box display="flex" justifyContent="center" alignItems="center" height="100%">
                    <Typography color="error">{error}</Typography>
                </Box>
            </UmapVisualizationContainer>
        );
    }

    return (
        <>
            <UmapVisualizationContainer ref={containerRef}>
                <UmapVisualizationControls>
                    <Typography variant="h6">
                        UMAP Visualization: {storageData?.name || `Storage ${storageId}`}
                    </Typography>
                    {storageData && (
                        <Typography variant="body2" color="textSecondary">
                            {storageData.json.metadata.n_cells} cells
                            {isSingleCellData && ` • ${(storageData as SingleCellStorageData).json.metadata.n_genes} genes`}
                            {selectedGenes.length > 0 && ` • Coloring by ${selectedGenes.length} selected gene${selectedGenes.length > 1 ? 's' : ''}`}
                        </Typography>
                    )}
                </UmapVisualizationControls>
                
                {umapPoints.length > 0 && (
                    <UmapPlotContainer>
                        <UmapScatterPlot
                            data={umapPoints}
                            selectedCellIds={selectedCells}
                            highlightedCellIds={highlightedCellIds}
                            onSelect={handlePlotOnSelect}
                            colorByExpression={isSingleCellData && selectedGenes.length > 0 && geneExpressionsData !== null}
                            ref={chartRef}
                        />
                    </UmapPlotContainer>
                )}
                
                {selectedCells.length > 0 && (
                    <Box mt={2}>
                        <Typography variant="body2">
                            Selected cells: {selectedCells.join(', ')}
                        </Typography>
                    </Box>
                )}
            </UmapVisualizationContainer>
        </>
    );
};

export default connector(UmapVisualization);
