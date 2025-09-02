import styled from 'styled-components';

export const UmapVisualizationContainer = styled.div`
    height: 100%;
    display: flex;
    flex-direction: column;
    padding: 16px;
`;

export const UmapControlsContainer = styled.div`
    display: flex;
    gap: 16px;
    margin-bottom: 16px;
    align-items: center;
`;

export const UmapChartContainer = styled.div`
    flex: 1;
    min-height: 400px;
    width: 100%;
    background: #ffffff;
    border-radius: 4px;
    border: 1px solid #e0e0e0;
`;

export const UmapStatusMessage = styled.div`
    padding: 20px;
    text-align: center;
    color: #666;
    font-size: 14px;
`;

export const UmapErrorMessage = styled.div`
    padding: 20px;
    text-align: center;
    color: #f44336;
    font-size: 14px;
    background-color: #ffebee;
    border-radius: 4px;
    border: 1px solid #f44336;
`;

export const UmapStatsContainer = styled.div`
    display: flex;
    gap: 20px;
    margin-bottom: 16px;
    padding: 12px;
    background-color: #f5f5f5;
    border-radius: 4px;
    font-size: 14px;
`;

export const UmapStat = styled.div`
    display: flex;
    flex-direction: column;
    
    .label {
        font-weight: bold;
        color: #333;
        margin-bottom: 4px;
    }
    
    .value {
        color: #666;
    }
`;

