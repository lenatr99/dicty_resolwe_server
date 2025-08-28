import styled from 'styled-components';

export const UmapVisualizationContainer = styled.div`
    display: flex;
    flex-flow: column nowrap;
    height: 100%;
    gap: ${({ theme }) => theme.spacing(1)};
`;

export const UmapVisualizationControls = styled.div`
    display: flex;
    flex-flow: row nowrap;
    align-items: center;
    justify-content: space-between;
    padding: ${({ theme }) => theme.spacing(1)};
`;

export const UmapPlotContainer = styled.div`
    flex-grow: 1;
    overflow: hidden;
    min-height: 400px;
`;
