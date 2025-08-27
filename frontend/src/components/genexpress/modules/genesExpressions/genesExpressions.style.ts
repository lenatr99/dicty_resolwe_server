import { Button } from '@mui/material';
import styled from 'styled-components';

export const GenesExpressionsContainer = styled.div`
    display: flex;
    flex-flow: column nowrap;
    height: 100%;
`;

export const GenesExpressionsControls = styled.div`
    display: flex;
    justify-content: space-between;
    flex-flow: row wrap;
`;

export const GenesExpressionsLineChartContainer = styled.div`
    flex-grow: 1;
    overflow: hidden;
`;

export const FindSimilarGenesButton = styled(Button)`
    align-self: flex-start;
`;
