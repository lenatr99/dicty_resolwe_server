import { Button } from '@mui/material';
import styled from 'styled-components';

export const MatchedCellContainer = styled.div`
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: center;
    position: relative;
    width: 100%;
    height: 100%;
`;

export const MatchedButton = styled(Button)`
    height: 26px;
    font-weight: bold;
`;
