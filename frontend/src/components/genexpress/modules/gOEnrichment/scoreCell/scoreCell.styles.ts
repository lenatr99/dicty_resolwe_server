import { LinearProgress } from '@mui/material';
import styled from 'styled-components';

export const ScoreCellContainer = styled.div`
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: center;
    position: relative;
    width: 100%;
    height: 100%;
`;

export const ScoreSpan = styled.span`
    position: absolute;
    top: 0;
    font-weight: bold;
`;

export const ScoreLinearProgress = styled(LinearProgress)`
    && {
        width: 100%;
        height: 26px;
        background-color: white;

        & > .MuiLinearProgress-bar1Determinate {
            background-color: rgb(185, 197, 211);
        }
    }
`;
