import React, { ReactElement } from 'react';
import { ICellRendererParams } from 'ag-grid-community';
import { ScoreCellContainer, ScoreSpan, ScoreLinearProgress } from './scoreCell.styles';
import { formatNumber } from 'utils/math';

const ScoreCell = ({ value, data }: Pick<ICellRendererParams, 'value' | 'data'>): ReactElement => (
    <ScoreCellContainer>
        <ScoreLinearProgress
            variant="determinate"
            value={data.score_percentage}
            color="secondary"
        />
        <ScoreSpan>{formatNumber(value, 'short')}</ScoreSpan>
    </ScoreCellContainer>
);

export default ScoreCell;
