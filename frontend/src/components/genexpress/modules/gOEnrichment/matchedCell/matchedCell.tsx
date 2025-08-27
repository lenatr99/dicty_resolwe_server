import React, { ReactElement } from 'react';
import { ICellRendererParams } from 'ag-grid-community';
import { MatchedButton, MatchedCellContainer } from './matchedCell.styles';
import { GOEnrichmentRow } from 'redux/models/internal';

const MatchedCell = ({
    value,
    data,
    onMatchedGenesClick,
}: Pick<ICellRendererParams, 'value' | 'data'> & {
    onMatchedGenesClick: (row: GOEnrichmentRow) => void;
}): ReactElement => (
    <MatchedCellContainer>
        <MatchedButton
            fullWidth
            color="secondary"
            onClick={(): void => {
                onMatchedGenesClick(data);
            }}
        >
            {`${value}/${data.total}`}
        </MatchedButton>
    </MatchedCellContainer>
);

export default MatchedCell;
