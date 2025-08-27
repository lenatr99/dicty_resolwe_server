import React, { ReactElement } from 'react';
import { ICellRendererParams } from 'ag-grid-community';
import { TermCellContainer, TermIndentation, TermCollapseIcon } from './termCell.styles';
import { GOEnrichmentRow } from 'redux/models/internal';

const TermCell = ({
    data,
    onToggleCollapseClick,
}: Pick<ICellRendererParams, 'data'> & {
    onToggleCollapseClick: (row: GOEnrichmentRow) => void;
}): ReactElement => (
    <TermCellContainer onClick={() => onToggleCollapseClick(data)}>
        <TermIndentation $depth={data.depth} />
        <TermCollapseIcon $collapsed={data.collapsed} $hasChildren={data.children?.length > 0} />
        <span>{data.term_name}</span>
    </TermCellContainer>
);

export default TermCell;
