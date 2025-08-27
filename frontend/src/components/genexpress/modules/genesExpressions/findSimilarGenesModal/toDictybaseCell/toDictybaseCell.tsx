import React, { ReactElement } from 'react';
import { ICellRendererParams } from 'ag-grid-community';
import { CellContainer, CellLogo } from './toDictybaseCell.styles';
import dictyBaseLogo from 'images/dictybase_logo2.jpg';
import sacgbLogo from 'images/sacgb_logo.png';

const ToDictybaseCell = ({
    value,
    data,
}: Pick<ICellRendererParams, 'value' | 'data'>): ReactElement => (
    <CellContainer>
        <span>{value}</span>
        <a
            href={`http://dictybase.org/gene/${data.feature_id}`}
            rel="noopener noreferrer"
            target="_blank"
        >
            <CellLogo src={dictyBaseLogo} alt="Open in dictyBase." />
        </a>
        <a
            href={`http://sacgb.fli-leibniz.de/cgi/freesearch.pl?ssi=free&word=${data.feature_id}`}
            rel="noopener noreferrer"
            target="_blank"
        >
            <CellLogo src={sacgbLogo} alt="Open in SACGB." />
        </a>
    </CellContainer>
);

export default ToDictybaseCell;
