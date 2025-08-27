import React, { ReactElement, useRef, useState } from 'react';
import Popper from '@mui/material/Popper';
import ClickAwayListener from '@mui/material/ClickAwayListener';
import { GeneChipWrapper } from './geneChip.styles';
import { Gene } from 'redux/models/internal';
import ConnectedGeneInformation from 'components/genexpress/modules/timeSeriesAndGeneSelector/geneSelector/selectedGenes/geneInformation/geneInformation';

export type GeneChipProps = {
    gene: Gene;
    highlighted: boolean;
    onHighlight: () => void;
    onUnhighlight: () => void;
    onRemove: () => void;
};
export const GeneChip = ({
    gene,
    onRemove,
    highlighted,
    onHighlight,
    onUnhighlight,
}: GeneChipProps): ReactElement => {
    const [infoPopperOpen, setInfoPopperOpen] = useState(false);
    const chipWrapperElement = useRef<HTMLDivElement>(null);

    const handleClickAway = (): void => {
        setInfoPopperOpen(false);
    };

    const handleOnClick = (): void => {
        setInfoPopperOpen((prev) => !prev);
    };

    return (
        <ClickAwayListener onClickAway={handleClickAway}>
            <div>
                <GeneChipWrapper
                    label={gene.name}
                    onDelete={onRemove}
                    sx={{ fontSize: '1rem' }}
                    ref={chipWrapperElement}
                    onClick={handleOnClick}
                    size="small"
                    color={highlighted ? 'secondary' : undefined}
                />

                <Popper anchorEl={chipWrapperElement.current} open={infoPopperOpen}>
                    <ConnectedGeneInformation
                        gene={gene}
                        highlighted={highlighted}
                        onHighlight={onHighlight}
                        onUnhighlight={onUnhighlight}
                    />
                </Popper>
            </div>
        </ClickAwayListener>
    );
};

export default GeneChip;
