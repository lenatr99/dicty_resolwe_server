import React, { ReactElement } from 'react';
import { connect, ConnectedProps } from 'react-redux';
import { Clear as ClearIcon, ContentCopy as ContentCopyIcon } from '@mui/icons-material';
import { GeneChip } from './geneChip/geneChip';
import { SelectedGenesContainer, ActionsContainer } from './selectedGenes.styles';
import { Gene } from 'redux/models/internal';
import {
    geneDeselected,
    geneUnhighlighted,
    geneHighlighted,
    allGenesDeselected,
} from 'redux/stores/genes';
import IconButtonWithTooltip from 'components/genexpress/common/iconButtonWithTooltip/iconButtonWithTooltip';
import { setClipboardText } from 'utils/documentHelpers';

const connector = connect(null, {
    connectedGeneDeselected: geneDeselected,
    connectedGeneHighlighted: geneHighlighted,
    connectedGeneUnhighlighted: geneUnhighlighted,
    connectedAllGenesDeselected: allGenesDeselected,
});

type PropsFromRedux = ConnectedProps<typeof connector>;

type SelectedGenesProps = {
    selectedGenes: Gene[];
    highlightedGenesIds: string[];
} & PropsFromRedux;

const SelectedGenes = ({
    selectedGenes,
    highlightedGenesIds,
    connectedGeneDeselected,
    connectedGeneHighlighted,
    connectedGeneUnhighlighted,
    connectedAllGenesDeselected,
}: SelectedGenesProps): ReactElement => {
    const handleOnRemove = (gene: Gene): void => {
        connectedGeneDeselected(gene.feature_id);
    };

    const handleOnHighlight = (gene: Gene): void => {
        connectedGeneHighlighted(gene.feature_id);
    };

    const handleOnUnhighlight = (gene: Gene): void => {
        connectedGeneUnhighlighted(gene.feature_id);
    };

    const handleCopyClick = (): void => {
        setClipboardText(selectedGenes.map((gene) => gene.name).join(', '));
    };

    const handleClearAllClick = (): void => {
        connectedAllGenesDeselected();
    };

    return (
        <>
            <SelectedGenesContainer>
                {selectedGenes.map((gene) => (
                    <GeneChip
                        key={gene.feature_id}
                        gene={gene}
                        highlighted={highlightedGenesIds.includes(gene.feature_id)}
                        onRemove={(): void => handleOnRemove(gene)}
                        onHighlight={(): void => handleOnHighlight(gene)}
                        onUnhighlight={(): void => handleOnUnhighlight(gene)}
                    />
                ))}
            </SelectedGenesContainer>
            <ActionsContainer>
                <IconButtonWithTooltip
                    title={`Copy ${selectedGenes.length} ${
                        selectedGenes.length === 1 ? 'gene' : 'genes'
                    } to clipboard`}
                    disabled={selectedGenes.length === 0}
                    onClick={handleCopyClick}
                >
                    <ContentCopyIcon />
                </IconButtonWithTooltip>
                <IconButtonWithTooltip
                    title="Clear all"
                    disabled={selectedGenes.length === 0}
                    onClick={handleClearAllClick}
                >
                    <ClearIcon />
                </IconButtonWithTooltip>
            </ActionsContainer>
        </>
    );
};

export default connector(SelectedGenes);
