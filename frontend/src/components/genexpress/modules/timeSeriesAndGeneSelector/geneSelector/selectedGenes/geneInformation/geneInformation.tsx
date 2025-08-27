import React, { ReactElement } from 'react';
import Button from '@mui/material/Button';
import {
    GeneInformationContainer,
    GeneInformationTitle,
    GeneInformationValue,
    GeneInformationLabel,
    DictyBaseLogo,
    GeneInformationHeader,
} from './geneInformation.styles';
import { Gene } from 'redux/models/internal';
import dictyBaseLogo from 'images/dictybase_logo2.jpg';

export type GeneInformationProps = {
    gene: Gene;
    highlighted?: boolean;
    onHighlight?: () => void;
    onUnhighlight?: () => void;
};
export const GeneInformation = ({
    gene,
    onHighlight,
    onUnhighlight,
    highlighted,
}: GeneInformationProps): ReactElement => {
    return (
        <GeneInformationContainer>
            <GeneInformationHeader>
                <GeneInformationTitle>Gene information</GeneInformationTitle>
                {highlighted ? (
                    <Button
                        onClick={onUnhighlight}
                        color="secondary"
                        variant="contained"
                        size="small"
                    >
                        Unhighlight
                    </Button>
                ) : (
                    <Button
                        onClick={onHighlight}
                        color="secondary"
                        variant="contained"
                        size="small"
                    >
                        Highlight
                    </Button>
                )}
            </GeneInformationHeader>
            <p>{gene.name}</p>
            <div>
                <GeneInformationLabel>Gene ID:</GeneInformationLabel>
                <GeneInformationValue>{gene.feature_id}</GeneInformationValue>
            </div>
            <div>
                <GeneInformationLabel>Alternatives:</GeneInformationLabel>
                <GeneInformationValue>{gene.aliases?.join(', ')}</GeneInformationValue>
            </div>
            <div>
                <GeneInformationLabel>Database link:</GeneInformationLabel>
                <GeneInformationValue>
                    <a
                        href={`http://dictybase.org/gene/${gene.feature_id}`}
                        rel="noopener noreferrer"
                        target="_blank"
                    >
                        Open in dictyBase
                        <DictyBaseLogo src={dictyBaseLogo} alt="Open in dictyBase." />
                    </a>
                </GeneInformationValue>
            </div>
            <div>
                <GeneInformationLabel>Description:</GeneInformationLabel>
                <GeneInformationValue>{gene.description}</GeneInformationValue>
            </div>
        </GeneInformationContainer>
    );
};

export default GeneInformation;
