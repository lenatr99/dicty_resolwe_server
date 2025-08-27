import React, { ReactElement, useState } from 'react';
import Button from '@mui/material/Button';
import { Checkbox, FormControlLabel } from '@mui/material';
import { connect, ConnectedProps } from 'react-redux';
import {
    ModalFooter,
    FooterControlsContainer,
} from 'components/genexpress/common/dictyModal/dictyModal.styles';
import { allGenesDeselected, genesSelected } from 'redux/stores/genes';

const connector = connect(null, {
    connectedGenesSelected: genesSelected,
    connectedAllGenesDeselected: allGenesDeselected,
});

type PropsFromRedux = ConnectedProps<typeof connector>;

type GeneSelectorModalControlsProps = {
    selectedGenesIds: string[];
    allGenesIds: string[];
    onClose: () => void;
} & PropsFromRedux;

const GeneSelectorModalControls = ({
    selectedGenesIds,
    allGenesIds,
    onClose,
    connectedGenesSelected,
    connectedAllGenesDeselected,
}: GeneSelectorModalControlsProps): ReactElement => {
    const [append, setAppend] = useState(true);

    const handleAppendCheckboxChange = (event: React.ChangeEvent<HTMLInputElement>): void => {
        setAppend(event.target.checked);
    };

    const setSelectedGenes = (genesIds: string[]): void => {
        connectedAllGenesDeselected();
        connectedGenesSelected(genesIds);
        onClose();
    };

    const appendSelectedGenes = (genesIds: string[]): void => {
        connectedGenesSelected(genesIds);
        onClose();
    };

    const handleSelectOnClick = (): void => {
        if (append) {
            appendSelectedGenes(selectedGenesIds);
        } else {
            setSelectedGenes(selectedGenesIds);
        }
    };

    const handleSelectAllOnClick = (): void => {
        if (append) {
            appendSelectedGenes(allGenesIds);
        } else {
            setSelectedGenes(allGenesIds);
        }
    };

    return (
        <ModalFooter>
            <FooterControlsContainer>
                <FormControlLabel
                    control={
                        <Checkbox
                            checked={append}
                            onChange={handleAppendCheckboxChange}
                            name="append"
                        />
                    }
                    label="Append selected genes to Genes module"
                />
                <div>
                    <Button onClick={handleSelectOnClick} disabled={selectedGenesIds.length === 0}>
                        Select
                    </Button>
                    <Button onClick={handleSelectAllOnClick} disabled={allGenesIds.length === 0}>
                        Select all {allGenesIds.length}
                    </Button>
                    <Button onClick={onClose}>Close</Button>
                </div>
            </FooterControlsContainer>
        </ModalFooter>
    );
};

export default connector(GeneSelectorModalControls);
