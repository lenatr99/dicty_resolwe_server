import React, { ReactElement, useState } from 'react';
import { Button } from '@mui/material';
import SaveIcon from '@mui/icons-material/Save';
import { useDispatch } from 'react-redux';
import ManageGeneSetsModal from './manageGeneSetsModal';
import { GeneSet, Gene } from 'redux/models/internal';
import { LocalStorageKey } from 'components/genexpress/common/constants';
import useLocalStorage from 'components/genexpress/common/useLocalStorage';
import IconButtonWithTooltip from 'components/genexpress/common/iconButtonWithTooltip/iconButtonWithTooltip';
import { addInfoSnackbar } from 'redux/stores/notifications';

type SelectGeneSetProps = {
    selectedGenes: Gene[];
    disabled?: boolean;
    onSelect: (genesNames: Gene['name'][]) => void;
};

const GeneSetSelector = ({
    selectedGenes,
    disabled,
    onSelect,
}: SelectGeneSetProps): ReactElement => {
    const [manageModalOpened, setManageModalOpened] = useState(false);
    const [localStorageGeneSets, localStorageSetGeneSets] = useLocalStorage<GeneSet[]>(
        LocalStorageKey.geneSets,
        [
            {
                dateTime: new Date(),
                genesNames: ['pkaR', 'pkgB', 'pkaC', 'ppk1', 'pks5', 'pks6'],
            },
        ],
    );
    const dispatch = useDispatch();

    /**
     * Saves gene set to local storage so user can later use it from history.
     */
    const handleSaveOnClick = (): void => {
        if (selectedGenes != null && selectedGenes.length > 0) {
            const newGeneSet = {
                dateTime: new Date(),
                genesNames: selectedGenes.map((gene) => gene.name),
            } as GeneSet;

            localStorageSetGeneSets((prevGeneSets) => [...prevGeneSets, newGeneSet]);

            dispatch(addInfoSnackbar('Gene set saved.'));
        }
    };

    const handleOnClick = (geneSet: GeneSet): void => {
        setManageModalOpened(false);
        onSelect(geneSet.genesNames);
    };

    const handleOnDelete = (geneSetsToDelete: GeneSet[]): void => {
        const updatedGeneSets = localStorageGeneSets.filter(
            (geneSet) => !geneSetsToDelete.includes(geneSet),
        );
        localStorageSetGeneSets(updatedGeneSets);
    };

    return (
        <>
            <div>
                <IconButtonWithTooltip
                    title="Save current gene set to local storage"
                    disabled={selectedGenes.length === 0}
                    onClick={handleSaveOnClick}
                >
                    <SaveIcon />
                </IconButtonWithTooltip>
                <Button
                    size="small"
                    onClick={(): void => setManageModalOpened(true)}
                    disabled={disabled}
                >
                    History
                </Button>
            </div>

            {manageModalOpened && (
                <ManageGeneSetsModal
                    geneSets={localStorageGeneSets}
                    onClick={handleOnClick}
                    onDelete={handleOnDelete}
                    onClose={(): void => setManageModalOpened(false)}
                />
            )}
        </>
    );
};

export default GeneSetSelector;
