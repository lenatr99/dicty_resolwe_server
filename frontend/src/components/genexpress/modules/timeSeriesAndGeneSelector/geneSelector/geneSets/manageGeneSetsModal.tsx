import React, { ReactElement, useState } from 'react';
import Button from '@mui/material/Button';
import { ColDef, ValueFormatterParams } from 'ag-grid-community';
import { ManageGeneSetsGridWrapper } from './manageGeneSetsModal.styles';
import {
    ModalFooter,
    ModalBody,
    ModalHeader,
    ModalContainer,
    CenteredModal,
} from 'components/genexpress/common/dictyModal/dictyModal.styles';
import DictyGrid from 'components/genexpress/common/dictyGrid/dictyGrid';
import { GeneSet } from 'redux/models/internal';
import useStateWithEffect from 'components/genexpress/common/useStateWithEffect';

type ManageGeneSetsModalProps = {
    geneSets: GeneSet[];
    onDelete: (geneSets: GeneSet[]) => void;
    onClick: (geneSet: GeneSet) => void;
    onClose: () => void;
};

const checkboxColumnDef = {
    headerCheckboxSelection: true,
    checkboxSelection: true,
    width: 25,
};

const geneSetColumnDefs = [
    {
        field: 'dateTime',
        headerName: 'Date',
        width: 90,
        sort: 'desc',
        valueFormatter: (params: ValueFormatterParams): string => {
            return params.value.toLocaleString('en-US');
        },
    },
    {
        field: 'genesNames',
        headerName: 'Genes',
        autoHeight: true,
        cellStyle: { whiteSpace: 'normal' },
        valueFormatter: (params: ValueFormatterParams): string => {
            return params.value.join(', ');
        },
    },
] as ColDef[];

const ManageGeneSetsModal = ({
    geneSets,
    onDelete,
    onClick,
    onClose,
}: ManageGeneSetsModalProps): ReactElement => {
    const [selectedGeneSets, setSelectedGeneSets] = useState<GeneSet[]>([]);
    const [isEditMode, setIsEditMode] = useState(false);
    const geneSetsSelectionChangedHandler = (newSelectedGeneSets: GeneSet[]): void => {
        setSelectedGeneSets(newSelectedGeneSets);
    };

    const handleOnDelete = (): void => {
        onDelete(selectedGeneSets);
    };

    const handleOnClick = (geneSet: GeneSet): void => {
        onClick(geneSet);
    };

    const columnDefs = useStateWithEffect(() => {
        return isEditMode ? [checkboxColumnDef, ...geneSetColumnDefs] : geneSetColumnDefs;
    }, [isEditMode]);

    return (
        <CenteredModal
            open
            aria-labelledby="modalTitle"
            aria-describedby="modalDescription"
            onClose={onClose}
        >
            <ModalContainer>
                <ModalHeader id="modalTitle">Gene List History</ModalHeader>
                <ModalBody>
                    <div id="modalDescription">
                        Choose a gene set from a list of locally saved gene sets
                    </div>
                    <ManageGeneSetsGridWrapper>
                        <DictyGrid
                            data={geneSets}
                            getRowId={(data): string => data.dateTime.toString()}
                            columnDefs={columnDefs}
                            selectionMode="multiple"
                            onSelectionChanged={geneSetsSelectionChangedHandler}
                            onRowClicked={handleOnClick}
                            suppressRowClickSelection
                            hideFilter
                        />
                    </ManageGeneSetsGridWrapper>
                </ModalBody>
                <ModalFooter>
                    {isEditMode ? (
                        <Button onClick={handleOnDelete} disabled={selectedGeneSets.length === 0}>
                            Delete selected
                        </Button>
                    ) : (
                        <Button
                            onClick={() => {
                                setIsEditMode(true);
                            }}
                            disabled={geneSets.length === 0}
                        >
                            Edit
                        </Button>
                    )}

                    <Button onClick={onClose}>Close</Button>
                </ModalFooter>
            </ModalContainer>
        </CenteredModal>
    );
};

export default ManageGeneSetsModal;
