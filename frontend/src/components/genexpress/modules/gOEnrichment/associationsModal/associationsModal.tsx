import React, { ReactElement, useEffect, useState } from 'react';
import { connect, ConnectedProps } from 'react-redux';
import { ColDef, ValueGetterParams } from 'ag-grid-community';
import { Tooltip } from '@mui/material';
import {
    TermInfo,
    TermName,
    AssociationsGridWrapper,
    AmigoLink,
    AmigoLinkImage,
} from './associationsModal.style';
import {
    ModalBody,
    ModalHeader,
    ModalContainer,
    CenteredModal,
} from 'components/genexpress/common/dictyModal/dictyModal.styles';
import DictyGrid from 'components/genexpress/common/dictyGrid/dictyGrid';
import GeneSelectorModalControls from 'components/genexpress/common/geneSelectorModalControls/geneSelectorModalControls';
import { BasketInfo, Gene, GOEnrichmentRow } from 'redux/models/internal';
import amigoLogo from 'images/amigo_logo.png';
import {
    getGenesById,
    getIsFetchingAssociationsGenes,
    getSelectedGenesIds,
} from 'redux/stores/genes';
import { RootState } from 'redux/rootReducer';
import { fetchAssociationsGenes } from 'redux/epics/epicsActions';
import { mapGeneIdsBetweenSources } from 'api/kbApi';
import { getBasketInfo } from 'redux/stores/timeSeries';

const mapStateToProps = (state: RootState) => {
    return {
        genesById: getGenesById(state.genes),
        selectedGenesIds: getSelectedGenesIds(state.genes),
        isFetchingAssociationsGenes: getIsFetchingAssociationsGenes(state.genes),
        basketInfo: getBasketInfo(state.timeSeries) as BasketInfo,
    };
};

const connector = connect(mapStateToProps, {
    connectedFetchAssociationsGenes: fetchAssociationsGenes,
});

type PropsFromRedux = ConnectedProps<typeof connector>;

type GOEnrichmentAssociationsModalProps = {
    gOEnrichmentRow: GOEnrichmentRow;
    handleOnClose: () => void;
} & PropsFromRedux;

const columnDefs = [
    {
        headerCheckboxSelection: true,
        checkboxSelection: true,
        width: 25,
    },
    {
        valueGetter: (params: ValueGetterParams): string => {
            return params.data.feature_id;
        },
        headerName: 'ID',
        width: 150,
    },
    {
        valueGetter: (params: ValueGetterParams): string => {
            return params.data.name;
        },
        headerName: 'Name',
        sort: 'asc',
        flex: 1,
    },
] as ColDef[];

const GOEnrichmentAssociationsModal = ({
    gOEnrichmentRow,
    genesById,
    selectedGenesIds,
    handleOnClose,
    connectedFetchAssociationsGenes,
    isFetchingAssociationsGenes,
    basketInfo,
}: GOEnrichmentAssociationsModalProps): ReactElement => {
    const [associatedGenes, setAssociatedGenes] = useState<Gene[]>([]);
    const [selectedAssociatedGenes, setSelectedAssociatedGenes] = useState<Gene[]>([]);
    const [genesMappingsSourceIds, setGenesMappingsSourceIds] = useState<Gene['feature_id'][]>();
    const [isFetchingMappings, setIsFetchingMappings] = useState(true);

    useEffect(() => {
        const mapGenesAndTriggerFetch = async () => {
            setIsFetchingMappings(true);

            setGenesMappingsSourceIds(
                (
                    await mapGeneIdsBetweenSources({
                        targetGenesIds: gOEnrichmentRow.gene_associations,
                        sourceDb: basketInfo.source,
                        sourceSpecies: basketInfo.species,
                    })
                ).map((geneMapping) => geneMapping.source_id),
            );

            setIsFetchingMappings(false);
        };

        void mapGenesAndTriggerFetch();
    }, [basketInfo.source, basketInfo.species, gOEnrichmentRow.gene_associations]);

    // When modal opens, fetch all genes that are associated with clicked gene ontology enrichment row.
    useEffect(() => {
        if (genesMappingsSourceIds != null) {
            connectedFetchAssociationsGenes({
                geneIds: genesMappingsSourceIds,
            });
        }
    }, [connectedFetchAssociationsGenes, genesMappingsSourceIds]);

    // Prepare data -> attach Gene to each volcano point.
    useEffect(() => {
        if (genesMappingsSourceIds != null) {
            setAssociatedGenes(genesMappingsSourceIds.flatMap((geneId) => genesById[geneId] ?? []));
        }
    }, [genesById, genesMappingsSourceIds]);

    // Propagate already selected genes to selected volcano points.
    useEffect(() => {
        setSelectedAssociatedGenes(
            associatedGenes.filter((gene) => selectedGenesIds.includes(gene.feature_id)),
        );
    }, [associatedGenes, selectedGenesIds]);

    return (
        <CenteredModal
            open
            aria-labelledby="modalTitle"
            aria-describedby="modalDescription"
            onClose={handleOnClose}
        >
            <ModalContainer>
                <ModalHeader id="modalTitle">GO Term Association</ModalHeader>
                <ModalBody>
                    <TermInfo id="modalDescription">
                        <TermName>{gOEnrichmentRow.term_name}</TermName>
                        <Tooltip title="Open in AmiGO">
                            <AmigoLink
                                href={`http://amigo.geneontology.org/amigo/term/${gOEnrichmentRow.term_id}`}
                                target="_blank"
                                rel="noopener noreferrer"
                            >
                                {gOEnrichmentRow.term_id}
                                <AmigoLinkImage src={amigoLogo} alt="AmiGO link" />
                            </AmigoLink>
                        </Tooltip>
                    </TermInfo>
                    <p>
                        {gOEnrichmentRow.total === 1
                            ? `1 gene is associated with ${gOEnrichmentRow.term_id} term`
                            : `${gOEnrichmentRow.total} genes are associated with ${gOEnrichmentRow.term_id} term`}
                    </p>
                    <AssociationsGridWrapper>
                        <DictyGrid
                            isFetching={isFetchingAssociationsGenes || isFetchingMappings}
                            data={associatedGenes}
                            getRowId={(data): string => data.feature_id}
                            filterLabel="Filter"
                            selectedData={selectedAssociatedGenes}
                            columnDefs={columnDefs}
                            selectionMode="multiple"
                            onSelectionChanged={setSelectedAssociatedGenes}
                        />
                    </AssociationsGridWrapper>
                </ModalBody>
                <GeneSelectorModalControls
                    allGenesIds={associatedGenes.map((gene) => gene.feature_id)}
                    selectedGenesIds={selectedAssociatedGenes.map((gene) => gene.feature_id)}
                    onClose={handleOnClose}
                />
            </ModalContainer>
        </CenteredModal>
    );
};

export default connector(GOEnrichmentAssociationsModal);
