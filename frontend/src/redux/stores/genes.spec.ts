import _ from 'lodash';
import genesReducer, {
    genesSelected,
    GenesState,
    geneDeselected,
    allGenesDeselected,
    geneHighlighted,
    genesHighlighted,
    geneUnhighlighted,
    genesFetchSucceeded,
} from './genes';
import { Gene, GenesById } from 'redux/models/internal';
import { generateGenesById, generateGene } from 'tests/mock';

describe('genes store', () => {
    let initialState: GenesState;
    let genesById: GenesById;
    let genes: Gene[];

    beforeEach(() => {
        genesById = generateGenesById(2);
    });

    describe('empty initial state', () => {
        beforeEach(() => {
            initialState = {
                byId: {},
                selectedGenesIds: [],
                highlightedGenesIds: [],
                isFetchingDifferentialExpressionGenes: false,
                isFetchingAssociationsGenes: false,
                isFetchingSimilarGenes: false,
                isFetchingBookmarkedGenes: false,
            };

            genes = _.flatMap(genesById);
        });

        it('should add a gene to byId with genesFetchSucceeded action', () => {
            const newState = genesReducer(initialState, genesFetchSucceeded(genes));
            const expectedState = {
                ...initialState,
                byId: genesById,
            };

            expect(newState).toEqual(expectedState);
        });

        it('should add selected genes feature_ids to selectedGenesId with genesSelected action', () => {
            const newState = genesReducer(
                initialState,
                genesSelected(genes.map((gene) => gene.feature_id)),
            );
            const expectedState = {
                ...initialState,
                selectedGenesIds: genes.map((gene) => gene.feature_id),
            };

            expect(newState).toEqual(expectedState);
        });

        it('should highlight a gene', () => {
            const newState = genesReducer(initialState, geneHighlighted(genes[0].name));
            const expectedState = { ...initialState, highlightedGenesIds: [genes[0].name] };

            expect(newState).toEqual(expectedState);
        });

        it('should highlight multiple genes', () => {
            const newState = genesReducer(
                initialState,
                genesHighlighted(genes.map((gene) => gene.feature_id)),
            );
            const expectedState = {
                ...initialState,
                highlightedGenesIds: genes.map((gene) => gene.feature_id),
            };

            expect(newState).toEqual(expectedState);
        });
    });

    describe('not empty initial state, both genes selected and highlighted', () => {
        beforeEach(() => {
            initialState = {
                byId: genesById,
                selectedGenesIds: genes.map((gene) => gene.feature_id),
                highlightedGenesIds: genes.map((gene) => gene.feature_id),
                isFetchingDifferentialExpressionGenes: false,
                isFetchingAssociationsGenes: false,
                isFetchingSimilarGenes: false,
                isFetchingBookmarkedGenes: false,
            };
        });

        it('should add genes to byId with genesFetchSucceeded action', () => {
            const newGene = generateGene(5);
            const newState = genesReducer(initialState, genesFetchSucceeded([newGene]));
            const expectedState = {
                ...initialState,
                byId: { ...initialState.byId, [newGene.feature_id]: newGene },
            };

            expect(newState).toEqual(expectedState);
        });

        it('should add selected genes feature_ids to selectedGenesId with genesSelected action', () => {
            const genesByIdToSelect = generateGenesById(2);

            const newState = genesReducer(
                initialState,
                genesSelected(Object.keys(genesByIdToSelect)),
            );
            const expectedState = {
                ...initialState,
                selectedGenesIds: [
                    ...initialState.selectedGenesIds,
                    ...Object.keys(genesByIdToSelect),
                ],
            };

            expect(newState).toEqual(expectedState);
        });

        it('should remove gene from selected and highlighted genes ids with geneDeselected action', () => {
            const newState = genesReducer(initialState, geneDeselected(genes[0].feature_id));
            const expectedState = {
                ...initialState,
                selectedGenesIds: [genes[1].feature_id],
                highlightedGenesIds: [genes[1].feature_id],
            };

            expect(newState).toEqual(expectedState);
        });

        it('should clear selected and highlighted genes ids with deselectAll action', () => {
            const newState = genesReducer(initialState, allGenesDeselected());
            const expectedState = {
                ...initialState,
                selectedGenesIds: [],
                highlightedGenesIds: [],
            };

            expect(newState).toEqual(expectedState);
        });

        it('should highlight a gene', () => {
            const newGene = generateGene(5);
            const newState = genesReducer(initialState, geneHighlighted(newGene.feature_id));
            const expectedState = {
                ...initialState,
                highlightedGenesIds: [...initialState.highlightedGenesIds, newGene.feature_id],
            };

            expect(newState).toEqual(expectedState);
        });

        it('should highlight multiple genes', () => {
            const genesByIdToHighlight = generateGenesById(2);
            const newState = genesReducer(
                initialState,
                genesHighlighted(Object.keys(genesByIdToHighlight)),
            );
            const expectedState = {
                ...initialState,
                highlightedGenesIds: Object.keys(genesByIdToHighlight),
            };

            expect(newState).toEqual(expectedState);
        });

        it('should unhighlight gene', () => {
            const newState = genesReducer(initialState, geneUnhighlighted(genes[0].feature_id));
            const expectedState = { ...initialState, highlightedGenesIds: [genes[1].feature_id] };

            expect(newState).toEqual(expectedState);
        });
    });
});
