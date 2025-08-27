import _ from 'lodash';
import genesSimilaritiesReducer, {
    genesSimilaritiesDistanceMeasureChanged,
    genesSimilaritiesQueryGeneSelected,
    genesSimilaritiesFetchSucceeded,
    GenesSimilaritiesState,
} from './genesSimilarities';
import { DistanceMeasure } from 'components/genexpress/common/constants';
import { generateGenesById, generateGeneSimilarity } from 'tests/mock';

const genesById = generateGenesById(5);
const genes = _.flatMap(genesById);
const genesSimilaritiesJson = genes.map((gene) => generateGeneSimilarity(gene.feature_id));

describe('genesSimilarities store', () => {
    let initialState: GenesSimilaritiesState;

    describe('empty initial state', () => {
        beforeEach(() => {
            initialState = {
                status: null,
                data: null,
                queryGeneId: null,
                distanceMeasure: DistanceMeasure.euclidean,
                isFetchingGenesSimilarities: false,
            };
        });

        it('should add fetched json to state with genesSimilaritiesFetchSucceeded action', () => {
            const newState = genesSimilaritiesReducer(
                initialState,
                genesSimilaritiesFetchSucceeded(genesSimilaritiesJson),
            );
            const expectedState = {
                ...initialState,
                data: genesSimilaritiesJson,
            };

            expect(newState).toEqual(expectedState);
        });

        it('should change distanceMeasure with genesSimilaritiesDistanceMeasureChanged action', () => {
            const newState = genesSimilaritiesReducer(
                initialState,
                genesSimilaritiesDistanceMeasureChanged(DistanceMeasure.pearson),
            );
            const expectedState = {
                ...initialState,
                distanceMeasure: DistanceMeasure.pearson,
            };

            expect(newState).toEqual(expectedState);
        });

        it('should change queryGeneId with genesSimilaritiesQueryGeneSelected action', () => {
            const newState = genesSimilaritiesReducer(
                initialState,
                genesSimilaritiesQueryGeneSelected('123'),
            );
            const expectedState = {
                ...initialState,
                queryGeneId: '123',
            };

            expect(newState).toEqual(expectedState);
        });
    });

    describe('not empty initial state', () => {
        beforeEach(() => {
            initialState = {
                status: null,
                data: genesSimilaritiesJson,
                queryGeneId: genes[1].feature_id,
                distanceMeasure: DistanceMeasure.spearman,
                isFetchingGenesSimilarities: false,
            };
        });

        it('should clear data on genesSimilaritiesDistanceMeasureChanged action', () => {
            const newState = genesSimilaritiesReducer(
                initialState,
                genesSimilaritiesDistanceMeasureChanged(DistanceMeasure.pearson),
            );
            const expectedState = {
                ...initialState,
                distanceMeasure: 'pearson',
                data: null,
            };

            expect(newState).toEqual(expectedState);
        });

        it('should clear data on genesSimilaritiesQueryGeneSelected action', () => {
            const newState = genesSimilaritiesReducer(
                initialState,
                genesSimilaritiesQueryGeneSelected('123'),
            );
            const expectedState = {
                ...initialState,
                queryGeneId: '123',
                data: null,
            };

            expect(newState).toEqual(expectedState);
        });
    });
});
