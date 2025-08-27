import _ from 'lodash';
import timeSeriesReducer, {
    timeSeriesSelected,
    TimeSeriesState,
    timeSeriesFetchSucceeded,
    addSamplesToBasketSucceeded,
    fetchBasketExpressionsIdsSucceeded,
    genesMappingsFetchSucceeded,
} from './timeSeries';
import {
    generateTimeSeriesById,
    generateBasketAddSamplesResponse,
    generateSingleTimeSeries,
    generateBasketExpression,
} from 'tests/mock';
import { BasketInfo, GeneMapping } from 'redux/models/internal';

const timeSeriesById = generateTimeSeriesById(2);
const timeSeries = _.flatMap(timeSeriesById);
const basket = generateBasketAddSamplesResponse('123');
const basketInfo: BasketInfo = {
    id: basket.id,
    source: basket.permitted_sources[0],
    species: basket.permitted_organisms[0],
    type: 'gene',
};
const basketExpressions = [generateBasketExpression(), generateBasketExpression()];
const basketExpressionsIds = basketExpressions.map((basketExpression) => basketExpression.id);
const genesMappings: GeneMapping[] = [
    { source_db: 'DICTYBASE', source_id: 'gene1', target_db: 'UniProtKB', target_id: 'Q55eR8' },
    { source_db: 'DICTYBASE', source_id: 'gene2', target_db: 'UniProtKB', target_id: 'Q54PU2' },
    { source_db: 'DICTYBASE', source_id: 'gene3', target_db: 'UniProtKB', target_id: 'Q54L63' },
];

describe('timeSeries store', () => {
    let initialState: TimeSeriesState;

    describe('empty initial state', () => {
        beforeEach(() => {
            initialState = {
                byId: {},
                selectedId: 0,
                comparisonIds: [],
                isFetching: false,
                isAddingToBasket: false,
                isFetchingGenesMappings: false,
                basketInfo: null,
                basketExpressionsIds: [],
            };
        });

        it('should add fetched timeSeries to state with timeSeriesFetchSucceeded action', () => {
            const newState = timeSeriesReducer(initialState, timeSeriesFetchSucceeded(timeSeries));
            const expectedState = {
                ...initialState,
                byId: timeSeriesById,
            };

            expect(newState).toEqual(expectedState);
        });

        it('should set timeSeries geneMappings with timeSeriesFetchSucceeded action', () => {
            const newState = timeSeriesReducer(
                { ...initialState, byId: timeSeriesById },
                genesMappingsFetchSucceeded({
                    timeSeriesId: timeSeries[0].id,
                    genesMappings,
                    basketInfo,
                }),
            );
            const expectedState = {
                ...initialState,
                byId: {
                    ...timeSeriesById,
                    [timeSeries[0].id]: { ...timeSeries[0], genesMappings, basketInfo },
                },
            };

            expect(newState).toEqual(expectedState);
        });

        it('should set selectedId with timeSeriesSelected action', () => {
            const newState = timeSeriesReducer(initialState, timeSeriesSelected(123));
            const expectedState = {
                ...initialState,
                selectedId: 123,
            };

            expect(newState).toEqual(expectedState);
        });

        it('should set basketInfo with addSamplesToBasketSucceeded action', () => {
            const newState = timeSeriesReducer(initialState, addSamplesToBasketSucceeded(basket));
            const expectedState = {
                ...initialState,
                basketInfo,
            };

            expect(newState).toEqual(expectedState);
        });

        it('should set basketExpressionsIds with fetchBasketExpressionsIdsSucceeded action', () => {
            const newState = timeSeriesReducer(
                initialState,
                fetchBasketExpressionsIdsSucceeded(basketExpressionsIds),
            );
            const expectedState = {
                ...initialState,
                basketExpressionsIds,
            };

            expect(newState).toEqual(expectedState);
        });
    });

    describe('not empty initial state', () => {
        beforeEach(() => {
            initialState = {
                byId: timeSeriesById,
                selectedId: 2,
                comparisonIds: [],
                isFetching: false,
                isAddingToBasket: false,
                isFetchingGenesMappings: false,
                basketInfo,
                basketExpressionsIds,
            };
        });

        it('should set new fetched timeSeries to state with timeSeriesFetchSucceeded action', () => {
            const newRelation = generateSingleTimeSeries(2);
            const newState = timeSeriesReducer(
                initialState,
                timeSeriesFetchSucceeded([newRelation]),
            );
            const expectedState = {
                ...initialState,
                byId: { [newRelation.id]: newRelation },
            };

            expect(newState).toEqual(expectedState);
        });

        it('should clear basketInfo, basketExpressionsIds on timeSeriesSelected action', () => {
            const newState = timeSeriesReducer(initialState, timeSeriesSelected(1));
            const expectedState = {
                ...initialState,
                selectedId: 1,
                basketInfo: null,
                basketExpressionsIds: [],
            };

            expect(newState).toEqual(expectedState);
        });
    });
});
