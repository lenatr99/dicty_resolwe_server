import samplesExpressionsReducer, {
    SamplesExpressionsState,
    samplesExpressionsFetchSucceeded,
    samplesExpressionsComparisonFetchSucceeded,
} from './samplesExpressions';
import { timeSeriesSelected } from './timeSeries';
import { generateSamplesExpressionsById } from 'tests/mock';

const samplesExpressionsById = generateSamplesExpressionsById(2, []);

describe('samplesExpressions store', () => {
    let initialState: SamplesExpressionsState;

    describe('empty initial state', () => {
        beforeEach(() => {
            initialState = {
                byId: {},
                isFetchingSamplesExpressions: false,
            };
        });

        it('should set fetched samplesExpressions to state with samplesExpressionsFetchSucceeded action', () => {
            const newState = samplesExpressionsReducer(
                initialState,
                samplesExpressionsFetchSucceeded(samplesExpressionsById),
            );
            const expectedState = {
                ...initialState,
                byId: samplesExpressionsById,
            };

            expect(newState).toEqual(expectedState);
        });

        it('should add fetched samplesExpressions to state with samplesExpressionsComparisonFetchSucceeded action', () => {
            const comparisonsSamplesExpressionsById = generateSamplesExpressionsById(2, []);
            const newState = samplesExpressionsReducer(
                {
                    ...initialState,
                    byId: samplesExpressionsById,
                },
                samplesExpressionsComparisonFetchSucceeded(comparisonsSamplesExpressionsById),
            );
            const expectedState = {
                ...initialState,
                byId: { ...samplesExpressionsById, ...comparisonsSamplesExpressionsById },
            };

            expect(newState).toEqual(expectedState);
        });
    });

    describe('not empty initial state', () => {
        beforeEach(() => {
            initialState = {
                byId: samplesExpressionsById,
                isFetchingSamplesExpressions: false,
            };
        });

        it('should clear mergedData on timeSeriesSelected action', () => {
            const newState = samplesExpressionsReducer(initialState, timeSeriesSelected(1));
            const expectedState = {
                ...initialState,
                byId: {},
            };

            expect(newState).toEqual(expectedState);
        });
    });
});
