import _ from 'lodash';
import differentialExpressionsReducer, {
    differentialExpressionSelected,
    differentialExpressionsFetchSucceeded,
    DifferentialExpressionsState,
    differentialExpressionStorageFetchSucceeded,
} from './differentialExpressions';
import {
    generateDifferentialExpressionJson,
    generateDifferentialExpressionsById,
    generateStorage,
} from 'tests/mock';

const differentialExpressionsById = generateDifferentialExpressionsById(2);
const differentialExpressions = _.flatMap(differentialExpressionsById);
const differentialExpressionStorageJson = generateDifferentialExpressionJson(5);

describe('differentialExpressions store', () => {
    let initialState: DifferentialExpressionsState;

    describe('differentialExpressions not fetched', () => {
        beforeEach(() => {
            initialState = {
                byId: {},
                selectedId: 0,
                isFetchingDifferentialExpressions: false,
                isFetchingDifferentialExpressionsData: false,
            };
        });

        it('should add fetched differentialExpressions to state with differentialExpressionsFetchSucceeded action', () => {
            const newState = differentialExpressionsReducer(
                initialState,
                differentialExpressionsFetchSucceeded(differentialExpressions),
            );
            const expectedState = {
                ...initialState,
                byId: differentialExpressionsById,
            };

            expect(newState).toEqual(expectedState);
        });

        it('should set selectedId with differentialExpressionSelected action', () => {
            const newState = differentialExpressionsReducer(
                initialState,
                differentialExpressionSelected(123),
            );
            const expectedState = {
                ...initialState,
                selectedId: 123,
            };

            expect(newState).toEqual(expectedState);
        });
    });

    describe('differentialExpressions fetched', () => {
        beforeEach(() => {
            initialState = {
                byId: differentialExpressionsById,
                selectedId: 0,
                isFetchingDifferentialExpressions: false,
                isFetchingDifferentialExpressionsData: false,
            };
        });

        it('should add fetched differentialExpressions to state with differentialExpressionsFetchSucceeded action', () => {
            const newState = differentialExpressionsReducer(
                initialState,
                differentialExpressionStorageFetchSucceeded(
                    generateStorage(
                        differentialExpressions[0].output.de_json,
                        differentialExpressionStorageJson,
                    ),
                ),
            );
            const expectedState = {
                ...initialState,
                byId: {
                    [differentialExpressions[0].id]: {
                        ...differentialExpressions[0],
                        json: differentialExpressionStorageJson,
                    },
                    [differentialExpressions[1].id]: differentialExpressions[1],
                },
            };

            expect(newState).toEqual(expectedState);
        });
    });
});
