import { Action } from '@reduxjs/toolkit';
import { Epic, combineEpics } from 'redux-observable';
import { map, mergeMap, startWith, endWith, catchError } from 'rxjs/operators';
import { of, from, forkJoin, EMPTY } from 'rxjs';
import { Data, Storage } from '@genialis/resolwe/dist/api/types/rest';
import { mapStateSlice } from './rxjsCustomFilters';
import { RootState } from 'redux/rootReducer';
import { handleError } from 'utils/errorUtils';
import { SamplesGenesExpressionsById } from 'redux/models/internal';
import { getDataBySamplesIds, getStorage } from 'api';
import {
    getSamplesExpressionsSamplesIds,
    samplesExpressionsFetchEnded,
    samplesExpressionsFetchStarted,
    samplesExpressionsFetchSucceeded,
} from 'redux/stores/samplesExpressions';

// Helper function to get cell expression storage
const getCellStorage = async (
    cellData: Data,
): Promise<{ cellId: string; storage: Storage }> => {
    const storage = await getStorage(cellData.output.exp_json);
    
    return {
        cellId: cellData.output.cell_id,
        storage,
    };
};

// Epic to fetch single-cell expressions (similar to time series pattern)
const fetchSingleCellExpressionsEpic: Epic<Action, Action, RootState> = (_action$, state$) => {
    return state$.pipe(
        mapStateSlice(
            // This would need to be implemented - get all single-cell sample IDs from the relation
            (state) => getAllSingleCellSamplesIds(state.singleCellSeries), // Similar to getAllTimeSeriesSamplesIds
        ),
        mergeMap((singleCellSamplesIds) => {
            const samplesExpressionsInStore = getSamplesExpressionsSamplesIds(
                state$.value.samplesExpressions,
            );

            const singleCellSamplesIdsToFetch = singleCellSamplesIds.filter(
                (sampleId) => !samplesExpressionsInStore.includes(sampleId),
            );

            if (singleCellSamplesIdsToFetch.length === 0) {
                return EMPTY;
            }

            return from(getDataBySamplesIds(singleCellSamplesIdsToFetch)).pipe(
                mergeMap((cellData) => {
                    // Once cell data is retrieved, use its output.exp_json to retrieve gene expressions
                    return forkJoin(cellData.map(getCellStorage)).pipe(
                        map((cellStorages) => {
                            const singleCellSamplesExpressions = {} as SamplesGenesExpressionsById;
                            cellStorages.forEach(({ cellId, storage }) => {
                                // Map cell ID to entity ID for consistency with time series pattern
                                const entityId = cellData.find(d => d.output.cell_id === cellId)?.entity?.id || 0;
                                singleCellSamplesExpressions[entityId] = storage.json.genes;
                            });

                            return samplesExpressionsFetchSucceeded(singleCellSamplesExpressions);
                        }),
                        catchError((error) =>
                            of(handleError(`Error retrieving single-cell storage data.`, error)),
                        ),
                    );
                }),
                catchError((error) =>
                    of(handleError(`Error retrieving single-cell data.`, error)),
                ),
                startWith(samplesExpressionsFetchStarted()),
                endWith(samplesExpressionsFetchEnded()),
            );
        }),
    );
};

// Helper function to get all single-cell samples IDs (would need to be implemented)
function getAllSingleCellSamplesIds(singleCellSeries: any): number[] {
    // This would extract entity IDs from the single-cell relation partitions
    // Similar to how getAllTimeSeriesSamplesIds works
    return [];
}

export default combineEpics(fetchSingleCellExpressionsEpic);
