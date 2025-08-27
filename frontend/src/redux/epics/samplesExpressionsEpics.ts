import { Action } from '@reduxjs/toolkit';
import { Epic, combineEpics } from 'redux-observable';
import { map, mergeMap, startWith, endWith, catchError } from 'rxjs/operators';
import { of, from, forkJoin, EMPTY } from 'rxjs';
import { Data, Storage } from '@genialis/resolwe/dist/api/types/rest';
import { mapStateSlice } from './rxjsCustomFilters';
import { getAllTimeSeriesSamplesIds } from 'redux/stores/timeSeries';
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

const getSampleStorage = async (
    sampleData: Data,
): Promise<{ sampleId: number; storage: Storage }> => {
    const storage = await getStorage(sampleData.output.exp_json);

    return {
        sampleId: sampleData.entity != null ? sampleData.entity.id : 0,
        storage,
    };
};

const fetchSamplesExpressionsEpic: Epic<Action, Action, RootState> = (_action$, state$) => {
    return state$.pipe(
        mapStateSlice(
            (state) => getAllTimeSeriesSamplesIds(state.timeSeries),
            (timeSeriesSamplesIds) => timeSeriesSamplesIds.length > 0,
        ),
        mergeMap((timeSeriesSamplesIds) => {
            const samplesExpressionsInStore = getSamplesExpressionsSamplesIds(
                state$.value.samplesExpressions,
            );

            const timeSeriesSamplesIdsToFetch = timeSeriesSamplesIds.filter(
                (sampleId) => !samplesExpressionsInStore.includes(sampleId),
            );

            if (timeSeriesSamplesIdsToFetch.length === 0) {
                return EMPTY;
            }

            return from(getDataBySamplesIds(timeSeriesSamplesIdsToFetch)).pipe(
                mergeMap((sampleData) => {
                    // Once samples data is retrieved use it's output.exp_json to retrieve genes expressions.
                    return forkJoin(sampleData.map(getSampleStorage)).pipe(
                        map((sampleStorages) => {
                            const timeSeriesSamplesExpressions = {} as SamplesGenesExpressionsById;
                            sampleStorages.forEach(({ sampleId, storage }) => {
                                timeSeriesSamplesExpressions[sampleId] = storage.json.genes;
                            });

                            return samplesExpressionsFetchSucceeded(timeSeriesSamplesExpressions);
                        }),
                        catchError((error) =>
                            of(handleError(`Error retrieving samples storage data.`, error)),
                        ),
                    );
                }),
                catchError((error) =>
                    of(handleError(`Error retrieving samples storage data.`, error)),
                ),
                startWith(samplesExpressionsFetchStarted()),
                endWith(samplesExpressionsFetchEnded()),
            );
        }),
    );
};

export default combineEpics(fetchSamplesExpressionsEpic);
