import { Action } from '@reduxjs/toolkit';
import { Epic, combineEpics } from 'redux-observable';
import { map, mergeMap, startWith, endWith, catchError, filter, switchMap } from 'rxjs/operators';
import { of, from, combineLatest } from 'rxjs';
import { isEmpty, compact } from 'lodash';
import { appStarted } from './epicsActions';
import { filterNullAndUndefined, mapStateSlice } from './rxjsCustomFilters';
import {
    getSelectedTimeSeriesSamplesIds,
    addToBasketStarted,
    addToBasketEnded,
    addSamplesToBasketSucceeded,
    timeSeriesFetchStarted,
    timeSeriesFetchEnded,
    timeSeriesFetchSucceeded,
    fetchBasketExpressionsIdsSucceeded,
    getBasketInfo,
    getBasketExpressionsIds,
    timeSeriesSelected,
    getTimeSeries,
    getComparisonTimeSeries,
    genesMappingsFetchSucceeded,
    genesMappingsFetchStarted,
    genesMappingsFetchEnded,
    getComparisonTimeSeriesIds,
    getSelectedTimeSeries,
} from 'redux/stores/timeSeries';
import { RootState } from 'redux/rootReducer';
import { handleError } from 'utils/errorUtils';
import { addToBasket, getTimeSeriesRelations, getBasketExpressions } from 'api';
import { mapGeneIdsBetweenSources } from 'api/kbApi';
import { getSelectedGenesIds } from 'redux/stores/genes';
import { BasketInfo } from 'redux/models/internal';

const fetchTimeSeriesEpic: Epic<Action, Action, RootState> = (action$) => {
    return action$.pipe(
        filter(appStarted.match),
        mergeMap(() => {
            return from(getTimeSeriesRelations()).pipe(
                map((response) => timeSeriesFetchSucceeded(response)),
                catchError((error) => of(handleError('Error retrieving time series.', error))),
                startWith(timeSeriesFetchStarted()),
                endWith(timeSeriesFetchEnded()),
            );
        }),
    );
};

const selectDefaultTimeSeriesEpic: Epic<Action, Action, RootState> = (action$, state$) =>
    state$.pipe(
        mapStateSlice((state) => getTimeSeries(state.timeSeries)),
        filter(
            (timeSeries) =>
                timeSeries.length > 0 && getSelectedTimeSeries(state$.value.timeSeries) == null,
        ),
        map((timeSeries) => {
            const defaultTimeSeriesSlug =
                typeof SELECTED_TIMESERIES_SLUG !== 'undefined' && SELECTED_TIMESERIES_SLUG != null
                    ? timeSeries.find((relation) => relation.slug === SELECTED_TIMESERIES_SLUG)
                    : undefined;
            return timeSeriesSelected((defaultTimeSeriesSlug ?? timeSeries[0]).id);
        }),
    );

const timeSeriesSelectedEpic: Epic<Action, Action, RootState> = (action$, state$) => {
    return state$.pipe(
        mapStateSlice(
            (state) => getSelectedTimeSeriesSamplesIds(state.timeSeries),
            (timeSeriesSamplesIds) => timeSeriesSamplesIds.length > 0,
        ),
        mergeMap((timeSeriesSamplesIds) => {
            return from(addToBasket(timeSeriesSamplesIds)).pipe(
                map((response) => {
                    if (
                        isEmpty(compact(response.permitted_sources)) ||
                        isEmpty(compact(response.permitted_organisms))
                    ) {
                        throw new Error(
                            'Selected time series are missing source and species information.',
                        );
                    }
                    return addSamplesToBasketSucceeded(response);
                }),
                catchError((error) =>
                    of(
                        handleError(
                            'Error adding time series samples to visualization basket.',
                            error,
                        ),
                    ),
                ),
                startWith(addToBasketStarted()),
                endWith(addToBasketEnded()),
            );
        }),
    );
};

const comparisonTimeSeriesSelectedEpic: Epic<Action, Action, RootState> = (action$, state$) => {
    return combineLatest([
        state$.pipe(
            mapStateSlice((state) => {
                return getComparisonTimeSeriesIds(state.timeSeries);
            }),
        ),
        state$.pipe(
            mapStateSlice((state) => {
                return getSelectedGenesIds(state.genes);
            }),
        ),
    ]).pipe(
        filter(([allComparisonTimeSeriesIds, selectedGenesIds]) => {
            return allComparisonTimeSeriesIds.length > 0 && selectedGenesIds.length > 0;
        }),
        switchMap(([, selectedGenesIds]) => {
            const allComparisonTimeSeries = getComparisonTimeSeries(state$.value.timeSeries);

            return from(allComparisonTimeSeries).pipe(
                filter(
                    (comparisonTimeSeries) =>
                        comparisonTimeSeries.genesMappings == null ||
                        selectedGenesIds.some(
                            (geneId) =>
                                comparisonTimeSeries.genesMappings?.find(
                                    (geneMapping) => geneMapping.source_id === geneId,
                                ) == null,
                        ),
                ),
                switchMap((comparisonTimeSeries) => {
                    if (comparisonTimeSeries.basketInfo == null) {
                        return from(
                            addToBasket(
                                comparisonTimeSeries?.partitions.map(
                                    (partition) => partition.entity,
                                ) ?? [],
                            ),
                        ).pipe(
                            map((response) => ({
                                comparisonTimeSeries,
                                comparisonTimeSeriesBasketInfo: {
                                    id: response.id,
                                    source: response.permitted_sources[0],
                                    species: response.permitted_organisms[0],
                                    type: 'gene',
                                } as BasketInfo,
                            })),
                        );
                    } else {
                        return of({
                            comparisonTimeSeries,
                            comparisonTimeSeriesBasketInfo: comparisonTimeSeries.basketInfo,
                        });
                    }
                }),
                switchMap(({ comparisonTimeSeries, comparisonTimeSeriesBasketInfo }) => {
                    // Fetch genes mappings only for genes that weren't mapped yet.
                    return from(
                        mapGeneIdsBetweenSources({
                            sourceGenesIds: selectedGenesIds.filter(
                                (geneId) =>
                                    comparisonTimeSeries.genesMappings?.find(
                                        (geneMapping) => geneMapping.source_id === geneId,
                                    ) == null,
                            ),
                            targetDb: comparisonTimeSeriesBasketInfo.source,
                            targetSpecies: comparisonTimeSeriesBasketInfo.species,
                        }),
                    ).pipe(
                        map((genesMappings) => {
                            return genesMappingsFetchSucceeded({
                                timeSeriesId: comparisonTimeSeries.id,
                                genesMappings,
                                basketInfo: comparisonTimeSeriesBasketInfo,
                            });
                        }),
                        catchError((error) =>
                            of(
                                handleError(
                                    'Error fetching comparison time series genes mappings.',
                                    error,
                                ),
                            ),
                        ),
                    );
                }),
                catchError((error) => {
                    return of(
                        handleError(
                            'Error adding comparison time series samples to visualization basket.',
                            error,
                        ),
                    );
                }),
                startWith(genesMappingsFetchStarted()),
                endWith(genesMappingsFetchEnded()),
            );
        }),
    );
};

const fetchBasketExpressionsEpic: Epic<Action, Action, RootState> = (action$, state$) => {
    return state$.pipe(
        mapStateSlice(
            (state) => getBasketInfo(state.timeSeries),
            () => getBasketExpressionsIds(state$.value.timeSeries).length === 0,
        ),
        filterNullAndUndefined(),
        mergeMap((basketInfo) => {
            return from(getBasketExpressions(basketInfo.id)).pipe(
                map((basketExpressions) => {
                    return fetchBasketExpressionsIdsSucceeded(
                        basketExpressions.map((basketExpression) => basketExpression.id),
                    );
                }),
                catchError((error) => of(handleError('Error fetching basket expressions.', error))),
            );
        }),
    );
};

export default combineEpics(
    timeSeriesSelectedEpic,
    fetchTimeSeriesEpic,
    selectDefaultTimeSeriesEpic,
    fetchBasketExpressionsEpic,
    comparisonTimeSeriesSelectedEpic,
);
