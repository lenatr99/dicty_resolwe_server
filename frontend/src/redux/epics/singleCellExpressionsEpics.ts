import { Action } from '@reduxjs/toolkit';
import { Epic, combineEpics } from 'redux-observable';
import { map, mergeMap, startWith, endWith, catchError, filter, distinctUntilChanged, switchMap, concatMap, reduce } from 'rxjs/operators';
import { of, from, forkJoin, EMPTY } from 'rxjs';
import { Data, Storage } from '@genialis/resolwe/dist/api/types/rest';
import { mapStateSlice } from './rxjsCustomFilters';
import { RootState } from 'redux/rootReducer';
import { handleError } from 'utils/errorUtils';
import { SamplesGenesExpressionsById } from 'redux/models/internal';
import { getDataBySamplesIds, getStorage, getRelationBySlug } from 'api';
import {
    getSamplesExpressionsSamplesIds,
    samplesExpressionsFetchEnded,
    samplesExpressionsFetchStarted,
    samplesExpressionsFetchSucceeded,
} from 'redux/stores/samplesExpressions';
import { getSelectedTimeSeries } from 'redux/stores/timeSeries';
import { getSelectedGenes } from 'redux/stores/genes';
import { singleCellSeriesFetchEnded, singleCellSeriesFetchStarted, singleCellSeriesFetchSucceeded } from 'redux/stores/singleCellSeries';
import { Relation } from '@genialis/resolwe/dist/api/types/rest';
import _ from 'lodash';

// Helper to get storage for a generic entity (by its Data item), mapping by entity id
const getEntityStorage = async (
    dataItem: Data,
): Promise<{ entityId: number; storage: Storage }> => {
    const storage = await getStorage(dataItem.output.exp_json);
    const entityId = dataItem.entity != null ? dataItem.entity.id : 0;
    return { entityId, storage };
};

// Fetch single-cell relation for current time series (slug + '_sc')
const fetchSingleCellRelationEpic: Epic<Action, Action, RootState> = (_action$, state$) => {
    return state$.pipe(
        mapStateSlice((state) => {
            const ts = getSelectedTimeSeries(state.timeSeries);
            const bySlug: Record<string, Relation> = (state.singleCellSeries as any).bySlug || {};
            const scSlugDash = ts ? `${ts.slug}-sc` : '';
            const has = !!(scSlugDash && bySlug[scSlugDash]);
            return { scSlugDash, has };
        }),
        distinctUntilChanged((prev, curr) => prev.scSlugDash === curr.scSlugDash && prev.has === curr.has),
        filter(({ scSlugDash, has }) => !!scSlugDash && !has),
        switchMap(({ scSlugDash }) => from(getRelationBySlug(scSlugDash)).pipe(
            filter((rel): rel is Relation => rel != null),
            map((relation) => singleCellSeriesFetchSucceeded(relation)),
            catchError((error) => of(handleError('Error fetching single-cell relation.', error))),
            startWith(singleCellSeriesFetchStarted()),
            endWith(singleCellSeriesFetchEnded()),
        )),
    );
};

// Epic to fetch single-cell expressions (similar to time series pattern)
const fetchSingleCellExpressionsEpic: Epic<Action, Action, RootState> = (_action$, state$) => {
    return state$.pipe(
        mapStateSlice((state) => {
            const ts = getSelectedTimeSeries(state.timeSeries);
            const genes = getSelectedGenes(state.genes);
            const sc: Record<string, Relation> = (state.singleCellSeries as any).bySlug || {};
            const rel = ts ? sc[`${ts.slug}-sc`] : undefined;
            if (!ts || !rel || genes.length === 0) return { key: '', ids: [] as number[] };
            const labels = new Set<string>();
            genes.forEach((g) => { if (g.name) labels.add(g.name); if (g.feature_id) labels.add(g.feature_id); });
            const allIds = _.uniq((rel.partitions as any[])
                .filter((p) => p.label && labels.has(p.label))
                .map((p) => p.entity))
                .sort((a: number, b: number) => a - b);
            // Always consider only the first 500 candidates overall; ignore the rest entirely
            const firstBatch = allIds.slice(0, 500);
            const existing = new Set<number>(getSamplesExpressionsSamplesIds(state.samplesExpressions));
            const need = firstBatch.filter((id: number) => !existing.has(id));
            // Key is based on the fixed firstBatch so it doesn't change after those are fetched
            return { key: `${ts.slug}|${rel.id}|${firstBatch.join(',')}`, ids: need };
        }),
        distinctUntilChanged((a, b) => a.key === b.key),
        filter(({ ids }) => ids.length > 0),
        switchMap(({ ids }) => from(getDataBySamplesIds(ids)).pipe(
            // Process storages in small chunks to limit parallel requests
            mergeMap((dataItems) => {
                const MAX_STORAGES = 500;
                const CHUNK_SIZE = 50;
                const limited = dataItems.slice(0, MAX_STORAGES);
                const chunks: Array<typeof limited> = _.chunk(limited, CHUNK_SIZE) as any;
                return from(chunks).pipe(
                    concatMap((chunk) => forkJoin(chunk.map(getEntityStorage))),
                    reduce((acc, results) => acc.concat(results), [] as Array<{ entityId: number; storage: Storage }>),
                    map((storages) => {
                        const byId: SamplesGenesExpressionsById = {};
                        storages.forEach(({ entityId, storage }) => { byId[entityId] = (storage as any).json.genes; });
                        return samplesExpressionsFetchSucceeded(byId);
                    }),
                    catchError((error) => of(handleError('Error retrieving single-cell storage data.', error))),
                );
            }),
            catchError((error) => of(handleError('Error retrieving single-cell data.', error))),
            startWith(samplesExpressionsFetchStarted()),
            endWith(samplesExpressionsFetchEnded()),
        )),
    );
};

export default combineEpics(fetchSingleCellRelationEpic, fetchSingleCellExpressionsEpic);
