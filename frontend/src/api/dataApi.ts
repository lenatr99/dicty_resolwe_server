import { Data, DONE_DATA_STATUS } from '@genialis/resolwe/dist/api/types/rest';
import { Action } from '@reduxjs/toolkit';
import { Observable } from 'rxjs';
import { deserializeResponse } from '../utils/apiUtils';
import { apiUrl } from './base';
import { get } from './fetch';
import {
    DisposeFunction as QueryObserverDisposeFunction,
    reactiveGet,
    IdObject,
} from 'managers/queryObserverManager';
import { handleError } from 'utils/errorUtils';

const baseUrl = `${apiUrl}/data`;

export const getGafs = async (): Promise<Data[]> => {
    const getGafsResponse = await get(baseUrl, {
        type: 'data:gaf',
        status: DONE_DATA_STATUS,
    });

    return deserializeResponse<Data[]>(getGafsResponse);
};

export const getOntologyObo = async (): Promise<Data> => {
    const getDataBySlugResponse = await get(baseUrl, {
        slug: 'gene-ontology-core',
    });

    return (await deserializeResponse<Data[]>(getDataBySlugResponse))[0];
};

export const getDataBySamplesIds = async (samplesIds: number[]): Promise<Data[]> => {
    if (samplesIds.length === 0) {
        return [] as Data[];
    }

    const getSamplesDataResponse = await get(baseUrl, {
        type: 'data:expression',
        entity__in: samplesIds.join(','),
    });

    return deserializeResponse<Data[]>(getSamplesDataResponse);
};

export const getDataReactive = async <T extends IdObject>(
    dataId: number,
    handleDataResponse: (items: T) => Observable<Action>,
): Promise<{ item: T; disposeFunction: QueryObserverDisposeFunction }> => {
    const webSocketMessageOutputReduxAction = (
        items: unknown[],
    ): Observable<ReturnType<typeof handleError> | Action> => {
        return handleDataResponse(items[0] as T);
    };

    const { items, disposeFunction } = await reactiveGet<T>(
        baseUrl,
        { id: dataId },
        webSocketMessageOutputReduxAction,
    );
    return { item: items[0], disposeFunction };
};
