import _ from 'lodash';
import { Action } from '@reduxjs/toolkit';
import { EMPTY, Observable, from, mergeMap } from 'rxjs';
import { sessionId } from '../api/base';
import { get, post, QueryParams } from '../api/fetch';
import { logError } from 'utils/errorUtils';
import { deserializeResponse } from 'utils/apiUtils';

const MESSAGE_CREATE = 'CREATE';
const MESSAGE_UPDATE = 'UPDATE';

export interface IdObject {
    id: number;
}

export type QueryObserver = {
    subscriptionId: string;
    baseUrl: string;
    refetch: (id: number) => Promise<unknown[]>;
    itemsUpdateHandler: (items: unknown[]) => Observable<Action>;
};

type QueryObserverResponse = {
    subscription_id: string;
};

export type WebSocketMessage = {
    change_type: string;
    subscription_id: string;
    object_id: number;
};
export type DisposeFunction = () => Promise<Response>;
export type ItemsAndDisposeFunction<T extends IdObject> = {
    items: T[];
    disposeFunction: DisposeFunction;
};

const observers: QueryObserver[] = [];

const unsubscribeObserver = (subscriptionId: string): Promise<Response> => {
    const removed = _.remove(observers, { subscriptionId });

    if (removed.length === 0) return new Promise((resolve) => resolve(new Response('')));

    return post(`${removed[0].baseUrl}/unsubscribe`, { subscription_id: subscriptionId });
};

export const clearObservers = async (): Promise<void> => {
    const unsubscribePromises = observers.map((observer) =>
        unsubscribeObserver(observer.subscriptionId),
    );

    try {
        await Promise.all(unsubscribePromises);
    } catch (error) {
        logError('Error unsubscribing from queryobserver.', error);
    }
};

export const reactiveGet = async <T extends IdObject>(
    baseUrl: string,
    params: QueryParams,
    webSocketMessageOutputReduxAction: (items: unknown[]) => Observable<Action>,
): Promise<ItemsAndDisposeFunction<T>> => {
    const items = await deserializeResponse<T[]>(await get(baseUrl, params));
    const ids = items.map((item) => item.id);

    const observerResponse = await post(`${baseUrl}/subscribe`, { ids, session_id: sessionId });
    const { subscription_id: subscriptionId } =
        await deserializeResponse<QueryObserverResponse>(observerResponse);

    const refetch = (id: number): Promise<unknown[]> =>
        new Promise((resolve, reject) => {
            get(baseUrl, { id })
                .then((resp) =>
                    deserializeResponse<T[]>(resp)
                        .then((objects) => resolve(objects))
                        .catch(reject),
                )
                .catch(reject);
        });

    observers.push({
        subscriptionId,
        baseUrl,
        refetch,
        itemsUpdateHandler: webSocketMessageOutputReduxAction,
    });

    return {
        items,
        disposeFunction: () => unsubscribeObserver(subscriptionId),
    };
};

const getObserver = (subscriptionId: string): QueryObserver | undefined => {
    return _.find(observers, { subscriptionId });
};

export const handleWebSocketMessage = (message: WebSocketMessage): Observable<Action> => {
    const observer = getObserver(message.subscription_id);
    if (observer == null) {
        return EMPTY;
    }
    const { object_id: objectId, change_type: changeType } = message;
    const { refetch } = observer;

    let refetchPromise: Promise<unknown[]> = new Promise((resolve) => resolve([]));

    if (changeType === MESSAGE_CREATE || changeType === MESSAGE_UPDATE) {
        refetchPromise = refetchPromise.then(async () => refetch(objectId));
    }

    return from(refetchPromise).pipe(mergeMap((items) => observer.itemsUpdateHandler(items)));
};
