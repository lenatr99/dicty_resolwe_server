import { Action } from '@reduxjs/toolkit';
import { Epic } from 'redux-observable';
import {
    catchError,
    retryWhen,
    delay,
    mergeMap,
    filter,
    mergeWith,
    endWith,
    switchMap,
    take,
} from 'rxjs/operators';
import { defer, iif, of, throwError } from 'rxjs';
import { webSocket, WebSocketSubject } from 'rxjs/webSocket';
import { filterNullAndUndefined } from './rxjsCustomFilters';
import { appFocused, appStarted, reconnectToServer } from './epicsActions';
import { RootState } from 'redux/rootReducer';
import { handleError } from 'utils/errorUtils';
import { sessionId, webSocketUrl } from 'api/base';
import { handleWebSocketMessage, WebSocketMessage } from 'managers/queryObserverManager';

const reconnectionTimeout = 6000;
const reconnectionMaxAttempts = 3;

class ConnectToWebSocketError extends Error {}

let webSocketSubject: WebSocketSubject<WebSocketMessage>;

const connectToServerEpic: Epic<Action, Action, RootState> = (action$) =>
    action$.pipe(
        filter(appStarted.match),
        mergeWith(
            action$.pipe(
                filter(reconnectToServer.match),
                switchMap(() => {
                    // Reconnect to server when page is / becomes visible.
                    if (document.visibilityState === 'visible') {
                        return of(true);
                    }

                    return defer(() => {
                        return action$.pipe(filter(appFocused.match), take(1));
                    });
                }),
            ),
        ),
        mergeMap(() => {
            webSocketSubject = webSocket({
                url: `${webSocketUrl}/${sessionId}`,
                WebSocketCtor: WebSocket,
                closeObserver: {
                    next: () => {
                        webSocketSubject?.complete();
                    },
                },
            });

            return webSocketSubject.pipe(
                retryWhen((errors) => {
                    return errors.pipe(
                        mergeMap((e, i) => {
                            return iif(
                                () => i + 1 >= reconnectionMaxAttempts,
                                throwError(() => new ConnectToWebSocketError()),
                                of(e).pipe(delay(reconnectionTimeout)),
                            );
                        }),
                    );
                }),
                mergeMap((message) => handleWebSocketMessage(message)),
                endWith(reconnectToServer()),
                catchError((err) => {
                    return of(
                        handleError(
                            err instanceof ConnectToWebSocketError
                                ? `WebSocket connection failed after ${reconnectionMaxAttempts} tries.`
                                : `WebSocket message handler error.`,
                            err,
                        ),
                    );
                }),
            );
        }),
        filterNullAndUndefined(),
    );

export default connectToServerEpic;
