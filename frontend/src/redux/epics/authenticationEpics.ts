import { Action } from '@reduxjs/toolkit';
import { Epic, combineEpics } from 'redux-observable';
import { map, startWith, endWith, catchError, mergeMap, filter } from 'rxjs/operators';
import { of, from } from 'rxjs';
import { appStarted } from './epicsActions';
import { RootState } from 'redux/rootReducer';
import { userFetchSucceeded, userFetchStarted, userFetchEnded } from 'redux/stores/authentication';
import { handleError } from 'utils/errorUtils';
import { getCurrentUser } from 'api';
import { setSentryUser } from 'utils/sentryUtils';

const getCurrentUserEpic: Epic<Action, Action, RootState> = (action$) =>
    action$.pipe(
        filter(appStarted.match),
        mergeMap(() => {
            return from(getCurrentUser()).pipe(
                map((user) => {
                    // Propagate user to sentry context, so future sentry logs will include
                    // user information.
                    setSentryUser(user);
                    return userFetchSucceeded(user);
                }),
                catchError((error) => of(handleError(`Error fetching user profile`, error))),
                startWith(userFetchStarted()),
                endWith(userFetchEnded()),
            );
        }),
    );

export default combineEpics(getCurrentUserEpic);
