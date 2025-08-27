import { configureStore, Action } from '@reduxjs/toolkit';
import { createEpicMiddleware } from 'redux-observable';
import { switchMap } from 'rxjs/operators';
import { BehaviorSubject } from 'rxjs';
import rootReducer, { RootState } from './rootReducer';
import rootEpic from './epics/rootEpic';

/* In production mode, immutable-state-invariant and serializable-state-invariant-middleware
 * middlewares are disabled by default. If you want the same performance in development, disable
 * immutableCheck and serializableCheck in getDefaultMiddleware function:
 * getDefaultMiddelwareOptions = { immutableCheck: false, serializableCheck: false }.
 */
type MiddleWareProps = {
    thunk?: boolean;
    immutableCheck?: boolean;
    serializableCheck?: boolean;
};

const getStore = (initialState?: RootState, defaultMiddlewareOptions?: MiddleWareProps) => {
    const epicMiddleware = createEpicMiddleware<Action, Action, RootState>();

    const store = configureStore({
        reducer: rootReducer,
        middleware: (getDefaultMiddleware) => {
            return getDefaultMiddleware({
                serializableCheck: {
                    // Ignore these action types
                    // ignoredActions: ['your/action/type'],
                    // Ignore these field paths in all actions
                    ignoredActionPaths: ['payload.navigate'],
                    // Ignore these paths in the state
                    // ignoredPaths: ['items.dates'],
                },
                ...defaultMiddlewareOptions,
            }).concat(epicMiddleware);
        },
        preloadedState: initialState,
        devTools: import.meta.env.DEV,
    });

    const epic$ = new BehaviorSubject(rootEpic);
    // Every time a new epic is given to epic$ it
    // will unsubscribe from the previous one then
    // call and subscribe to the new one because of
    // how switchMap works
    epicMiddleware.run((actionIn, actionOut, state) =>
        epic$.pipe(switchMap((epic) => epic(actionIn, actionOut, state))),
    );

    if (import.meta.hot) {
        import.meta.hot.accept('./rootReducer', () => {
            store.replaceReducer(rootReducer);
        });

        import.meta.hot.accept('./epics/rootEpic', () => {
            epic$.next(rootEpic);
        });
    }

    return store;
};

export default getStore;
