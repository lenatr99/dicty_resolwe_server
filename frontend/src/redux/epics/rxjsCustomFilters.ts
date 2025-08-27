import { Observable } from 'rxjs';
import { distinctUntilChanged, filter, map } from 'rxjs/operators';
import { RootState } from 'redux/rootReducer';

const inputIsNotNullOrUndefined = <T>(input: null | undefined | T): input is T => {
    return input !== null && input !== undefined;
};

/**
 * TypeScript doesn't know how to correctly infer type with filter operator:
 * - https://github.com/Microsoft/TypeScript/issues/10734
 * - https://github.com/Microsoft/TypeScript/issues/16069.
 *
 * That's why this custom operator is used instead.
 */
export const filterNullAndUndefined = <T>() => {
    return (source$: Observable<null | undefined | T>): Observable<T> =>
        source$.pipe(filter(inputIsNotNullOrUndefined));
};

export const mapStateSlice = <T>(
    selector: (state: RootState) => T,
    // Filter if target state is already set (e.g. unit tests initialState).
    existenceFilter?: (value: T) => boolean,
) => {
    return (state$: Observable<RootState>): Observable<T> =>
        state$.pipe(
            map((state) => selector(state)),
            distinctUntilChanged(),
            filter((value) => existenceFilter?.(value) ?? true),
        );
};
