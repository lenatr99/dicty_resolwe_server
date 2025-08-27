import { Dispatch, SetStateAction, useEffect, useState } from 'react';
import { BookmarkStatePath } from './constants';
import * as bookmarkStateManager from 'managers/bookmarkStateManager';

const useBookmarkableState = <T>(
    initialValue: T | (() => T),
    bookmarkStatePath: BookmarkStatePath,
): [T, Dispatch<SetStateAction<T>>] => {
    const [value, setValue] = useState<T>(initialValue);

    useEffect(
        () =>
            bookmarkStateManager.register({
                bookmarkStatePath,
                get: () => value,
                set: (bookmarkedValue) => {
                    setValue(bookmarkedValue);
                },
            }),
        [bookmarkStatePath, value],
    );

    return [value, setValue];
};

export default useBookmarkableState;
