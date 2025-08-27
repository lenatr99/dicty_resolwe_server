import { useEffect, useState, useCallback } from 'react';
import {
    getLocalStorageKey,
    writeToLocalStorage,
    readFromLocalStorage,
    deleteFromLocalStorage,
} from '../../../utils/localStorageUtils';

/**
 * React hook to enable updates to state via localStorage.
 * This updates when the {writeStorage} function is used, when the returned function
 * is called.
 * This function takes an optional default value to start off with.
 *
 * @example
 * ```js
 * const MyComponent = () => {
 *   const [myStoredItem, setMyStoredItem] = useLocalStorage('myStoredItem', 'initialValue');
 *   return (
 *     <p>{myStoredItem}</p>
 *   );
 * };
 * ```
 *
 * @export
 * @template TValue The type of the given initial value.
 * @param {string} targetKey The key in the localStorage that you wish to watch.
 * @param {TValue} initialValue Initial value to start with.
 * @returns {[TValue, Dispatch<TValue>, Dispatch<void>]} An array containing the value
 * associated with the key in position 0, a function to set the value in position 1,
 * and a function to delete the value from localStorage in position 2.
 */
const useLocalStorage = <TValue = string>(
    targetKey: string,
    initialValue: TValue,
): [TValue, (value: TValue | ((prevState: TValue) => TValue)) => void, () => void] => {
    const key = getLocalStorageKey(targetKey);

    const [localState, updateLocalState] = useState<TValue>(
        readFromLocalStorage(key) ?? initialValue,
    );

    // Each time local state changes, write it's new value to local storage.
    useEffect(() => {
        writeToLocalStorage(key, localState);
    }, [key, localState]);

    useEffect(() => {
        const itemExists = readFromLocalStorage(key) != null;

        // Write initial value to the local storage if it's not present or contains invalid JSON data.
        if (initialValue !== undefined && !itemExists) {
            writeToLocalStorage(key, initialValue);
        }
    }, [initialValue, key]);

    const deleteState = useCallback(() => {
        updateLocalState(initialValue);
        deleteFromLocalStorage(key);
    }, [initialValue, key]);

    return [localState === null ? initialValue : localState, updateLocalState, deleteState];
};

export default useLocalStorage;
