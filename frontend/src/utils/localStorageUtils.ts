import LZString from 'lz-string';
import { dateTimeReviver } from './dateTimeUtils';

/**
 * localStorageUtils handle communication with user local storage. All objects are stringified and
 * stored as zipped strings. These strings are unzipped and parsed when we read from local storage.
 */

/**
 * Each item in local storage has to have a unique key. Key is build prefix 'dictyexpress_'
 * and item key.
 * @param key - Item key.
 */
export const getLocalStorageKey = (key: string): string => {
    return `dictyexpress_${key}`;
};

/**
 * Stringifies object and saves it zipped to local storage.
 * @param key - Item key.
 * @param value - Object to save.
 */
export const writeToLocalStorage = (key: string, value: unknown): void => {
    let stringifiedObject = '';

    if (typeof value === 'object') {
        stringifiedObject = JSON.stringify(value);
    } else if (typeof value === 'string') {
        stringifiedObject = `${value}`;
    }

    return localStorage.setItem(getLocalStorageKey(key), LZString.compress(stringifiedObject));
};

/**
 * Reads zipped string from local storage, unzips it and returns it as object.
 * @param key - Item key.
 */
// eslint-disable-next-line @typescript-eslint/no-explicit-any
export const readFromLocalStorage = (key: string): any => {
    const zippedString = localStorage.getItem(getLocalStorageKey(key));
    if (zippedString == null || zippedString === '') {
        return null;
    }

    const unzippedString = LZString.decompress(zippedString);

    if (unzippedString == null) {
        return null;
    }

    return JSON.parse(unzippedString, dateTimeReviver);
};

/**
 * Removes item from local storage.
 * @param key - Item key.
 */
export const deleteFromLocalStorage = (key: string): void => {
    localStorage.removeItem(key);
};
