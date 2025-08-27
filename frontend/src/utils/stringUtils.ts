import _ from 'lodash';

/**
 * Converts commas and newlines into spaces. Only keeps a-z A-Z 0-9 _.-().
 * @param {string} dirtyString - Original string (user input).
 */
export const splitAndCleanGenesString = (dirtyString: string): string[] => {
    let str = dirtyString;
    if (!str) return [];
    str = str.replace(/[,\n\t]/g, ' ');

    // Remove other signs.
    // eslint-disable-next-line no-useless-escape
    str = str.replace(/[^a-zA-Z0-9\s_\.\-\(\)\/]/g, '');

    // Don't split by punctuation space (\u2008), other & chars are already removed.
    str = str.replace(/\u2008/g, '&');

    let arr = str.split(/\s/) || [];
    arr = _.map(arr, (v) => {
        return v.replace(/&/g, '\u2008').trim();
    });
    return _.uniq(_.compact(arr));
};

export const generateRandomString = (length: number): string => {
    return Math.random().toString(36).substr(2, length);
};

export const generateRandomStrings = (arrayLength: number, stringLength = 5): string[] =>
    _.times(arrayLength, () => generateRandomString(stringLength));

export const pluralize = (word: string, count: number): string => {
    return count === 1 ? word : `${word}s`;
};
