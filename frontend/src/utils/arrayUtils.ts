import _ from 'lodash';

/**
 * Joins an array into "value1, value2, and value3"
 */
export const advancedJoin = (
    values: (string | undefined)[],
    delimiter = ', ',
    lastDelimiter = ', and ',
): string => {
    const nonNullValues = _.compact(values);
    if (_.size(nonNullValues) === 0) {
        return '';
    }
    if (_.size(nonNullValues) === 1) {
        return nonNullValues[0];
    }
    return `${_.initial(nonNullValues).join(delimiter)}${lastDelimiter}${
        _.last(nonNullValues) as string
    }`;
};
