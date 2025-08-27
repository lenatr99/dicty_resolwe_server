import _ from 'lodash';
import { assertExhaustive } from './utils';

/**
 * Calculates logarithm of given base.
 */
export const logOfBase = (value: number, base: number): number => {
    if (value === null) return NaN;
    return Math.log(value) / Math.log(base);
};

/**
 * Arithmetic mean of the provided array. Forgiving with blank values (null, NaN, undefined).
 * simply ignores them both from the sum and the length.
 *
 * @param values An array containing sample data
 */
export const mean = (values: number[]): number => {
    const finiteArray = values.filter((value) => _.isFinite(value));
    return finiteArray.length > 0 ? _.sum(finiteArray) / finiteArray.length : NaN;
};

/**
 * Returns minimum and maximum value from array.
 * @param values An array with numberic values
 */
export const getMinMax = (values: number[]): [number, number] => {
    const min = Math.min(...values);
    const max = Math.max(...values);

    return [min, max];
};

export type NumberFormat = 'short' | 'long';
/**
 * Formats a number using fixed-point or exponential notation.
 */
export const formatNumber = (value: number, format: NumberFormat): string => {
    let precision = 0;

    if (format === 'short') {
        precision = 2;
    } else if (format === 'long') {
        precision = 5;
    } else {
        assertExhaustive({ format });
    }

    if (value === 0) return value.toFixed(precision);

    const boundary = 10 ** -precision;
    return -boundary < value && value < boundary
        ? value.toExponential(precision)
        : value.toFixed(precision);
};
