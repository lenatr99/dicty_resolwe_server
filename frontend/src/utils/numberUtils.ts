import _ from 'lodash';

/**
 * Generates an array of random numbers. If "randomNumberGenerator" function isn't given,
 * Math.random() is used.
 * @param arrayLength - How many numbers are in generated array.
 * @param randomNumberGenerator - Number generator function.
 */
export const generateRandomNumbers = (
    arrayLength: number,
    randomNumberGenerator: () => number = (): number => Math.random(),
): number[] => _.times(arrayLength, () => randomNumberGenerator());
