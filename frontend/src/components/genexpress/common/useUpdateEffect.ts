import { useEffect, useRef } from 'react';

/**
 * A custom useEffect hook that only triggers on updates, not on initial mount
 * Idea taken from: https://stackoverflow.com/a/55075818/1526448
 * @param effect - Function to execute on dependencies updates.
 * @param dependencies - Dependencies to rerun update.
 */
const useUpdateEffect = (effect: () => void, dependencies: unknown[]): void => {
    const isInitialMount = useRef(true);

    useEffect(() => {
        if (isInitialMount.current) {
            isInitialMount.current = false;
        } else {
            effect();
        }
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, dependencies);
};

export default useUpdateEffect;
