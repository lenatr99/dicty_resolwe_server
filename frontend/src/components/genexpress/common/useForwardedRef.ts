import { MutableRefObject, useEffect, useRef } from 'react';

/**
 * If ref of an underlying element has to be exposed, use ref from this hook to handle
 * cases when parent component doesn't need the ref..
 * @param ref ForwardedRef. If it's null, a new ref object is created.
 */
const useForwardedRef = <T>(
    ref: ((instance: T | null) => void) | React.MutableRefObject<T | null> | null | undefined,
): MutableRefObject<T | null> => {
    const innerRef = useRef<T>(null);
    useEffect(() => {
        if (!ref) return;
        if (typeof ref === 'function') {
            ref(innerRef.current);
        } else {
            ref.current = innerRef.current;
        }
    });

    return innerRef;
};

export default useForwardedRef;
