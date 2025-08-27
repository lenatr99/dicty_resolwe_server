import useResizeObserver, { UseResizeObserverCallback } from '@react-hook/resize-observer';
import { debounce, pick } from 'lodash';
import { RefObject, useEffect, useLayoutEffect, useRef, useState } from 'react';

const useSize = (targetRef: RefObject<HTMLElement>, refreshRate = 16) => {
    const [size, setSize] = useState<{ width: number; height: number }>({ width: 0, height: 0 });
    const isUnmounted = useRef(false);
    const debouncedSetSize = useRef(
        debounce<UseResizeObserverCallback>((entry) => {
            if (!isUnmounted.current) {
                setSize(entry.contentRect);
            }
        }, refreshRate),
    );

    useLayoutEffect(() => {
        if (targetRef.current != null) {
            setSize(pick(targetRef.current?.getBoundingClientRect(), ['width', 'height']));
        }
    }, [targetRef]);

    const resizeObserverRef = useResizeObserver(targetRef, debouncedSetSize.current);

    useEffect(() => {
        isUnmounted.current = false;
        const targetRefCurrent = targetRef.current;
        const debouncedSetSizeCurrent = debouncedSetSize.current;

        return () => {
            if (targetRefCurrent != null) {
                resizeObserverRef.unobserve(targetRefCurrent);
            }
            debouncedSetSizeCurrent.cancel();
            isUnmounted.current = true;
        };
    }, [resizeObserverRef, targetRef]);

    return size;
};

export default useSize;
