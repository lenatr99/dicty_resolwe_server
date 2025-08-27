import { useEffect } from 'react';

export const useBrowserVisibility = ({
    onShow,
    onHide,
}: {
    onShow?: () => void;
    onHide?: () => void;
}) => {
    useEffect(() => {
        const handleVisibilityChange = () => {
            const visibilityState = document.visibilityState;
            if (visibilityState === 'visible') {
                onShow?.();
            } else if (visibilityState === 'hidden') {
                onHide?.();
            }
        };

        document.addEventListener('visibilitychange', handleVisibilityChange);

        return () => {
            document.removeEventListener('visibilitychange', handleVisibilityChange);
        };
    }, [onShow, onHide]);
};

export default useBrowserVisibility;
