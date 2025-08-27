import { useEffect, useState } from 'react';

const useStateWithEffect = <T>(getValue: () => T, dependencies: unknown[]): T => {
    const [value, setValue] = useState<T>(getValue);

    useEffect(() => {
        setValue(getValue());
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, dependencies);

    return value;
};

export default useStateWithEffect;
