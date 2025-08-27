import { createSlice, Slice } from '@reduxjs/toolkit';

const createIsFetchingSlice = (
    sliceName: string,
): Slice<
    boolean,
    {
        started: () => boolean;
        ended: () => boolean;
    },
    string
> => {
    return createSlice({
        name: `fetch/${sliceName}`,
        initialState: false as boolean,
        reducers: {
            started: (): boolean => true,
            ended: (): boolean => false,
        },
    });
};

export default createIsFetchingSlice;
