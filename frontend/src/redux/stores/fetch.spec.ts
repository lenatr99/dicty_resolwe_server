import createIsFetchingSlice from './fetch';

const fetchingSlice = createIsFetchingSlice('test');

describe('fetch store', () => {
    it('should set isFetching to false', () => {
        const newState = fetchingSlice.reducer(true, fetchingSlice.actions.ended());
        expect(newState).toBeFalsy();
    });

    it('should set isFetching to true', () => {
        const newState = fetchingSlice.reducer(true, fetchingSlice.actions.started());
        expect(newState).toBeTruthy();
    });
});
