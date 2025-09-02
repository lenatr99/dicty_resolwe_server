import { createSlice, PayloadAction, createSelector, combineReducers } from '@reduxjs/toolkit';
import { Relation } from '@genialis/resolwe/dist/api/types/rest';
import _ from 'lodash';

type SingleCellRelationsBySlug = Record<string, Relation>;

const bySlugInitialState: SingleCellRelationsBySlug = {};
const bySlugSlice = createSlice({
    name: 'singleCellSeriesBySlug',
    initialState: bySlugInitialState,
    reducers: {
        fetchSucceeded: (state, action: PayloadAction<Relation>): SingleCellRelationsBySlug => {
            return { ...state, [action.payload.slug]: action.payload };
        },
        reset: () => bySlugInitialState,
    },
});

const isFetchingSlice = createSlice({
    name: 'singleCellSeriesFetching',
    initialState: false,
    reducers: {
        started: () => true,
        ended: () => false,
    },
});

const singleCellSeriesReducer = combineReducers({
    bySlug: bySlugSlice.reducer,
    isFetching: isFetchingSlice.reducer,
});

export type SingleCellSeriesState = ReturnType<typeof singleCellSeriesReducer>;
export default singleCellSeriesReducer;

export const { fetchSucceeded: singleCellSeriesFetchSucceeded, reset: singleCellSeriesReset } =
    bySlugSlice.actions;
export const { started: singleCellSeriesFetchStarted, ended: singleCellSeriesFetchEnded } =
    isFetchingSlice.actions;

// Selectors
const getBySlug = (state: SingleCellSeriesState): SingleCellRelationsBySlug => state.bySlug;
export const getSingleCellRelationBySlug = (state: SingleCellSeriesState, slug: string): Relation | undefined =>
    getBySlug(state)[slug];

export const getSingleCellRelationForTimeSeriesSlug = createSelector(
    getBySlug,
    (_: SingleCellSeriesState, timeSeriesSlug?: string | null) => timeSeriesSlug,
    (bySlug, timeSeriesSlug) => {
        if (!timeSeriesSlug) return undefined;
        return bySlug[`${timeSeriesSlug}_sc`];
    },
);



