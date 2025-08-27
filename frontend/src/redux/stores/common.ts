import { ActionCreator, ActionReducerMapBuilder } from '@reduxjs/toolkit';

export const clearStateOnActions = <T>(
    builder: ActionReducerMapBuilder<T>,
    actions: ActionCreator<unknown>[],
    initialState: T,
): void => {
    actions.forEach((action) => {
        builder.addCase(action.toString(), () => initialState);
    });
};
