import { combineReducers, createSlice, PayloadAction } from '@reduxjs/toolkit';
import { SnackbarKey } from 'notistack';
import { SnackbarNotification, SnackbarNotificationContent } from 'redux/models/internal';

// State slices.
const notificationsInitialState = [] as SnackbarNotification[];
const notificationsSlice = createSlice({
    name: 'notifications',
    initialState: notificationsInitialState,
    reducers: {
        addSnackbar: (
            state,
            action: PayloadAction<SnackbarNotificationContent>,
        ): SnackbarNotification[] => {
            const { key } = action.payload;
            return [
                ...state,
                {
                    ...action.payload,
                    key: key || Math.ceil(Math.random() * 1000),
                },
            ];
        },
        removeSnackbar: (state, action: PayloadAction<SnackbarKey>): SnackbarNotification[] => {
            return state.filter((notification) => notification.key !== action.payload);
        },
    },
});

const notificationsReducer = combineReducers({
    notifications: notificationsSlice.reducer,
});

// Export actions.
export const { addSnackbar, removeSnackbar } = notificationsSlice.actions;
export const addErrorSnackbar = (
    message: string,
): PayloadAction<SnackbarNotificationContent, string> => addSnackbar({ message, variant: 'error' });

export const addInfoSnackbar = (
    message: string,
): PayloadAction<SnackbarNotificationContent, string> =>
    addSnackbar({ message, variant: 'success' });

export type NotificationsState = ReturnType<typeof notificationsReducer>;

export default notificationsReducer;

// Selectors (exposes the store to containers).
export const getNotifications = (state: NotificationsState): SnackbarNotification[] =>
    state.notifications;
