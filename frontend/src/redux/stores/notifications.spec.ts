import notificationsReducer, {
    NotificationsState,
    addSnackbar,
    removeSnackbar,
} from './notifications';
import { SnackbarNotificationContent } from 'redux/models/internal';
import { generateNotification } from 'tests/mock';

describe('notifications store', () => {
    let initialState: NotificationsState;

    beforeEach(() => {
        initialState = {
            notifications: [generateNotification(1)],
        };
    });

    it('should add a notification to state', () => {
        const newNotification: SnackbarNotificationContent = generateNotification(2);
        const newState = notificationsReducer(initialState, addSnackbar(newNotification));
        const expectedState = {
            notifications: [...initialState.notifications, newNotification],
        };

        expect(newState).toEqual(expectedState);
    });

    it('should remove notification from state', () => {
        const newState = notificationsReducer(
            initialState,
            removeSnackbar(initialState.notifications[0].key),
        );
        const expectedState = {
            notifications: [],
        };

        expect(newState).toEqual(expectedState);
    });
});
