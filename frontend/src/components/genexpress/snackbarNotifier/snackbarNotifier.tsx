import { useEffect } from 'react';
import { connect, ConnectedProps } from 'react-redux';
import { useSnackbar, SnackbarKey } from 'notistack';
import { RootState } from 'redux/rootReducer';
import { getNotifications, removeSnackbar } from 'redux/stores/notifications';

let displayedKeys: SnackbarKey[] = [];

const mapStateToProps = (state: RootState) => {
    return {
        notifications: getNotifications(state.notifications),
    };
};

const connector = connect(mapStateToProps, {
    connectedRemoveSnackbar: removeSnackbar,
});

type PropsFromRedux = ConnectedProps<typeof connector>;

const SnackbarNotifier = ({ notifications, connectedRemoveSnackbar }: PropsFromRedux): null => {
    const { enqueueSnackbar } = useSnackbar();

    const storeDisplayed = (key: SnackbarKey): void => {
        displayedKeys = [...displayedKeys, key];
    };

    const removeDisplayed = (key: SnackbarKey): void => {
        displayedKeys = displayedKeys.filter((displayedKey) => displayedKey !== key);
    };

    useEffect(() => {
        notifications.forEach(({ key, message, variant, action }) => {
            // Do nothing if snackbar is already displayed.
            if (displayedKeys.includes(key)) return;

            // Display snackbar using notistack.
            enqueueSnackbar(message, {
                key,
                variant,
                action,
                onExited: (event, myKey) => {
                    // Remove this snackbar from redux store.
                    connectedRemoveSnackbar(myKey);
                    removeDisplayed(myKey);
                },
            });

            // Keep track of displayed snackbars.
            storeDisplayed(key);
        });
    }, [connectedRemoveSnackbar, enqueueSnackbar, notifications]);

    return null;
};

export default connector(SnackbarNotifier);
