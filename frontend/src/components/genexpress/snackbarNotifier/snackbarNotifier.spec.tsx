import React from 'react';
import SnackbarNotifier from './snackbarNotifier';
import { customRender } from 'tests/test-utils';
import { mockStore, testState } from 'tests/mock';
import { SnackbarNotifications } from 'redux/models/internal';

const initialTestState = testState();
const testNotifications: SnackbarNotifications = [
    { key: '1', message: 'Test message', variant: 'default' },
    { key: '2', message: 'Test message 2', variant: 'error' },
    { key: '3', message: 'Test message 3', variant: 'success' },
    { key: '4', message: 'Test message 4', variant: 'info' },
    { key: '5', message: 'Test message 5', variant: 'warning' },
];

initialTestState.notifications.notifications = testNotifications;

const mockedStore = mockStore(initialTestState);

describe('snackbarNotifier', () => {
    it('should render as snapshot with notifications in all different variants', () => {
        const { asFragment } = customRender(<SnackbarNotifier />, {
            mockedStore,
        });

        expect(asFragment()).toMatchSnapshot();
    });
});
