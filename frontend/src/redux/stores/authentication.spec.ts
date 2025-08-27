import { User } from '@genialis/resolwe/dist/api/types/rest';
import authenticationReducer, { AuthenticationState, userFetchSucceeded } from './authentication';
import { generateUser } from 'tests/mock';

const user = generateUser(1);

describe('authentication store', () => {
    let initialState: AuthenticationState;

    describe('not logged in', () => {
        beforeEach(() => {
            initialState = {
                isFetchingUser: false,
                isLoggedIn: false,
                isLoggingIn: false,
                isLoggingOut: false,
                user: {} as User,
            };
        });

        it('should add fetched user to state with fetchUserSucceeded action', () => {
            const newState = authenticationReducer(initialState, userFetchSucceeded(user));
            const expectedState = {
                ...initialState,
                isLoggedIn: true,
                user,
            };

            expect(newState).toEqual(expectedState);
        });
    });

    describe('logged in', () => {
        beforeEach(() => {
            initialState = {
                isFetchingUser: false,
                isLoggedIn: true,
                isLoggingIn: false,
                isLoggingOut: false,
                user,
            };
        });

        it('should clear user from state fetchUserSucceeded(null) action', () => {
            const newState = authenticationReducer(initialState, userFetchSucceeded(undefined));
            const expectedState = {
                ...initialState,
                isLoggedIn: false,
                user: {},
            };

            expect(newState).toEqual(expectedState);
        });
    });
});
