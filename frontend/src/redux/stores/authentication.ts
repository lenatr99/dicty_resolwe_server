import { combineReducers, createSlice, PayloadAction } from '@reduxjs/toolkit';
import { User } from '@genialis/resolwe/dist/api/types/rest';
import _ from 'lodash';
import createIsFetchingSlice from './fetch';

const userInitialState = {} as User;
const userSlice = createSlice({
    name: 'authentication',
    initialState: userInitialState,
    reducers: {
        fetchSucceeded: (_state, action: PayloadAction<User | undefined>): User => {
            console.log('👤 USER AUTHENTICATION STATE CHANGE:');
            console.log('  User fetched:', action.payload ? `${action.payload.first_name} ${action.payload.last_name}` : 'Anonymous');
            console.log('  User ID:', action.payload?.id);
            console.log('─────────────────────────────────────────');
            
            return action.payload ?? userInitialState;
        },
    },
});

const isLoggedInInitialState = false;
const isLoggedInSlice = createSlice({
    name: 'authentication',
    initialState: isLoggedInInitialState,
    reducers: {},
    extraReducers: (builder) => {
        builder.addCase(userSlice.actions.fetchSucceeded, (_state, action): boolean => {
            return !_.isEmpty(action.payload);
        });
    },
});

const isLoggingInSlice = createIsFetchingSlice('authentication/isLoggingIn');
const isLoggingOutSlice = createIsFetchingSlice('authentication/isLoggingOut');
const isFetchingUserSlice = createIsFetchingSlice('authentication/isFetchingUser');

const authenticationReducer = combineReducers({
    user: userSlice.reducer,
    isLoggedIn: isLoggedInSlice.reducer,
    isLoggingIn: isLoggingInSlice.reducer,
    isLoggingOut: isLoggingOutSlice.reducer,
    isFetchingUser: isFetchingUserSlice.reducer,
});

export const { fetchSucceeded: userFetchSucceeded } = userSlice.actions;
export const { started: loginStarted, ended: loginEnded } = isLoggingInSlice.actions;
export const { started: logoutStarted, ended: logoutEnded } = isLoggingOutSlice.actions;
export const { started: userFetchStarted, ended: userFetchEnded } = isLoggingInSlice.actions;

export type AuthenticationState = ReturnType<typeof authenticationReducer>;

export default authenticationReducer;

export const getUser = (state: AuthenticationState): User => state.user;
export const getIsLoggedIn = (state: AuthenticationState): boolean => state.isLoggedIn;
export const getIsLoggingIn = (state: AuthenticationState): boolean => state.isLoggingIn;
export const getIsLoggingOut = (state: AuthenticationState): boolean => state.isLoggingOut;
export const getIsFetchingUser = (state: AuthenticationState): boolean => state.isFetchingUser;
