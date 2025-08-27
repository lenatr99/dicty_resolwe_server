import { BackendAppState } from '../redux/models/rest';
import { deserializeResponse } from '../utils/apiUtils';
import { apiUrl } from './base';
import { get, post } from './fetch';
import { BookmarkReduxState } from 'redux/rootReducer';
import { BookmarkComponentsState } from 'redux/models/internal';

const baseUrl = `${apiUrl}/app-state`;

export const createAppState = async (
    appState: BackendAppState<BookmarkReduxState & BookmarkComponentsState>,
): Promise<BackendAppState<BookmarkReduxState & BookmarkComponentsState>> => {
    const createAppStateResponse = await post(`${baseUrl}`, appState);

    return deserializeResponse<BackendAppState<BookmarkReduxState & BookmarkComponentsState>>(
        createAppStateResponse,
    );
};

export const getAppState = async (
    appStateId: string,
): Promise<BookmarkReduxState & BookmarkComponentsState> => {
    const url = `${baseUrl}/${appStateId}`;

    const getAppStateJsonResponse = await get(url);

    return (
        await deserializeResponse<BackendAppState<BookmarkReduxState & BookmarkComponentsState>>(
            getAppStateJsonResponse,
        )
    ).state;
};
