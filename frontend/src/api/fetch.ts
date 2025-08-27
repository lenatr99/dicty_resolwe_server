import { getCookie } from '../utils/documentHelpers';
import { ResponseError } from 'redux/models/internal';

export type QueryParams = { [key: string]: string | number | string[] | number[] };
type BodyParams = Record<string, unknown>;

/**
 * These methods do not require CSRF protection.
 * @param method - Method to check.
 * @returns True if the given method is safe (does not require CSRF cookie), false otherwise.
 */
const csrfSafeMethod = (method: string): boolean => {
    return /^(GET|HEAD|OPTIONS|TRACE)$/.test(method);
};

/**
 * The fetch() API only rejects a promise when a “network error is encountered, although this usually means permissions issues or similar.”.
 * Basically fetch() will only reject a promise if the user is offline, or some unlikely networking error occurs, such a DNS lookup failure.
 * That's why if response isn't ok, an error is thrown and it should be handled in epics / calling component.
 * @param response - Received response.
 */
const throwErrorIfResponseNotOk = (response: Response): Response => {
    if (!response.ok) {
        throw new ResponseError(response);
    }

    return response;
};

const request = async (
    url: string,
    params?: QueryParams | BodyParams,
    method: 'GET' | 'POST' | 'PUT' | 'DELETE' = 'GET',
): Promise<Response> => {
    const headers = new Headers({ 'Content-Type': 'application/json' });

    /* If method requires csrf protection, extra header X-CSRFToken is needed:
     *  https://docs.djangoproject.com/en/3.0/ref/csrf/#ajax
     */
    if (!csrfSafeMethod(method)) {
        headers.append('X-CSRFToken', getCookie('csrftoken'));
    }

    const options: RequestInit = {
        method,
        headers,
    };

    /* Development is usually done via proxy, that's why cookies aren't sent through (not same origin).
     *  That's why credentials are set to "include" in development mode.
     */
    if (import.meta.env.DEV) {
        options.credentials = 'include';
    }

    const fullUrl = new URL(url);
    // If params exists and method is GET, add query string to url
    // otherwise, just add params as a "body" property to the options object.
    if (params) {
        if (method === 'GET') {
            Object.keys(params).forEach((key) =>
                fullUrl.searchParams.append(key, (params as QueryParams)[key]?.toString()),
            );
        } else {
            options.body = JSON.stringify(params as BodyParams);
        }
    }

    const test = fetch(fullUrl.toString(), options);
    return test.then(throwErrorIfResponseNotOk);
};

export const get = (url: string, params?: QueryParams): Promise<Response> => {
    return request(url, params, 'GET');
};

export const post = (url: string, params?: BodyParams): Promise<Response> => {
    return request(url, params, 'POST');
};

export const put = (url: string, params?: BodyParams): Promise<Response> => {
    return request(url, params, 'PUT');
};

export const remove = (url: string, params?: BodyParams): Promise<Response> => {
    return request(url, params, 'DELETE');
};
