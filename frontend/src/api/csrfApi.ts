import { apiUrl } from './base';
import { get } from './fetch';

const baseUrl = `${apiUrl}/base/csrf`;

export const getCSRFCookie = (): Promise<unknown> => {
    return get(baseUrl);
};
