import { authApiUrl } from './base';

export const login = (nextPath: string) => {
    window.location.href = `${authApiUrl}/login/?next=${nextPath}`;
};

export const logout = (nextPath: string) => {
    window.location.href = `${authApiUrl}/logout/?next=${nextPath}`;
};
