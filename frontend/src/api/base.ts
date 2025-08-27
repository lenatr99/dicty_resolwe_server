import { v4 as uuidv4 } from 'uuid';

export const baseUrl = window.location.origin;
export const apiUrl = API_URL ?? `${baseUrl}/api`;
export const authApiUrl = SAML_AUTH_URL ?? `${baseUrl}/saml-auth`;

const websocketHost =
    WEBSOCKET_URL ??
    `ws${window.location.hostname === 'localhost' ? '' : 's'}://${window.location.host}/ws/v2`;
export const webSocketUrl = import.meta.env.VITE_APP_TEST_WSS ?? websocketHost;

export const sessionId = uuidv4();
