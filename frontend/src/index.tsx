import { createRoot } from 'react-dom/client';
import React, { StrictMode } from 'react';
import { createBrowserHistory } from 'history';
import { initializeSentry } from 'utils/sentryUtils';
import App from 'components/app/app';

const history = createBrowserHistory();

if (import.meta.env.PROD) {
    if (import.meta.env.VITE_APP_SENTRY_URL) {
        initializeSentry(import.meta.env.VITE_APP_SENTRY_URL);
    }
}
const rootElement = document.getElementById('root');

if (!rootElement) {
    throw new Error('container not found!');
}
const root = createRoot(rootElement);

/*
 * S3 bucket has limitations for client side routing. A special rule was
 * added to remove 404 / redirect error.
 * Based on https://via.studio/journal/hosting-a-reactjs-app-with-routing-on-aws-s3
 *
 * Basically S3 bucket redirects all not found paths to a path with a "#!" prefix.
 * Application has to remove this special character and fix URL in browser.
 */

const path = (/#!(.*)$/.exec(location.hash) || [])[1];
if (path) {
    history.replace(path);
}

root.render(
    <StrictMode>
        <App />
    </StrictMode>,
);

if (import.meta.hot) {
    import.meta.hot.accept('./components/app/app', () => {
        // eslint-disable-next-line @typescript-eslint/no-require-imports
        const NextApp = require('./components/app/app').default;
        root.render(<NextApp />);
    });
}
