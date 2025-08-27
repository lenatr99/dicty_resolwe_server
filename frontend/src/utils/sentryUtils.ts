import { captureException, captureMessage, init, setUser, SeverityLevel } from '@sentry/browser';
import { User } from '@genialis/resolwe/dist/api/types/rest';
import _ from 'lodash';
import { getUsername } from './user';
import { ResponseError } from 'redux/models/internal';

export const initializeSentry = (sentryUrl: string): void => {
    init({
        dsn: sentryUrl,
        // Send stack trace for non-error messages too.
        attachStacktrace: true,
    });
};

export const setSentryUser = (user: User | undefined): void => {
    if (user != null) {
        setUser({
            id: `${user.id}`,
            username: user.username,
            name: getUsername(user),
        });
    } else {
        setUser(null);
    }
};

export const sentryCapture = (
    message: string,
    associatedObject: Record<string, unknown> | Error | ResponseError = {},
    severity: SeverityLevel = 'error',
) => {
    if (severity === 'error' || _.isError(associatedObject)) {
        if (_.isError(associatedObject)) {
            associatedObject.message = `${message} ${associatedObject.message}`;
        }

        captureException(associatedObject, {
            level: severity,
            extra: {
                message,
                responseData: associatedObject,
            },
        });
    } else {
        captureMessage(message, {
            extra: {
                associatedObject,
            },
            level: severity,
        });
    }
};
