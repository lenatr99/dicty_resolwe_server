import { sentryCapture } from './sentryUtils';
import { ResponseError } from 'redux/models/internal';
import { addErrorSnackbar } from 'redux/stores/notifications';

export const handleError = (
    message: string,
    associatedObject: Record<string, unknown> | Error | ResponseError = {},
    suppressSentryCapture = false,
): ReturnType<typeof addErrorSnackbar> => {
    if (!suppressSentryCapture) {
        sentryCapture(message, associatedObject, 'error');
    }
    return addErrorSnackbar(message);
};

export const logError = (
    message: string,
    associatedObject: Record<string, unknown> | Error = {},
): void => {
    sentryCapture(message, associatedObject, 'error');
};

export const logInfo = (message: string, associatedObject: Record<string, unknown>): void => {
    sentryCapture(message, associatedObject, 'info');
};

export const logWarning = (message: string, associatedObject: Record<string, unknown>): void => {
    sentryCapture(message, associatedObject, 'warning');
};
