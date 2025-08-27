import { useTheme } from '@mui/material';
import React, { ReactElement } from 'react';
import ReactCookieConsent from 'react-cookie-consent';

const CookieConsent = (): ReactElement => {
    const theme = useTheme();

    return (
        <ReactCookieConsent
            buttonStyle={{
                background: theme.palette.secondary.main,
                color: 'white',
                borderRadius: theme.shape.borderRadius,
            }}
        >
            This website uses cookies to ensure you get the best experience on our website.{' '}
            <a
                href="https://www.genialis.com/privacy-policy/"
                rel="noopener noreferrer"
                target="_blank"
            >
                Learn more
            </a>
        </ReactCookieConsent>
    );
};

export default CookieConsent;
