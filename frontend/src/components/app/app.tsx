import React, { ReactElement, useEffect } from 'react';
import { Route, Routes, BrowserRouter } from 'react-router-dom';
import './app.scss';
import '../../../node_modules/react-grid-layout/css/styles.css';
import '../../../node_modules/react-resizable/css/styles.css';
import styled, { ThemeProvider as StyledComponentsThemeProvider } from 'styled-components';
import { Provider as ReduxProvider } from 'react-redux';
import {
    ThemeProvider as MuiThemeProvider,
    createTheme,
    CssBaseline,
    StyledEngineProvider,
} from '@mui/material';
import { SnackbarProvider } from 'notistack';
import appStore from '../../redux/appStore';
import { GlobalStyle } from './globalStyle';
import theme from './theme';
import LandingPage from 'components/landing/landingPage';
import PageNotFound from 'components/pageNotFound';
import GeneExpressGrid from 'components/genexpress/geneExpressGrid';
import { RendererContext } from 'components/common/rendererContext';
import { getCSRFCookie } from 'api';
import CookieConsent from 'components/common/cookieConsent';

void getCSRFCookie();

const appTheme = createTheme(theme);

const StyledSnackbarProvider = styled(SnackbarProvider)`
    &.SnackbarItem-variantSuccess {
        background: ${({ theme: styledTheme }) => styledTheme.palette.secondary.main};
    }
`;

const App = (): ReactElement => {
    useEffect(() => {
        window.dataLayer = window.dataLayer || [];

        function gtag(..._args: unknown[]) {
            // eslint-disable-next-line prefer-rest-params
            window.dataLayer.push(arguments);
        }

        if (import.meta.env.PROD) {
            gtag('js', new Date());
            gtag('config', 'G-VLE4FY5V5F');
        }
    }, []);

    return (
        <RendererContext.Provider value="canvas">
            <StyledEngineProvider injectFirst>
                <ReduxProvider store={appStore}>
                    <MuiThemeProvider theme={appTheme}>
                        <CssBaseline />
                        <StyledComponentsThemeProvider theme={appTheme}>
                            <GlobalStyle />
                            <StyledSnackbarProvider maxSnack={3}>
                                <BrowserRouter>
                                    <Routes>
                                        <Route path="/" element={<LandingPage />} />
                                        <Route path="/landing" element={<LandingPage />} />
                                        <Route path="/bcm" element={<GeneExpressGrid />} />
                                        <Route path="*" element={<PageNotFound />} />
                                    </Routes>
                                </BrowserRouter>
                                <CookieConsent />
                            </StyledSnackbarProvider>
                        </StyledComponentsThemeProvider>
                    </MuiThemeProvider>
                </ReduxProvider>
            </StyledEngineProvider>
        </RendererContext.Provider>
    );
};

export default App;
