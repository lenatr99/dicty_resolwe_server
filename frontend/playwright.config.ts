import { PlaywrightTestConfig } from '@playwright/test';

const config: PlaywrightTestConfig = {
    use: {
        viewport: { width: 1600, height: 700 },
        ignoreHTTPSErrors: true,
        baseURL: 'http://localhost:3000',
        bypassCSP: true,
    },
    webServer:
        process.env.CI === 'true'
            ? {
                  command: 'yarn serve',
                  port: 3000,
              }
            : undefined,
};

export default config;
