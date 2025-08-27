import { defineConfig } from 'vitest/config';
import tsconfigPaths from 'vite-tsconfig-paths';

export default defineConfig({
    plugins: [tsconfigPaths()],
    test: {
        globals: true,
        environment: 'jsdom',
        root: 'src',
        setupFiles: ['./src/setupTests.ts'],
        browser: {
            name: 'chrome',
            viewport: {
                width: 1920,
                height: 1080,
            },
        },
    },
});
