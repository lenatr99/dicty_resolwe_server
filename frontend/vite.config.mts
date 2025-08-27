import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react-swc';
import svgr from 'vite-plugin-svgr';
import viteTsconfigPaths from 'vite-tsconfig-paths';
import type { PluginOption, ProxyOptions } from 'vite';

const targetDomain = 'localhost:8000';
const secure = false;

const wsProxyConfig: ProxyOptions = {
    // QA websocket is on a .ws subdomain, that's why the target is hardcoded.
    // Do not forget to change it when proxying to other environments.
    // target: `wss://${targetDomain}`,
    target: 'wss://ws.qa.genialis.io',
    ws: true,
    headers: {
        Host: 'ws.qa.genialis.io',
    },
    secure: false,
};

const proxyConfig: ProxyOptions = {
    target: `http://${targetDomain}`,  // Use HTTP for local Django backend
    secure,
    changeOrigin: true,
    rewrite: (path) => path.replace(/\/(\?|$)/, '$1'), // Remove trailing slash before query params or end of path
    headers: {
        Host: targetDomain,
        Origin: `http://${targetDomain}`,  // Use HTTP for local Django backend
        // Connection: 'keep-alive',
    },
};

const injectConfigScriptToIndex = (): PluginOption => {
    return {
        name: 'build-html',
        apply: 'build',
        transformIndexHtml: (html) => {
            return {
                html,
                tags: [
                    {
                        tag: 'script',
                        attrs: {
                            src: '/config.js?version=3.0',
                            charset: 'utf-8',
                        },
                        injectTo: 'head',
                    },
                ],
            };
        },
    };
};

// https://vitejs.dev/config/
export default defineConfig({
    plugins: [react(), svgr(), viteTsconfigPaths(), injectConfigScriptToIndex()],
    build: {
        sourcemap: true,
        outDir: 'build',
        // https://github.com/vitejs/vite/issues/15012
        rollupOptions: {
            onwarn(warning, defaultHandler) {
                if (warning.code === 'SOURCEMAP_ERROR') {
                    return;
                }

                defaultHandler(warning);
            },
        },
    },
    define: { global: 'window' },
    server: {
        proxy: {
            '/ws': wsProxyConfig,
            '/api': proxyConfig,
            '/saml-auth': proxyConfig,
        },
    },
});
