// eslint-disable-next-line spaced-comment
/// <reference types="vite/client" />
// eslint-disable-next-line spaced-comment
/// <reference types="vite-plugin-svgr/client" />

interface ImportMetaEnv {
    readonly VITE_APP_VERSION: string;
    readonly VITE_APP_SENTRY_URL: string;
}

interface ImportMeta {
    readonly env: ImportMetaEnv;
}

declare module '*.png' {
    const src: string;
    export default src;
}

declare module '*.svg' {
    const src: string;
    export default src;
}
