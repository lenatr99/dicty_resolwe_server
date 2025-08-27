declare namespace NodeJS {
    interface Global {
        baseURL: string;
        port: number;
    }
}
interface Process {
    env: {
        PORT?: number;
        LOGIN_EMAIL?: string;
        LOGIN_PASSWORD?: string;
    };
}
