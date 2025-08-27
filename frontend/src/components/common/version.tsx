import React, { ReactElement } from 'react';

export const Version = (): ReactElement => (
    <span className="version">{import.meta.env.VITE_APP_VERSION}</span>
);
