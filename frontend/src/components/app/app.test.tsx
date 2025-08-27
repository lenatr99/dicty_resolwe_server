import React from 'react';
import { render } from '@testing-library/react';
import App from './app';

describe('app', () => {
    it('app renders', async () => {
        const { findByText } = render(<App />);
        expect(await findByText('dictyExpress')).toBeInTheDocument();
    });
});
