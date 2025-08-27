import React from 'react';
import { screen, fireEvent } from '@testing-library/react';
import { vi } from 'vitest';
import { GeneInformation } from './geneInformation';
import { customRender } from 'tests/test-utils';
import { generateGene } from 'tests/mock';

const gene = generateGene(0);

describe('geneInformation', () => {
    it('should render as default snapshot', () => {
        const { asFragment } = customRender(
            <GeneInformation
                gene={gene}
                highlighted={false}
                onHighlight={(): void => {}}
                onUnhighlight={(): void => {}}
            />,
        );

        expect(asFragment()).toMatchSnapshot();
    });

    it('should show Unhighlight button if gene is highlighted', () => {
        customRender(
            <GeneInformation
                gene={gene}
                highlighted
                onHighlight={(): void => {}}
                onUnhighlight={(): void => {}}
            />,
        );

        screen.getByText('Unhighlight');
    });

    it('should show Unhighlight button if gene is highlighted', () => {
        customRender(
            <GeneInformation
                gene={gene}
                highlighted={false}
                onHighlight={(): void => {}}
                onUnhighlight={(): void => {}}
            />,
        );

        screen.getByText('Highlight');
    });

    it('should call onHighlight when user clicks Highlight', () => {
        const mockedOnHighlight = vi.fn();
        customRender(
            <GeneInformation
                gene={gene}
                highlighted={false}
                onHighlight={mockedOnHighlight}
                onUnhighlight={(): void => {}}
            />,
        );

        fireEvent.click(screen.getByText('Highlight'));

        expect(mockedOnHighlight).toBeCalled();
    });

    it('should call onUnhighlight when user clicks Unhighlight', () => {
        const mockedOnUnhighlight = vi.fn();
        customRender(
            <GeneInformation
                gene={gene}
                highlighted
                onHighlight={(): void => {}}
                onUnhighlight={mockedOnUnhighlight}
            />,
        );

        fireEvent.click(screen.getByText('Unhighlight'));

        expect(mockedOnUnhighlight).toBeCalled();
    });
});
