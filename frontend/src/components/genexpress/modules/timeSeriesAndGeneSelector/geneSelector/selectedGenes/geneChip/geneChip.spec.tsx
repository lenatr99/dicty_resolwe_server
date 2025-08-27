import React from 'react';
import { screen, fireEvent, waitFor } from '@testing-library/react';
import { GeneChip } from './geneChip';
import { customRender } from 'tests/test-utils';
import { generateGene } from 'tests/mock';

const gene = generateGene(0);

describe('geneChip', () => {
    it('should render the same as default snapshot', () => {
        const { asFragment } = customRender(
            <GeneChip
                gene={gene}
                highlighted={false}
                onHighlight={(): void => {}}
                onUnhighlight={(): void => {}}
                onRemove={(): void => {}}
            />,
        );

        expect(asFragment()).toMatchSnapshot();
    });

    it('should render as highlighted snapshot', () => {
        const { asFragment } = customRender(
            <GeneChip
                gene={gene}
                highlighted
                onHighlight={(): void => {}}
                onUnhighlight={(): void => {}}
                onRemove={(): void => {}}
            />,
        );

        expect(asFragment()).toMatchSnapshot();
    });

    it('should open gene information on click and close it on click away', async () => {
        customRender(
            <GeneChip
                gene={gene}
                highlighted
                onHighlight={(): void => {}}
                onUnhighlight={(): void => {}}
                onRemove={(): void => {}}
            />,
        );

        fireEvent.click(screen.getByText(gene.name));

        expect(await screen.findByText(gene.description));

        fireEvent.click(document);

        await waitFor(() => {
            expect(screen.queryByText(gene.description)).not.toBeInTheDocument();
        });
    });
});
