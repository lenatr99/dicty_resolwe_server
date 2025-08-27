import React from 'react';
import { screen, fireEvent } from '@testing-library/react';
import { MockStoreEnhanced } from 'redux-mock-store';
import { vi } from 'vitest';
import SelectedGenes from './selectedGenes';
import { customRender } from 'tests/test-utils';
import { generateBasketInfo, generateGene, mockStore, testState } from 'tests/mock';
import { allGenesDeselected } from 'redux/stores/genes';
import { RootState } from 'redux/rootReducer';
import { AppDispatch } from 'redux/appStore';

const initialState = testState();
const genes = [generateGene(0), generateGene(1)];

initialState.timeSeries.basketInfo = generateBasketInfo('1');
initialState.genes.byId = {
    [genes[0].name]: genes[0],
    [genes[1].name]: genes[1],
};

Object.assign(navigator, {
    clipboard: {
        writeText: (): void => {},
    },
});

describe('selectedGenes', () => {
    let mockedStore: MockStoreEnhanced<RootState, AppDispatch>;

    beforeEach(() => {
        mockedStore = mockStore(initialState);
    });

    it('should show genes', () => {
        customRender(<SelectedGenes selectedGenes={genes} highlightedGenesIds={[]} />, {
            mockedStore,
        });

        // All genes should be visible.
        genes.forEach((gene) => screen.getByText(gene.name));
    });

    it('should copy selected genes names', () => {
        customRender(<SelectedGenes selectedGenes={genes} highlightedGenesIds={[]} />, {
            mockedStore,
        });

        const clipboardSpy = vi.spyOn(navigator.clipboard, 'writeText');

        fireEvent.click(screen.getByLabelText(`Copy ${genes.length} genes to clipboard`));

        // All genes should be visible.
        expect(clipboardSpy).toHaveBeenCalledWith(genes.map((gene) => gene.name).join(', '));
    });

    it('should clear all genes', () => {
        customRender(<SelectedGenes selectedGenes={genes} highlightedGenesIds={[]} />, {
            mockedStore,
        });

        fireEvent.click(screen.getByLabelText('Clear all'));

        expect(mockedStore.getActions()).toEqual([allGenesDeselected()]);
    });
});
