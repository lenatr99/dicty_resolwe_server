import React from 'react';
import { screen, fireEvent, waitFor } from '@testing-library/react';
import { Mock, vi } from 'vitest';
import ManageGeneSetsModal from './manageGeneSetsModal';
import { generateGeneSets } from 'tests/mock';
import { customRender, waitForButtonEnabled } from 'tests/test-utils';

// GeneSets have to be sorted so that first element in array is also first row in grid.
// Needed to assert selected item with '~toggle Row Selection' (first row).
const testGeneSets = generateGeneSets(5).sort(
    (a, b) => b.dateTime.getTime() - a.dateTime.getTime(),
);

describe('manageGeneSets', () => {
    let mockedOnDelete: Mock;
    let mockedOnClick: Mock;
    let mockedOnClose: Mock;

    beforeEach(() => {
        mockedOnDelete = vi.fn();
        mockedOnClick = vi.fn();
        mockedOnClose = vi.fn();

        customRender(
            <ManageGeneSetsModal
                geneSets={testGeneSets}
                onDelete={mockedOnDelete}
                onClick={mockedOnClick}
                onClose={mockedOnClose}
            />,
        );
    });

    it('should render all gene sets in a grid', () => {
        testGeneSets.forEach(async (testGeneSet) => {
            expect(await screen.findByText(testGeneSet.genesNames[0], { exact: false }));
        });
    });

    it('should call onDelete with selected gene sets', async () => {
        fireEvent.click(screen.getByRole('button', { name: 'Edit' }));
        // Simulate select first row (gene set).
        fireEvent.click(
            screen.getAllByLabelText('Press Space to toggle row selection (unchecked)')[0],
        );
        await waitForButtonEnabled(() => screen.getByText('Delete selected'));
        fireEvent.click(await screen.findByText('Delete selected'));

        await waitFor(() => {
            expect(mockedOnDelete.mock.calls.length).toBe(1);
            expect(mockedOnDelete.mock.calls[0][0]).toEqual([testGeneSets[0]]);
        });
    });

    it('should call onClick with clicked gene set', async () => {
        // Simulate click on first gene set.
        fireEvent.click(screen.getByText(testGeneSets[0].genesNames.join(', ')));

        // Grid function 'onRowClicked' is async, that's why we must wait for it to execute.
        await waitFor(() => {
            expect(mockedOnClick.mock.calls.length).toBe(1);
        });
        expect(mockedOnClick.mock.calls[0][0]).toEqual(testGeneSets[0]);
    });

    it('should call onClose when user clicks close button', () => {
        fireEvent.click(screen.getByText('Close'));

        expect(mockedOnClose.mock.calls.length).toBe(1);
    });
});
