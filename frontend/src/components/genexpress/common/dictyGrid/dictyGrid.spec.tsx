import React from 'react';
import { screen, fireEvent, waitFor } from '@testing-library/react';
import _ from 'lodash';
import DictyGrid from './dictyGrid';
import { generateTimeSeriesById } from 'tests/mock';
import { customRender } from 'tests/test-utils';

const timeSeriesById = generateTimeSeriesById(2);
const timeSeries = _.flatMap(timeSeriesById);

describe('dictyGrid', () => {
    let container: HTMLElement;
    beforeEach(() => {
        ({ container } = customRender(
            <DictyGrid
                data={timeSeries}
                columnDefs={[
                    { field: 'id', headerName: 'Id', width: 20 },
                    { field: 'collection.name', headerName: 'Name', width: 50 },
                ]}
                getRowId={(data): string => data.id.toString(10)}
                filterLabel="Filter"
                selectedData={[timeSeries[0]]}
            />,
        ));
    });

    it('should display results', () => {
        // All timeSeries are displayed in ag-grid.
        timeSeries.forEach((ts) => expect(screen.getByText(ts.collection.name)));
    });

    it('should filter results', async () => {
        // Simulate filter change.
        fireEvent.change(screen.getByLabelText('Filter'), {
            target: { value: timeSeries[1].collection.name },
        });

        // Check if results are filtered (first time series isn't shown).
        await waitFor(() => {
            expect(screen.queryByText(timeSeries[0].collection.name)).toBeNull();
        });
        expect(screen.queryByText(timeSeries[0].collection.name)).toBeNull();
    });

    it('should mark already selected', async () => {
        await waitFor(() => {
            expect(container.querySelector("div[role='row'].ag-row-selected")).toBeInTheDocument();
        });
    });
});
