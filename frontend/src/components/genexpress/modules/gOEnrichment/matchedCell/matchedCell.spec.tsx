import React from 'react';
import { fireEvent, screen } from '@testing-library/react';
import { Mock, vi } from 'vitest';
import GOEnrichmentMatchedCell from './matchedCell';
import { customRender } from 'tests/test-utils';
import { generateGOEnrichmentRow } from 'tests/mock';
import { GOEnrichmentRow } from 'redux/models/internal';

describe('gOEnrichmentMatchedCell', () => {
    let goEnrichmentRow: GOEnrichmentRow;
    let asFragment: () => DocumentFragment;
    let mockedOnMatchedGenesClick: Mock;
    const value = 12;

    beforeEach(() => {
        mockedOnMatchedGenesClick = vi.fn();
        goEnrichmentRow = generateGOEnrichmentRow(1);

        ({ asFragment } = customRender(
            <GOEnrichmentMatchedCell
                value={value}
                data={goEnrichmentRow}
                onMatchedGenesClick={mockedOnMatchedGenesClick}
            />,
        ));
    });

    it('should render as default snapshot', () => {
        expect(asFragment()).toMatchSnapshot();
    });

    it('should call onMatchedGenesClick when user clicks on the link', () => {
        fireEvent.click(screen.getByText(`${value}/${goEnrichmentRow.total}`));

        expect(mockedOnMatchedGenesClick.mock.calls.length).toBe(1);
        expect(mockedOnMatchedGenesClick.mock.calls[0][0]).toBe(goEnrichmentRow);
    });
});
