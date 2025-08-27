import React from 'react';
import { Mock, vi } from 'vitest';
import GOEnrichmentTermCell from './termCell';
import { customRender } from 'tests/test-utils';
import { generateGOEnrichmentRow } from 'tests/mock';
import { GOEnrichmentRow } from 'redux/models/internal';

describe('gOEnrichmentTermCell', () => {
    let goEnrichmentRow: GOEnrichmentRow;
    let asFragment: () => DocumentFragment;
    let mockedOnToggleCollapseClick: Mock;

    beforeEach(() => {
        mockedOnToggleCollapseClick = vi.fn();
        goEnrichmentRow = generateGOEnrichmentRow(1);
        goEnrichmentRow.term_name = 'Test term';
        goEnrichmentRow.depth = 3;

        ({ asFragment } = customRender(
            <GOEnrichmentTermCell
                data={goEnrichmentRow}
                onToggleCollapseClick={mockedOnToggleCollapseClick}
            />,
        ));
    });

    it('should render as default snapshot', () => {
        expect(asFragment()).toMatchSnapshot();
    });
});
