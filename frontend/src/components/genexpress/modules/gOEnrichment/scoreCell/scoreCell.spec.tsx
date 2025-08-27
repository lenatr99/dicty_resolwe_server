import React from 'react';
import { fireEvent, screen } from '@testing-library/react';
import GOEnrichmentScoreCell from './scoreCell';
import { customRender } from 'tests/test-utils';
import { generateGOEnrichmentRow } from 'tests/mock';
import { GOEnrichmentRow } from 'redux/models/internal';
import { formatNumber } from 'utils/math';

describe('gOEnrichmentScoreCell', () => {
    let goEnrichmentRow: GOEnrichmentRow;
    let asFragment: () => DocumentFragment;
    const value = 12;

    beforeEach(() => {
        goEnrichmentRow = generateGOEnrichmentRow(1);

        ({ asFragment } = customRender(
            <GOEnrichmentScoreCell value={value} data={goEnrichmentRow} />,
        ));
    });

    it('should render as default snapshot', () => {
        expect(asFragment()).toMatchSnapshot();
    });

    it('should call onMatchedGenesClick when user clicks on the link', () => {
        fireEvent.click(screen.getByText(formatNumber(value, 'short')));
    });
});
