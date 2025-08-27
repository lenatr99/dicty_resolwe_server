import React from 'react';
import { screen } from '@testing-library/react';
import _ from 'lodash';
import ConnectedSingleCellSeriesSelector from './SingleCellSeriesSelector';
import { generateDescriptorSchema, testState } from 'tests/mock';
import {
    customRender,
    handleCommonRequests,
    resolveStringifiedObjectPromise,
    validateExportFile,
} from 'tests/test-utils';
import { RootState } from 'redux/rootReducer';
import { getSelectedTimeSeries } from 'redux/stores/timeSeries';
import { UrlDescriptor } from 'types/application';
import { DescriptorSchemaSlug } from 'components/genexpress/common/constants';

const dictyDescriptorSchema = {
    ...generateDescriptorSchema(1),
    schema: [
        { name: 'project', type: 'basic:string:', label: 'Project', required: false },
        { name: 'citation', type: 'basic:url:view:', label: 'Citation', required: false },
        { name: 'details', type: 'basic:string:', label: 'Details', required: false },
        { name: 'strain', type: 'basic:string:', label: 'Strain', required: false },
        { name: 'growth', type: 'basic:string:', label: 'Growth', required: false },
        { name: 'treatment', type: 'basic:string:', label: 'Treatment', required: false },
    ],
    slug: DescriptorSchemaSlug.DictyTimeSeries,
};

describe('SingleCellSeriesSelector', () => {
    let initialState: RootState;

    beforeEach(() => {
        fetchMock.resetMocks();

        fetchMock.mockResponse((req) => {
            if (req.url.includes('descriptorschema')) {
                return resolveStringifiedObjectPromise([dictyDescriptorSchema]);
            }

            return handleCommonRequests(req) ?? Promise.reject(new Error(`bad url: ${req.url}`));
        });

        initialState = testState();
        _.flatMap(initialState.timeSeries.byId)[1].descriptor = {};

        customRender(<ConnectedSingleCellSeriesSelector />, {
            initialState,
        });
    });

    it('should render all fields from descriptor schema', async () => {
        // Tests for displaying only collection name if descriptor schema is empty are in integration test files.
        if (_.flatMap(initialState.timeSeries.byId)[0].descriptor != null) {
            await screen.findByRole('button', {
                name: (
                    _.flatMap(initialState.timeSeries.byId)[0].descriptor.citation as UrlDescriptor
                ).name,
            });
        }
        await screen.findByText(_.flatMap(initialState.timeSeries.byId)[1].collection.name);
    });

    it('should export Collection/selectedCollection.tsv file', async () => {
        const selectedTimeSeries = getSelectedTimeSeries(initialState.timeSeries);
        await validateExportFile('Collection/selectedCollection.tsv', (exportFile) => {
            expect(exportFile?.content).toContain(selectedTimeSeries?.collection.name);
        });
    });
});
