import React from 'react';
import { screen, fireEvent, waitFor, RenderResult, configure } from '@testing-library/react';
import _ from 'lodash';
import {
    DataGafAnnotation,
    DONE_DATA_STATUS,
    WAITING_DATA_STATUS,
} from '@genialis/resolwe/dist/api/types/rest';
import { v4 as uuidv4 } from 'uuid';
import { Client, Server } from 'mock-socket';
import { aspectOptions } from './gOEnrichment';
import GeneExpressGrid from 'components/genexpress/geneExpressGrid';
import {
    customRender,
    getFetchMockCallsWithUrl,
    handleCommonRequests,
    resolveStringifiedObjectPromise,
    validateCreateStateRequest,
} from 'tests/test-utils';
import {
    testState,
    generateGenesById,
    generateGeneOntologyStorageJson,
    generateGaf,
    generateData,
    generateBackendBookmark,
    generateBookmarkQueryParameter,
} from 'tests/mock';
import { RootState } from 'redux/rootReducer';
import { sessionId, webSocketUrl } from 'api/base';
import { appendMissingAttributesToJson } from 'utils/gOEnrichmentUtils';
import { ProcessSlug, BookmarkStatePath } from 'components/genexpress/common/constants';
import { pValueThresholdsOptions } from 'redux/stores/gOEnrichment';
import { gOEnrichmentProcessDebounceTime } from 'redux/epics/gOEnrichmentEpics';

const genesById = generateGenesById(2);
const genes = _.flatMap(genesById);
const gOEnrichmentJson = generateGeneOntologyStorageJson(genes.map((gene) => gene.feature_id));
appendMissingAttributesToJson(gOEnrichmentJson, genes[0].source, genes[0].species);
const differentGOEnrichmentJson = generateGeneOntologyStorageJson(
    genes.map((gene) => gene.feature_id),
);
appendMissingAttributesToJson(differentGOEnrichmentJson, genes[0].source, genes[0].species);
const { humanGaf, mouseMGIGaf, mouseUCSCGaf } = generateGaf(1);

const backendBookmark = generateBackendBookmark(undefined, [genes[0].feature_id]);

backendBookmark.state.gOEnrichment.pValueThreshold = pValueThresholdsOptions[2];
backendBookmark.state.GOEnrichment = { selectedAspect: aspectOptions[1] };

const ontologyObo = generateData(123);

const dataObjectId = 123;

// GO Enrichment process is triggered after 3s debounce time.
configure({ asyncUtilTimeout: gOEnrichmentProcessDebounceTime + 1000 });

describe('goEnrichment integration', () => {
    let initialState: RootState;

    describe('process data exists', () => {
        beforeAll(() => {
            fetchMock.resetMocks();

            fetchMock.mockResponse((req) => {
                if (req.url.includes('get_or_create')) {
                    return resolveStringifiedObjectPromise({
                        id: dataObjectId,
                    });
                }

                if (decodeURIComponent(req.url).includes('data:gaf')) {
                    return resolveStringifiedObjectPromise([humanGaf, mouseMGIGaf, mouseUCSCGaf]);
                }

                if (req.url.includes('data') && req.url.includes('slug')) {
                    return resolveStringifiedObjectPromise(ontologyObo);
                }

                if (req.url.includes('data') && req.url.includes(dataObjectId.toString())) {
                    return resolveStringifiedObjectPromise([
                        {
                            ...generateData(1),
                            ...{
                                status: DONE_DATA_STATUS,
                                output: {
                                    terms: 1,
                                    species: genes[0].species,
                                    source: genes[0].source,
                                },
                            },
                        },
                    ]);
                }

                if (req.url.includes('storage')) {
                    return resolveStringifiedObjectPromise({
                        json: gOEnrichmentJson,
                    });
                }

                if (req.url.includes('feature')) {
                    return resolveStringifiedObjectPromise({
                        results: [genes[0], genes[1]],
                    });
                }

                if (req.url.includes('list_by_ids')) {
                    return resolveStringifiedObjectPromise(
                        backendBookmark.state.genes.selectedGenesIds.map(
                            (geneId) => genesById[geneId],
                        ),
                    );
                }

                if (req.url.includes('app-state')) {
                    return resolveStringifiedObjectPromise(backendBookmark);
                }

                return (
                    handleCommonRequests(req, genes) ??
                    Promise.reject(new Error(`bad url: ${req.url}`))
                );
            });
        });

        describe('gaf not fetched yet', () => {
            let unmount: RenderResult['unmount'];

            beforeEach(() => {
                initialState = testState();
                initialState.gOEnrichment.gaf = {} as DataGafAnnotation;

                ({ unmount } = customRender(<GeneExpressGrid />, {
                    initialState,
                }));
            });

            it('should show enrichment terms in a data grid after genes are chosen', async () => {
                fireEvent.change(screen.getByPlaceholderText('Search for a gene'), {
                    target: { value: genes[0].name.slice(0, 2) },
                });

                fireEvent.click(await screen.findByText(genes[0].name));

                await waitFor(() => {
                    screen.getByText(gOEnrichmentJson.tree.BP[0].term_name);
                });
            });

            it('should should enrichment terms if bookmark is loaded', async () => {
                unmount();

                customRender(<GeneExpressGrid />, {
                    initialState,
                    route: generateBookmarkQueryParameter(),
                });

                await screen.findByText(aspectOptions[1].label);
                await screen.findByText(pValueThresholdsOptions[2].toString());
                await screen.findByText(gOEnrichmentJson.tree[aspectOptions[1].value][0].term_name);
            });
        });
        describe('gaf already in store', () => {
            beforeEach(() => {
                initialState = testState();
                initialState.gOEnrichment.gaf = humanGaf;
            });

            it('should show enrichment terms in a data grid after genes are chosen', async () => {
                customRender(<GeneExpressGrid />, {
                    initialState,
                });

                fireEvent.change(screen.getByPlaceholderText('Search for a gene'), {
                    target: { value: genes[0].name.slice(0, 2) },
                });

                fireEvent.click(await screen.findByText(genes[0].name));

                await waitFor(() => {
                    screen.getByText(gOEnrichmentJson.tree.BP[0].term_name);
                });
            });
        });

        describe('gaf and data already in store', () => {
            beforeEach(() => {
                fetchMock.mockClear();

                initialState = testState();
                initialState.gOEnrichment.gaf = humanGaf;
                initialState.genes.byId = genesById;
                initialState.genes.selectedGenesIds = [genes[0].feature_id];
                initialState.gOEnrichment.json = differentGOEnrichmentJson;
                initialState.gOEnrichment.ontologyObo = ontologyObo;

                customRender(<GeneExpressGrid />, {
                    initialState,
                });
            });

            it('should show different enrichment terms after user changes p-value', async () => {
                // Click on dropdown. MouseDown event has to be used, because material-ui Select component
                // listens to mouseDown event to expand options menu.
                await waitFor(() => expect(screen.getByLabelText('p-value')).toBeEnabled());
                fireEvent.mouseDown(screen.getByLabelText('p-value'));
                fireEvent.click(await screen.findByText(pValueThresholdsOptions[2].toString()));

                await screen.findByText(gOEnrichmentJson.tree.BP[0].term_name);
            });

            it('should show different enrichment terms in a data grid user changes aspect', async () => {
                // Loop through all remaining aspect options and check if data grid is refreshing accordingly.
                for (let i = 1; i < aspectOptions.length; i += 1) {
                    const aspectOption = aspectOptions[i];

                    // Click on dropdown. MouseDown event has to be used, because material-ui Select component
                    // listens to mouseDown event to expand options menu.
                    fireEvent.mouseDown(await screen.findByLabelText('Aspect'));
                    fireEvent.click(await screen.findByText(aspectOption.label));

                    await screen.findByText(
                        differentGOEnrichmentJson.tree[aspectOption.value][0].term_name,
                    );
                }
            });

            it('should save selected time series, genes, highlighted genes and all component bookmarkable state to app-state api', async () => {
                await waitFor(() => {
                    expect(screen.getByLabelText('Bookmark')).toBeEnabled();
                });
                fireEvent.click(screen.getByLabelText('Bookmark'));

                await validateCreateStateRequest((bookmarkState) => {
                    expect(bookmarkState.gOEnrichment.pValueThreshold).toEqual(
                        initialState.gOEnrichment.pValueThreshold,
                    );
                    expect(
                        _.get(bookmarkState, BookmarkStatePath.gOEnrichmentSelectedAspect),
                    ).toEqual(aspectOptions[0]);
                });
            });
        });
    });

    describe("process data doesn't exist", () => {
        const observerId = uuidv4();
        let webSocketMock: Client;

        beforeAll(() => {
            fetchMock.resetMocks();

            fetchMock.mockResponse(async (req) => {
                if (
                    req.url.includes('get_or_create') &&
                    (await req.json()).process.slug === ProcessSlug.goEnrichment
                ) {
                    return resolveStringifiedObjectPromise({
                        id: dataObjectId,
                    });
                }

                if (req.url.includes('data:gaf')) {
                    return resolveStringifiedObjectPromise([humanGaf, mouseMGIGaf, mouseUCSCGaf]);
                }

                if (req.url.includes('subscribe')) {
                    return resolveStringifiedObjectPromise({
                        subscription_id: observerId,
                    });
                }

                if (req.url.includes('data') && req.url.includes('slug')) {
                    return resolveStringifiedObjectPromise(ontologyObo);
                }

                if (req.url.includes('data') && req.url.includes(dataObjectId.toString())) {
                    return resolveStringifiedObjectPromise([
                        {
                            ...generateData(1),
                            ...{
                                status: WAITING_DATA_STATUS,
                                output: {},
                            },
                        },
                    ]);
                }

                if (req.url.includes('storage')) {
                    return resolveStringifiedObjectPromise({
                        json: gOEnrichmentJson,
                    });
                }

                if (req.url.includes('feature')) {
                    return resolveStringifiedObjectPromise({
                        results: [genes[0], genes[1]],
                    });
                }

                return (
                    handleCommonRequests(req) ?? Promise.reject(new Error(`bad url: ${req.url}`))
                );
            });
        });

        beforeEach(() => {
            initialState = testState();
            initialState.gOEnrichment.gaf = humanGaf;

            const mockServer = new Server(`${webSocketUrl}/${sessionId}`);
            mockServer.on('connection', (socket) => {
                webSocketMock = socket;
            });

            customRender(<GeneExpressGrid />, {
                initialState,
            });
        });

        it('should display gOEnrichment data fetched via websocket, after gene is chosen', async () => {
            fireEvent.change(screen.getByPlaceholderText('Search for a gene'), {
                target: { value: genes[0].name.slice(0, 2) },
            });

            fireEvent.click(await screen.findByText(genes[0].name));

            // Wait for data object with 'waiting' status is returned.
            await waitFor(() => {
                expect(getFetchMockCallsWithUrl(`api/data?id=${dataObjectId}`)).toHaveLength(1);
            });

            await screen.findByTestId('ScheduleIcon');

            fetchMock.doMockOnce(() =>
                resolveStringifiedObjectPromise([
                    {
                        ...generateData(1),
                        ...{
                            status: DONE_DATA_STATUS,
                            output: {
                                terms: 1,
                                species: genes[0].species,
                                source: genes[0].source,
                            },
                        },
                    },
                ]),
            );

            webSocketMock.send(
                JSON.stringify({
                    change_type: 'UPDATE',
                    subscription_id: observerId,
                    object_id: 'id',
                }),
            );

            await screen.findByText(gOEnrichmentJson.tree.BP[0].term_name);
        });
    });
});
