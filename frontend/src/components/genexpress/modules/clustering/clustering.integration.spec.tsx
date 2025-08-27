import { screen, fireEvent, waitFor } from '@testing-library/react';
import _ from 'lodash';
import { vi } from 'vitest';
import { distanceMeasureOptions, linkageFunctionOptions } from './clustering';
import { highlightedColor } from './clusteringChart';
import GeneExpressGrid from 'components/genexpress/geneExpressGrid';
import {
    customRender,
    handleCommonRequests,
    hoverOverVegaSymbol,
    resolveStringifiedObjectPromise,
    validateCreateStateRequest,
    validateExportFile,
} from 'tests/test-utils';
import {
    testState,
    generateSamplesExpressionsById,
    generateTimeSeriesById,
    generateRelationPartitions,
    generateBasketExpression,
    generateBookmarkQueryParameter,
    generateBackendBookmark,
    generateGenesById,
} from 'tests/mock';
import * as reportBuilder from 'components/genexpress/common/reportBuilder/reportBuilder';
import { RootState } from 'redux/rootReducer';
import { Gene } from 'redux/models/internal';
import {
    BookmarkStatePath,
    ClusteringLinkageFunction,
    DistanceMeasure,
} from 'components/genexpress/common/constants';

const timeSeriesById = generateTimeSeriesById(1);
const selectedTimeSeriesId = parseInt(_.keys(timeSeriesById)[0], 10);

const basketExpressions = [generateBasketExpression(), generateBasketExpression()];

const genesById = generateGenesById(10);
const genes = _.flatMap(genesById);

const samplesExpressionsById = generateSamplesExpressionsById(
    5,
    genes.map((gene) => gene.feature_id),
);
const samplesExpressionsIds = _.map(_.keys(samplesExpressionsById), Number);
timeSeriesById[selectedTimeSeriesId].partitions = generateRelationPartitions(samplesExpressionsIds);

const backendBookmark = generateBackendBookmark(selectedTimeSeriesId, [
    genes[0].feature_id,
    genes[1].feature_id,
]);
backendBookmark.state.Clustering = {
    distanceMeasure: distanceMeasureOptions[2].value,
    linkageFunction: linkageFunctionOptions[2].value,
};

vi.setConfig({ testTimeout: 10000 });
describe('clustering integration', () => {
    let initialState: RootState;
    let container: HTMLElement;

    const validateChart = async (genesToValidate: Gene[]): Promise<void> =>
        waitFor(
            () => {
                expect(
                    container.querySelector("g[role='graphics-symbol'].horizontalLines > path"),
                ).toBeInTheDocument();
                expect(
                    container.querySelector("g[role='graphics-symbol'].verticalLines > path"),
                ).toBeInTheDocument();
                expect(
                    container.querySelector("g[role='graphics-object'].genesNames > text"),
                ).toBeInTheDocument();

                genesToValidate.forEach((gene) => {
                    expect(
                        Array.from(
                            container.querySelectorAll(
                                "g[role='graphics-object'].genesNames > text",
                            ),
                        ).filter((element) => element.textContent === gene.name),
                    ).toHaveLength(1);
                });
            },
            { timeout: 7000 },
        );

    beforeEach(() => {
        fetchMock.mockClear();

        /* InitialState must have all the data needed for gene expressions (used in
         * clustering heatmap). Fetching / loading geneExpressions data is tested in
         * genesExpressions module.
         */
        initialState = testState();
        initialState.timeSeries.byId = timeSeriesById;
        initialState.timeSeries.basketExpressionsIds = basketExpressions.map(
            (basketExpression) => basketExpression.id,
        );
        initialState.timeSeries.selectedId = selectedTimeSeriesId;
        initialState.genes.byId = genesById;
        initialState.samplesExpressions.byId = samplesExpressionsById;
    });

    beforeAll(() => {
        fetchMock.resetMocks();

        fetchMock.mockResponse((req) => {
            if (req.url.includes('add_samples')) {
                return resolveStringifiedObjectPromise({
                    id: 1,
                    permitted_organisms: ['Organism'],
                    permitted_sources: ['Source'],
                });
            }

            if (req.url.includes('basket_expressions')) {
                return resolveStringifiedObjectPromise(basketExpressions);
            }

            if (req.url.includes('feature')) {
                return resolveStringifiedObjectPromise({
                    results: [genes[0], genes[1]],
                });
            }
            if (req.url.includes('data')) {
                return resolveStringifiedObjectPromise(
                    [...samplesExpressionsIds].map((sampleExpressionId) => ({
                        id: sampleExpressionId,
                        output: {
                            exp_json: 123,
                        },
                        entity: {
                            id: sampleExpressionId,
                        },
                    })),
                );
            }

            if (req.url.includes('app-state')) {
                return resolveStringifiedObjectPromise(backendBookmark);
            }

            return handleCommonRequests(req) ?? Promise.reject(new Error(`bad url: ${req.url}`));
        });
    });

    describe('no genes selected', () => {
        beforeEach(() => {
            initialState.genes.selectedGenesIds = [];

            ({ container } = customRender(<GeneExpressGrid />, {
                initialState,
            }));
        });

        it('input controls should be disabled', () => {
            expect(screen.getByLabelText('Distance Measure')).toHaveClass('Mui-disabled');
            expect(screen.getByLabelText('Clustering Linkage')).toHaveClass('Mui-disabled');
        });

        it('should show a warning that no genes are selected', async () => {
            await screen.findByText('Select two or more genes.');
        });

        it('should show a warning that only one gene is selected (after it is selected)', async () => {
            fireEvent.change(screen.getByPlaceholderText('Search for a gene'), {
                target: { value: genes[0].name.slice(0, 2) },
            });

            fireEvent.click(await screen.findByText(genes[0].name));
            await screen.findByText(
                'Correlation between samples cannot be computed on a single gene.',
                { exact: false },
            );
        });

        it('should show hierarchical clustering plot after two genes are chosen', async () => {
            fireEvent.change(screen.getByPlaceholderText('Search for a gene'), {
                target: { value: genes[0].name.slice(0, 2) },
            });

            fireEvent.click(await screen.findByText(genes[0].name));

            fireEvent.change(screen.getByPlaceholderText('Search for a gene'), {
                target: { value: genes[1].name.slice(0, 3) },
            });

            fireEvent.click(await screen.findByText(genes[1].name));

            await validateChart([genes[0], genes[1]]);
        });

        it('should not export anything', async () => {
            const files = await reportBuilder.getRegisteredComponentsExportFiles();
            expect(
                files.filter((file) => file.path.startsWith('Sample Hierarchical Clustering')),
            ).toHaveLength(0);
        });
    });

    describe('two genes selected', () => {
        beforeEach(async () => {
            initialState.timeSeries.basketExpressionsIds = basketExpressions.map(
                (basketExpression) => basketExpression.id,
            );
            initialState.genes.selectedGenesIds = [genes[0].feature_id, genes[1].feature_id];

            ({ container } = customRender(<GeneExpressGrid />, {
                initialState,
            }));

            await waitFor(() => {
                expect(
                    container.querySelector("g[role='graphics-symbol'].genesExpressionsHeatmap"),
                ).toBeInTheDocument();
            });
        });

        it('should show different hierarchical clustering after user changes distance measure', async () => {
            expect(
                container.querySelector("g[role='graphics-symbol'].horizontalLines > path"),
            ).toHaveAttribute('d', 'M0,0h1v2h-1Z');
            // Click on dropdown. MouseDown event has to be used, because material-ui Select component
            // listens to mouseDown event to expand options menu.
            fireEvent.mouseDown(screen.getByLabelText('Distance Measure'));
            fireEvent.click(await screen.findByText(distanceMeasureOptions[2].label));

            await waitFor(() => {
                expect(
                    container.querySelector("g[role='graphics-symbol'].horizontalLines > path"),
                ).toHaveAttribute('d', 'M-151,0h151v2h-151Z');
            });
        });

        it('should show different hierarchical clustering after user changes linkage function', async () => {
            expect(
                container.querySelector(
                    "g[role='graphics-symbol'].horizontalLines > path:nth-child(1)",
                ),
            ).toHaveAttribute('d', 'M0,0h1v2h-1Z');
            // Click on dropdown. MouseDown event has to be used, because material-ui Select component
            // listens to mouseDown event to expand options menu.
            fireEvent.mouseDown(screen.getByLabelText('Clustering Linkage'));
            fireEvent.click(await screen.findByText(linkageFunctionOptions[2].label));

            await waitFor(() => {
                expect(
                    container.querySelector(
                        "g[role='graphics-symbol'].horizontalLines > path:nth-child(1)",
                    ),
                ).toHaveAttribute('d', 'M-151,0h151v2h-151Z');
            });
        });

        it('should display tooltip on gene expressions heatmap hover', async () => {
            await hoverOverVegaSymbol(container, 'genesExpressionsHeatmap');

            await screen.findByText('Time');
            await screen.findByText('Level');
        });

        it('should change line color on gene expressions heatmap hover', async () => {
            // When no genes are highlighted, all lines should have default color.
            expect(
                container.querySelector(
                    `g[role='graphics-symbol'].horizontalLines > path[fill='${highlightedColor}']`,
                ),
            ).not.toBeInTheDocument();
            expect(
                container.querySelector(
                    `g[role='graphics-symbol'].verticalLines > path[fill='${highlightedColor}']`,
                ),
            ).not.toBeInTheDocument();

            fireEvent.mouseOver(
                container.querySelector(
                    "g[role='graphics-symbol'].genesExpressionsHeatmap > path",
                ) as Element,
            );

            // On mouseOver at least one highlighted line should be in the document.
            await waitFor(() => {
                expect(
                    container.querySelector(
                        `g[role='graphics-symbol'].horizontalLines > path[fill='${highlightedColor}']`,
                    ),
                ).toBeInTheDocument();
                expect(
                    container.querySelector(
                        `g[role='graphics-symbol'].verticalLines > path[fill='${highlightedColor}']`,
                    ),
                ).toBeInTheDocument();
            });
        });

        it('should change line color on horizontalLines hover', async () => {
            // When no genes are highlighted, all lines should have default color.
            expect(
                container.querySelector(
                    `g[role='graphics-symbol'].horizontalLines > path[fill='${highlightedColor}']`,
                ),
            ).not.toBeInTheDocument();
            expect(
                container.querySelector(
                    `g[role='graphics-symbol'].verticalLines > path[fill='${highlightedColor}']`,
                ),
            ).not.toBeInTheDocument();

            fireEvent.mouseOver(
                container.querySelector(
                    "g[role='graphics-symbol'].horizontalLines > path",
                ) as Element,
            );

            // On mouseOver at least one highlighted line should be in the document.
            await waitFor(() => {
                expect(
                    container.querySelector(
                        `g[role='graphics-symbol'].horizontalLines > path[fill='${highlightedColor}']`,
                    ),
                ).toBeInTheDocument();
                expect(
                    container.querySelector(
                        `g[role='graphics-symbol'].verticalLines > path[fill='${highlightedColor}']`,
                    ),
                ).toBeInTheDocument();
            });
        });

        it('should change line color on verticalLines hover', async () => {
            // When no genes are highlighted, all lines should have default color.
            expect(
                container.querySelector(
                    `g[role='graphics-symbol'].horizontalLines > path[fill='${highlightedColor}']`,
                ),
            ).not.toBeInTheDocument();
            expect(
                container.querySelector(
                    `g[role='graphics-symbol'].verticalLines > path[fill='${highlightedColor}']`,
                ),
            ).not.toBeInTheDocument();

            fireEvent.mouseOver(
                container.querySelector(
                    "g[role='graphics-symbol'].verticalLines > path",
                ) as Element,
            );

            // On mouseOver at least one highlighted line should be in the document.
            await waitFor(() => {
                expect(
                    container.querySelector(
                        `g[role='graphics-symbol'].horizontalLines > path[fill='${highlightedColor}']`,
                    ),
                ).toBeInTheDocument();
                expect(
                    container.querySelector(
                        `g[role='graphics-symbol'].verticalLines > path[fill='${highlightedColor}']`,
                    ),
                ).toBeInTheDocument();
            });
        });

        it('should change line color after gene is highlighted', async () => {
            // When no genes are highlighted, all lines should have default color.
            expect(
                container.querySelector(
                    `g[role='graphics-symbol'].horizontalLines > path[fill='${highlightedColor}']`,
                ),
            ).not.toBeInTheDocument();

            fireEvent.click(screen.getByRole('button', { name: genes[0].name }));
            fireEvent.click(await screen.findByText('Highlight'));

            await waitFor(() =>
                expect(
                    container.querySelector(
                        `g[role='graphics-symbol'].horizontalLines > path[fill='${highlightedColor}']`,
                    ),
                ).toBeInTheDocument(),
            );
        });

        it('should highlight a gene when user clicks on a line', async () => {
            // When no genes are highlighted, all lines should have default color.
            expect(
                container.querySelector(
                    `g[role='graphics-symbol'].horizontalLines > path[fill='${highlightedColor}']`,
                ),
            ).not.toBeInTheDocument();

            fireEvent.click(
                container.querySelector(
                    "g[role='graphics-symbol'].genesExpressionsHeatmap > path",
                ) as Element,
            );

            await waitFor(() =>
                expect(
                    container.querySelector(
                        `g[role='graphics-symbol'].horizontalLines > path[fill='${highlightedColor}']`,
                    ),
                ).toBeInTheDocument(),
            );
        });

        it('should export visualization Sample Hierarchical Clustering/image.png file', async () => {
            await validateExportFile('Sample Hierarchical Clustering/image.png', (exportFile) => {
                expect(exportFile).toBeDefined();
            });
        });

        it('should export visualization Sample Hierarchical Clustering/image.svg file', async () => {
            await validateExportFile('Sample Hierarchical Clustering/image.svg', (exportFile) => {
                expect(exportFile).toBeDefined();
            });
        });

        it('should export visualization Sample Hierarchical Clustering/caption.txt file', async () => {
            await validateExportFile('Sample Hierarchical Clustering/caption.txt', (exportFile) => {
                expect(exportFile?.content).toContain(
                    'Hierarchical Clustering plotted as a dendrogram is showing clusters of samples',
                );
            });
        });

        it('should save selected time series, genes, highlighted genes and all component bookmarkable state to app-state api', async () => {
            fireEvent.click(screen.getByLabelText('Bookmark'));

            await validateCreateStateRequest((bookmarkState) => {
                expect(bookmarkState.timeSeries.selectedId).toEqual(
                    initialState.timeSeries.selectedId,
                );
                expect(_.get(bookmarkState, BookmarkStatePath.clusteringDistanceMeasure)).toEqual(
                    DistanceMeasure.euclidean,
                );
                expect(_.get(bookmarkState, BookmarkStatePath.clusteringLinkageFunction)).toEqual(
                    ClusteringLinkageFunction.average,
                );
            });
        });
    });

    it('should load selected time series, genes and genesExpressions controls values from bookmark', async () => {
        ({ container } = customRender(<GeneExpressGrid />, {
            initialState,
            route: generateBookmarkQueryParameter(),
        }));

        await validateChart(
            backendBookmark.state.genes.selectedGenesIds.map(
                (selectedGeneId) => genesById[selectedGeneId],
            ),
        );

        expect(screen.getByLabelText('Distance Measure')).toHaveTextContent(
            distanceMeasureOptions[2].label,
        );
        expect(screen.getByLabelText('Clustering Linkage')).toHaveTextContent(
            linkageFunctionOptions[2].label,
        );
    });
});
