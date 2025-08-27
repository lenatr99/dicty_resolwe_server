import _ from 'lodash';
import {
    appendMissingAttributesToJson,
    ontologyJsonToOntologyRows,
    ontologyJsonToTermsTable,
} from './gOEnrichmentUtils';
import { AspectValue } from 'components/genexpress/common/constants';
import { GOEnrichmentTerm } from 'redux/models/internal';
import { generateGeneOntologyStorageJson, generateGenes } from 'tests/mock';
import { resolveStringifiedObjectPromise } from 'tests/test-utils';

const genesAA = generateGenes(['0', '1'], 'A', 'A');
const genesAB = generateGenes(['0', '1'], 'A', 'B');
const genesBB = generateGenes(['0', '1'], 'B', 'B');
const genes = [...genesAA, ...genesAB, ...genesBB];
const genesById = _.keyBy(genes, 'feature_id');

const gOEnrichmentJson = generateGeneOntologyStorageJson(genes.map((gene) => gene.feature_id));
appendMissingAttributesToJson(gOEnrichmentJson, genesAA[0].source, genesAA[0].species);

gOEnrichmentJson.tree.BP[0].source = 'A';
gOEnrichmentJson.tree.BP[0].species = 'A';
gOEnrichmentJson.tree.BP[1].source = 'A';
gOEnrichmentJson.tree.BP[1].species = 'B';
gOEnrichmentJson.tree.BP[2].source = 'B';
gOEnrichmentJson.tree.BP[2].species = 'B';
gOEnrichmentJson.tree.BP[3].source = 'B';
gOEnrichmentJson.tree.BP[3].species = 'B';
gOEnrichmentJson.tree.BP[4].source = 'B';
gOEnrichmentJson.tree.BP[4].species = 'B';

describe('gOEnrichmentUtils', () => {
    beforeAll(() => {
        fetchMock.resetMocks();

        fetchMock.mockResponse(async (req) => {
            const requestJson = await req.json();
            if (req.url.includes('list_by_ids')) {
                if (requestJson.source === 'A' && requestJson.species === 'A') {
                    return resolveStringifiedObjectPromise({ results: genesAA });
                }
                if (requestJson.source === 'A' && requestJson.species === 'B') {
                    return resolveStringifiedObjectPromise({ results: genesAB });
                }
                if (requestJson.source === 'B' && requestJson.species === 'B') {
                    return resolveStringifiedObjectPromise({ results: genesBB });
                }
            }

            return Promise.reject(new Error(`bad url: ${req.url}`));
        });
    });

    it('should make one request for all genes of same source/species', async () => {
        await ontologyJsonToTermsTable(gOEnrichmentJson, AspectValue.bp);

        // One request per unique source/species combination.
        expect(fetchMock.mock.calls.length).toEqual(3);

        const filterCallsBySourceAndSpecies = (
            array: [string | Request | undefined, RequestInit | undefined][],
            source: string,
            species: string,
        ): [string | Request | undefined, RequestInit | undefined][] => {
            return array.filter((call) => {
                const body = JSON.parse(call[1]?.body as string);
                return body.source === source && body.species === species;
            });
        };

        // Genes of same source/species to be gathered in same request.
        expect(filterCallsBySourceAndSpecies(fetchMock.mock.calls, 'A', 'A').length).toEqual(1);
        expect(filterCallsBySourceAndSpecies(fetchMock.mock.calls, 'A', 'B').length).toEqual(1);
        expect(filterCallsBySourceAndSpecies(fetchMock.mock.calls, 'B', 'B').length).toEqual(1);
    });

    it('should append names', async () => {
        const result = await ontologyJsonToTermsTable(gOEnrichmentJson, AspectValue.bp);

        const expected: GOEnrichmentTerm[] = _.sortBy(
            _.map(ontologyJsonToOntologyRows(gOEnrichmentJson, AspectValue.bp), (row) => {
                return {
                    p_value: row.pval,
                    score: row.score,
                    selected_gene_associations_number: row.matched,
                    all_gene_associations_number: row.total,
                    term_id: row.term_id,
                    term_name: row.term_name,
                    source: row.source,
                    species: row.species,
                    selected_gene_associations: row.gene_ids,
                    selected_gene_associations_names: row.gene_ids.map(
                        (geneId) => genesById[geneId].name,
                    ),
                };
            }),
            (row) => row.p_value,
        );

        expect(result).toEqual(expected);
    });
});
