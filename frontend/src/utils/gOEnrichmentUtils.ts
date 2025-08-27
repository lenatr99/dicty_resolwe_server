import { GOEnrichmentJson } from '@genialis/resolwe/dist/api/types/rest';
import _ from 'lodash';
import {
    EnhancedGOEnrichmentJson,
    Gene,
    GOEnrichmentRow,
    GOEnrichmentTerm,
} from 'redux/models/internal';
import { listByIds } from 'api';
import { AspectValue } from 'components/genexpress/common/constants';

export const appendMissingAttributesToJson = (
    gOEnrichmentJson: EnhancedGOEnrichmentJson,
    source: string,
    species: string,
): void => {
    const walk = (
        items: GOEnrichmentRow[] | undefined,
        fn: (item: GOEnrichmentRow, depth: number) => void,
        depth = 0,
    ): void => {
        _.each(items, (item) => {
            fn(item, depth);
            walk(item.children, fn, depth + 1);
        });
    };

    _.each(gOEnrichmentJson.tree, (aspectRows) => {
        let maxScore = -1;
        walk(aspectRows, (item) => {
            maxScore = Math.max(item.score, maxScore);
        });

        walk(aspectRows, (item, depth) => {
            item.depth = depth;
            item.source = source;
            item.species = species;
            item.score_percentage = (item.score * 100) / maxScore;
            item.gene_associations = gOEnrichmentJson.gene_associations[item.term_id];
            item.collapsed = false;
        });
    });
};

export const flattenGoEnrichmentTree = (
    items: GOEnrichmentRow[] | undefined,
    path: string[],
    disableCollapseCheck = false,
): GOEnrichmentRow[] => {
    return (
        items?.flatMap((item) => {
            const itemPath = [...path, item.term_name];
            return [
                { ...item, path: itemPath },
                ...flattenGoEnrichmentTree(
                    !disableCollapseCheck && item.collapsed ? undefined : item.children,
                    [...itemPath],
                    disableCollapseCheck,
                ),
            ];
        }) ?? []
    );
};

export const getRowId = (data: GOEnrichmentRow, path: string[]) => {
    // Path must also be included as a unique identifier because the same
    // term_id can be found in multiple tree branches.
    return data.term_id + String(path);
};

export const ontologyJsonToOntologyRows = (
    gOEnrichmentJson: GOEnrichmentJson,
    aspect: AspectValue,
    disableCollapseCheck = false,
): GOEnrichmentRow[] => {
    if (_.isEmpty(gOEnrichmentJson)) return [];

    const filterByAspect = (json: GOEnrichmentJson): GOEnrichmentRow[] => {
        return json.tree[aspect] as GOEnrichmentRow[];
    };

    return flattenGoEnrichmentTree(filterByAspect(gOEnrichmentJson), [], disableCollapseCheck);
};

export const ontologyJsonToTermsTable = async (
    gOEnrichmentJson: GOEnrichmentJson,
    aspect: AspectValue,
): Promise<GOEnrichmentTerm[]> => {
    const gOEnrichmentRows = ontologyJsonToOntologyRows(gOEnrichmentJson, aspect);
    const uniqueTerms = _.uniqBy(gOEnrichmentRows, 'term_id');
    let table: GOEnrichmentTerm[] = _.map(uniqueTerms, (row) => {
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
        };
    });

    // Fetch genes for all rows (grouped by source and species).
    const getSourceSpeciesKey = (row: GOEnrichmentTerm | Gene): string =>
        `${row.species || ''}&${row.source || ''}`;
    const rowsBySourceSpecies = _.groupBy(table, (row) => getSourceSpeciesKey(row));

    const sourceSpeciesGenesPromises = _.map(rowsBySourceSpecies, async (rows, key) => {
        const geneIds = _.uniq(rows.flatMap((row) => row.selected_gene_associations));
        const { source, species } = rows[0];

        const genes = await listByIds(source ?? '', geneIds, species ?? '');
        return { sourceSpeciesKey: key, values: _.keyBy(genes, 'feature_id') };
    });

    // Wait until all genes are fetched (one request for each source and species).
    const sourceSpeciesDatas = await Promise.all(sourceSpeciesGenesPromises);
    const genesBySourceSpecies = _.zipObject(
        sourceSpeciesDatas.map((sourceSpeciesData) => sourceSpeciesData.sourceSpeciesKey),
        sourceSpeciesDatas.map((sourceSpeciesData) => sourceSpeciesData.values),
    );

    // Append genes names to GOEnrichmentTerms table.
    table = _.map(table, (row) => {
        const genesById = genesBySourceSpecies[getSourceSpeciesKey(row)];
        const selectedGeneAssociationsNames = _.map(row.selected_gene_associations, (geneId) => {
            return genesById[geneId] ? genesById[geneId].name : 'N/A';
        });

        return {
            ...row,
            selected_gene_associations_names: selectedGeneAssociationsNames,
        };
    });

    return _.sortBy(table, (row) => row.p_value);
};
