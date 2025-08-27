import { Gene } from '../redux/models/internal';
import { deserializeResponse } from '../utils/apiUtils';
import { post } from './fetch';
import { apiUrl } from './base';

const baseUrl = `${apiUrl}/_modules/gene_list/list_by_ids?hydrate_full_feature=1`;

export const listByIds = async (
    source: string,
    geneIds: string[],
    species: string,
): Promise<Gene[]> => {
    if (geneIds.length === 0) {
        return [];
    }

    const payload = {
        source,
        gene_ids: geneIds,
        species,
    };

    const getGenesResponse = await post(baseUrl, payload);
    return deserializeResponse<Gene[]>(getGenesResponse);
};
