import { Relation } from '@genialis/resolwe/dist/api/types/rest';
import { deserializeResponse } from '../utils/apiUtils';
import { apiUrl } from './base';
import { get } from './fetch';

const baseUrl = `${apiUrl}/relation`;

export const getTimeSeriesRelations = async (): Promise<Relation[]> => {
    const getRelationsResponse = await get(baseUrl, {
        category: 'Time series',
        tags: `community:${COMMUNITY_SLUG}`,
    });

    return deserializeResponse<Relation[]>(getRelationsResponse);
};

export const getRelationBySlug = async (slug: string): Promise<Relation | null> => {
    const response = await get(baseUrl, { slug });
    const relations = await deserializeResponse<Relation[]>(response);
    return relations.length > 0 ? relations[0] : null;
};
