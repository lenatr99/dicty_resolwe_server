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
