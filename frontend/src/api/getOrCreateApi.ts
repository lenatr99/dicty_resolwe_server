import { Data } from '@genialis/resolwe/dist/api/types/rest';
import { deserializeResponse } from '../utils/apiUtils';
import { post } from './fetch';
import { apiUrl } from './base';
import { ProcessInfo } from 'redux/models/internal';
import { ClusteringData } from 'redux/models/rest';

const baseUrl = `${apiUrl}/data/get_or_create`;

export const getOrCreateData = async <DataType extends Data>(
    input: Record<string, unknown>,
    slug: ProcessInfo['slug'],
): Promise<ClusteringData> => {
    const payload = {
        process: {
            slug,
        },
        input,
    };

    return deserializeResponse<DataType>(await post(baseUrl, payload));
};
