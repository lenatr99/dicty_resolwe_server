import { BasketAddSamplesResponse, BasketUpdateResponse } from '../redux/models/rest';
import { deserializeResponse } from '../utils/apiUtils';
import { apiUrl } from './base';
import { post } from './fetch';

const baseUrl = `${apiUrl}/basket`;

export const addToBasket = async (samplesIds: number[]): Promise<BasketAddSamplesResponse> => {
    if (samplesIds == null || samplesIds.length === 0) {
        throw new Error('No samples to add to basket.');
    }

    const addSamplesResponse = await post(`${baseUrl}/_/add_samples`, {
        samples: samplesIds,
        only_existing_organism: true,
        annotation_version: 'v2',
    });

    return deserializeResponse<BasketAddSamplesResponse>(addSamplesResponse);
};

export const makeBasketReadOnly = async (basketId: string): Promise<BasketUpdateResponse> => {
    const makeBasketReadOnlyResponse = await post(`${baseUrl}/${basketId}/make_read_only`);

    return deserializeResponse<BasketUpdateResponse>(makeBasketReadOnlyResponse);
};
