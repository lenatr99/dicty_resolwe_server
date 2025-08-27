import { Storage } from '@genialis/resolwe/dist/api/types/rest';
import { deserializeResponse } from '../utils/apiUtils';
import { apiUrl } from './base';
import { get } from './fetch';

const baseUrl = `${apiUrl}/storage`;

export const getStorage = async (storageId: number): Promise<Storage> => {
    const url = `${baseUrl}/${storageId}`;

    const getStorageJsonResponse = await get(url);

    return deserializeResponse<Storage>(getStorageJsonResponse);
};
