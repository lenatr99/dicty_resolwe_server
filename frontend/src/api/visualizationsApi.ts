import { BasketExpressionsResponse } from '../redux/models/rest';
import { deserializeResponse } from '../utils/apiUtils';
import { get } from './fetch';
import { apiUrl } from './base';

const baseUrl = `${apiUrl}/_modules/visualizations`;

export const getBasketExpressions = async (
    basketId: string,
): Promise<BasketExpressionsResponse> => {
    const getBasketExpressionsResponse = await get(`${baseUrl}/basket_expressions`, {
        basket: basketId,
    });

    return deserializeResponse<BasketExpressionsResponse>(getBasketExpressionsResponse);
};
