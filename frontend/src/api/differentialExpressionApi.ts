import { DifferentialExpression } from '../redux/models/internal';
import { deserializeResponse } from '../utils/apiUtils';
import { apiUrl } from './base';
import { get } from './fetch';

const baseUrl = `${apiUrl}/_modules/differential_expression/list`;

export const getDifferentialExpressions = async (): Promise<DifferentialExpression[]> => {
    console.log("Getting differential expressions");
    console.log("baseUrl: ", baseUrl);
    const getDifferentialExpressionsDataResponse = await get(baseUrl, {
        tags: `community:${COMMUNITY_SLUG}`,
    });
    console.log("getDifferentialExpressionsDataResponse: ", getDifferentialExpressionsDataResponse);

    return deserializeResponse<DifferentialExpression[]>(getDifferentialExpressionsDataResponse);
};
