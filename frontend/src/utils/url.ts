import { DictyUrlQueryParameter } from 'components/genexpress/common/constants';

export const getUrlQueryParameter = (
    url: string,
    parameterName: DictyUrlQueryParameter,
): string | null => {
    return new URLSearchParams(url).get(parameterName);
};
