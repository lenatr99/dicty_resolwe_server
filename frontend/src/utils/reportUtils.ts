import _, { isArray } from 'lodash';
import * as vega from 'vega';

export const objectsArrayToTsv = (data: Record<string, unknown>[]): string => {
    if (_.isEmpty(data) || _.isEmpty(data[0])) {
        return '';
    }

    const headers = Object.keys(data[0]);
    return [headers, ...data]
        .map((item) => {
            if (isArray(item)) {
                return item.join('\t');
            }

            return headers.map((header) => item[header]).join('\t');
        })
        .join('\n');
};

export const vegaToPng = async (chartView: vega.View): Promise<string> => {
    const pngBase64 = (await chartView.toCanvas()).toDataURL('image/png');

    return pngBase64.substring(pngBase64.indexOf(',') + 1);
};

export const vegaToSvg = async (chartView: vega.View): Promise<Blob> => {
    return new Blob([await chartView.toSVG()], { type: 'image/svg+xml' });
};
