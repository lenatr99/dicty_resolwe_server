import * as vega from 'vega';

/**
 * Default Vega config defines visual values to set a visualization's theme.
 * Any configuration provided within the specification itself will take
 * precedence over external configurations passed to the parser.
 */
export const vegaTheme: vega.Config = {
    legend: {
        padding: 10,
        layout: {
            bottom: {
                anchor: 'middle',
                direction: 'vertical',
                center: true,
                margin: 2,
            },
        },
    },
    range: {
        category: [
            '#7291CB',
            '#E78414',
            '#A1B900',
            '#D25886',
            '#3CA5B4',
            '#E3BA22',
            '#183457',
            '#C8D89F',
            '#DDBDCF',
            '#BF2B21',
            '#E3593B',
            '#008D79',
            '#B18C00',
            '#852654',
            '#FF9E80',
            '#8D3B00',
            '#8F6C8B',
            '#057B81',
            '#F3DB4A',
        ],
    },
};
