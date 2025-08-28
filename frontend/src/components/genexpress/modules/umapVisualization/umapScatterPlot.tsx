import React, { ReactElement, useRef, useCallback, forwardRef, useState, useEffect } from 'react';
import { Spec } from 'vega';
import Chart, { ChartHandle, DataDefinition, SignalDefinition } from '../../common/chart/chart';
import { getMinMax } from 'utils/math';
import { GEN_CYAN, GEN_GREY, GEN_ORANGE } from 'components/genexpress/common/theming/theming';
import useStateWithEffect from 'components/genexpress/common/useStateWithEffect';
import { UmapPoint } from 'types/application';

type UmapScatterPlotProps = {
    data: UmapPoint[];
    selectedCellIds: string[];
    highlightedCellIds: string[];
    onSelect: ((cellIds: string[]) => void) | undefined;
    colorByExpression?: boolean;
};

type Range = {
    minX: number;
    maxX: number;
    minY: number;
    maxY: number;
};

const getVegaSpecification = (
    data: UmapPoint[],
    getRange: () => Range,
    selectedCellPoints: UmapPoint[],
    highlightedCellPoints: UmapPoint[],
    colorByExpression: boolean = false,
): Spec => {
    const range = getRange();

    return {
        signals: [
            {
                name: 'range',
                value: range,
            },
            {
                name: 'disableSelection',
                value: false,
            },

        ],
        data: [
            {
                name: 'table',
                values: data,
            },
            {
                name: 'highlighted',
                values: highlightedCellPoints,
            },
            {
                name: 'selected',
                values: selectedCellPoints,
            },

        ],
        scales: [
            {
                name: 'xscale',
                type: 'linear',
                range: 'width',
                domain: [range.minX, range.maxX],
                nice: true,
            },
            {
                name: 'yscale',
                type: 'linear',
                range: 'height',
                domain: [range.minY, range.maxY],
                nice: true,
            },
            ...(colorByExpression ? [{
                name: 'colorscale',
                type: 'linear' as const,
                range: ['#f7f7f7', '#d1e5f0', '#92c5de', '#4393c3', '#2166ac', '#053061'],
                domain: { data: 'table', field: 'expressionValue' },
                zero: false,
            }] : []),
        ],
        axes: [
            {
                orient: 'bottom',
                scale: 'xscale',
                title: 'UMAP 1',
                titleFontSize: 14,
                titleFontWeight: 'bold',
                labelFontSize: 12,
            },
            {
                orient: 'left',
                scale: 'yscale',
                title: 'UMAP 2',
                titleFontSize: 14,
                titleFontWeight: 'bold',
                labelFontSize: 12,
            },
        ],
        marks: [
            {
                name: 'cellPointsRemaining',
                type: 'symbol',
                from: { data: 'table' },
                encode: {
                    enter: {
                        size: { value: 3 * 3 },
                        fill: colorByExpression 
                            ? { field: 'expressionValue', scale: 'colorscale' }
                            : { value: GEN_GREY['500'] },
                        tooltip: colorByExpression
                            ? {
                                signal: `{'Cell ID': datum.cellId, 'UMAP 1': datum.umap_1, 'UMAP 2': datum.umap_2, 'Expression': datum.expressionValue}`,
                            }
                            : {
                                signal: `{'Cell ID': datum.cellId, 'UMAP 1': datum.umap_1, 'UMAP 2': datum.umap_2}`,
                            },
                    },
                    update: {
                        x: { field: 'umap_1', scale: 'xscale' },
                        y: { field: 'umap_2', scale: 'yscale' },
                    },
                },
            },
            {
                name: 'cellPointsSelected',
                type: 'symbol',
                from: { data: 'selected' },
                encode: {
                    enter: {
                        size: { value: 9 * 9 },
                        fill: { value: GEN_GREY['700'] },
                        tooltip: {
                            signal: `{'Cell ID': datum.cellId, 'UMAP 1': datum.umap_1, 'UMAP 2': datum.umap_2}`,
                        },
                    },
                    update: {
                        x: { field: 'umap_1', scale: 'xscale' },
                        y: { field: 'umap_2', scale: 'yscale' },
                    },
                },
            },
            {
                name: 'cellPointsHighlighted',
                type: 'symbol',
                from: { data: 'highlighted' },
                encode: {
                    enter: {
                        size: { value: 9 * 9 },
                        fill: { value: GEN_CYAN['500'] },
                        tooltip: {
                            signal: `{'Cell ID': datum.cellId, 'UMAP 1': datum.umap_1, 'UMAP 2': datum.umap_2}`,
                        },
                    },
                    update: {
                        x: { field: 'umap_1', scale: 'xscale' },
                        y: { field: 'umap_2', scale: 'yscale' },
                    },
                },
            },

        ],
    };
};

const UmapScatterPlot = forwardRef<ChartHandle, UmapScatterPlotProps>(
    ({ data, selectedCellIds, highlightedCellIds, onSelect, colorByExpression = false }, ref): ReactElement => {
        const [selectedCellPoints, setSelectedCellPoints] = useState<UmapPoint[]>([]);
        const [highlightedCellPoints, setHighlightedCellPoints] = useState<UmapPoint[]>([]);

        useEffect(() => {
            setSelectedCellPoints(
                data.filter((point) => selectedCellIds.includes(point.cellId)),
            );
        }, [data, selectedCellIds]);

        useEffect(() => {
            setHighlightedCellPoints(
                data.filter((point) => highlightedCellIds.includes(point.cellId)),
            );
        }, [data, highlightedCellIds]);

        const getRange = useCallback((): Range => {
            if (data.length === 0) {
                return { minX: 0, maxX: 1, minY: 0, maxY: 1 };
            }

            const xValues = data.map((point) => point.umap_1);
            const yValues = data.map((point) => point.umap_2);

            const [minX, maxX] = getMinMax(xValues);
            const [minY, maxY] = getMinMax(yValues);

            // Expand the range by 10% like in differential expression
            const xRange = maxX - minX;
            const yRange = maxY - minY;
            const expandedX = xRange * 0.1;
            const expandedY = yRange * 0.1;

            return {
                minX: Math.floor((minX - expandedX) * 100) / 100,
                maxX: Math.ceil((maxX + expandedX) * 100) / 100,
                minY: Math.floor((minY - expandedY) * 100) / 100,
                maxY: Math.ceil((maxY + expandedY) * 100) / 100,
            };
        }, [data]);

        const updatableDataDefinitions: DataDefinition[] = useStateWithEffect(
            () => [
                {
                    name: 'highlighted',
                    data: highlightedCellPoints,
                },
                {
                    name: 'selected',
                    data: selectedCellPoints,
                },
                {
                    name: 'table',
                    data,
                },
            ],
            [data, highlightedCellPoints, selectedCellPoints],
        );

        const updatableSignalDefinitions: SignalDefinition[] = useStateWithEffect(
            () => [
                {
                    name: 'range',
                    value: getRange(),
                },
                {
                    name: 'disableSelection',
                    value: onSelect == null,
                },
            ],
            [getRange, onSelect],
        );

        // Data handlers - simplified for now
        const dataHandlers = useStateWithEffect(() => [], []);

        const renderSpecification = useRef<Spec>(
            getVegaSpecification(data, getRange, selectedCellPoints, highlightedCellPoints, colorByExpression),
        );

        return (
            <Chart
                updatableDataDefinitions={updatableDataDefinitions}
                updatableSignalDefinitions={updatableSignalDefinitions}
                dataHandlers={dataHandlers}
                vegaSpecification={renderSpecification.current}
                ref={ref}
            />
        );
    },
);

UmapScatterPlot.displayName = 'UmapScatterPlot';

export default UmapScatterPlot;
