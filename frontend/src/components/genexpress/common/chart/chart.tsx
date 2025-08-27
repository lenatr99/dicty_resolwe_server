import React, {
    ReactElement,
    useRef,
    useEffect,
    useCallback,
    useContext,
    forwardRef,
    useImperativeHandle,
    ForwardRefRenderFunction,
} from 'react';
import * as vega from 'vega';
import * as vegaTooltip from 'vega-tooltip';
import { SignalValue, Spec } from 'vega';
import { useDispatch } from 'react-redux';
import { useTheme } from '@mui/material';
import { vegaTheme } from '../theming/vegaTheme';
import useUpdateEffect from '../useUpdateEffect';
import useSize from '../useSize';
import { RendererContext } from 'components/common/rendererContext';
import { handleError } from 'utils/errorUtils';

export type DataHandler = {
    name: string;
    handler: (name: string, value: unknown) => void;
};

export type SignalHandler = {
    name: string;
    handler: (name: string, value: SignalValue) => void;
};

export type DataDefinition = {
    name: string;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    data: Array<any>;
};

export type SignalDefinition = {
    name: string;
    value: unknown;
};

type ChartProps = {
    vegaSpecification: Spec;
    updatableDataDefinitions?: DataDefinition[];
    updatableSignalDefinitions?: SignalDefinition[];
    dataHandlers?: DataHandler[];
    signalHandlers?: SignalHandler[];
    onChartResized?: (width: number, height: number) => void;
};

export type ChartHandle = {
    getChartView: () => vega.View | null;
    getSvgImage: () => Promise<Blob | null>;
    getPngImage: () => Promise<string | null>;
    getChartDivElement: () => HTMLDivElement | null;
};

const chartPadding = { top: 20, left: 15, bottom: 15, right: 15 };

/**
 * Chart component is the basis for all charts that want to use Vega charting library.
 *
 * For performance reasons, chart will get reinitialized only if vegaSpecification changes.
 * All other updates are handled with API calls (set data -> run). That's why we must pass
 * updatableDataDefinitions for each data that can get updated.
 *
 * Add all data handlers (that execute after data object changes) via 'dataHandlers'
 * variable. All handlers will be automatically removed after view is removed from DOM.
 *
 * Setting `renderer` to 'svg' will help with debugging, because it's more descriptive.
 */
const Chart: ForwardRefRenderFunction<ChartHandle, ChartProps> = (
    {
        vegaSpecification,
        updatableDataDefinitions,
        updatableSignalDefinitions,
        dataHandlers,
        signalHandlers,
        onChartResized,
    },
    ref,
): ReactElement => {
    const dispatch = useDispatch();
    const renderer = useContext(RendererContext);

    const theme = useTheme();

    const addedDataHandlers = useRef<DataHandler[]>([]);
    const addedSignalHandlers = useRef<DataHandler[]>([]);

    const chartElement = useRef<HTMLDivElement>(null);
    const chartView = useRef<vega.View | null>(null);

    const { width, height } = useSize(chartElement);

    useImperativeHandle(ref, () => ({
        getChartView: (): vega.View | null => {
            return chartView.current;
        },
        getSvgImage: async (): Promise<Blob | null> => {
            if (chartView.current == null) {
                return Promise.resolve(null);
            }
            const blob = await chartView.current.toSVG();
            return new Blob([blob], { type: 'image/svg+xml' });
        },
        getPngImage: async (): Promise<string | null> => {
            if (chartView.current == null) {
                return Promise.resolve(null);
            }
            const pngBase64 = (await chartView.current.toCanvas()).toDataURL('image/png');

            return pngBase64.substring(pngBase64.indexOf(',') + 1);
        },
        getChartDivElement: (): HTMLDivElement | null => chartElement.current,
    }));

    const getChartElementWidth = useCallback((): number => {
        if (!chartElement.current) {
            return 0;
        }

        return chartElement.current.clientWidth;
    }, []);

    const getChartElementHeight = useCallback((): number => {
        if (!chartElement.current) {
            return 0;
        }

        return chartElement.current.clientHeight;
    }, []);

    const runChart = useCallback((): void => {
        try {
            void chartView.current?.runAsync();
        } catch (error) {
            dispatch(handleError('Error generating chart.', error));
        }
    }, [chartView, dispatch]);

    const resizeChart = useCallback((): void => {
        const newWidth = getChartElementWidth();
        const newHeight = getChartElementHeight();
        chartView.current?.width(newWidth);
        chartView.current?.height(newHeight);

        runChart();

        onChartResized?.(newWidth, newHeight);
    }, [getChartElementHeight, getChartElementWidth, onChartResized, runChart]);

    useUpdateEffect(() => {
        resizeChart();
    }, [width, height, resizeChart]);

    useEffect((): undefined | (() => void) => {
        const defaultSpecification: vega.Spec = {
            width: getChartElementWidth(),
            height: getChartElementHeight(),
            background: theme.palette.background.default,
            padding: chartPadding,
            autosize: {
                type: 'fit',
                contains: 'padding',
            },
            signals: [],
            data: [],
            legends: [],
            marks: [],
            scales: [],
            axes: [],
        };
        const renderSpecification = { ...defaultSpecification, ...vegaSpecification };

        const generateChart = (): void => {
            chartView.current = new vega.View(vega.parse(renderSpecification, vegaTheme), {
                renderer,
                container: chartElement.current != null ? chartElement.current : undefined,
            });

            const tooltipHandler = new vegaTooltip.Handler({});
            chartView.current.tooltip(tooltipHandler.call);

            runChart();
        };

        generateChart();

        // Vega Cleanup.
        return (): void => {
            // Unregister any timers and remove any event listeners the visualization has
            // registered on external DOM elements.
            chartView.current?.finalize();
        };
    }, [
        chartView,
        getChartElementHeight,
        getChartElementWidth,
        renderer,
        runChart,
        theme.palette.background.default,
        vegaSpecification,
    ]);

    useEffect(() => {
        if (dataHandlers != null && chartView.current != null) {
            // If any handlers are already attached, first remove them all.
            for (let i = 0; i < addedDataHandlers.current.length; i += 1) {
                const addedDataHandler = addedDataHandlers.current[i];
                chartView.current.removeDataListener(
                    addedDataHandler.name,
                    addedDataHandler.handler,
                );
            }

            // Add new handlers.
            for (let i = 0; i < dataHandlers.length; i += 1) {
                const newDataHandler = dataHandlers[i];
                chartView.current.addDataListener(newDataHandler.name, newDataHandler.handler);
            }

            addedDataHandlers.current = dataHandlers;
        }
    }, [chartView, dataHandlers]);

    useEffect(() => {
        if (signalHandlers != null && chartView.current != null) {
            // If any handlers are already attached, first remove them all.
            for (let i = 0; i < addedSignalHandlers.current.length; i += 1) {
                const addedSignalHandler = addedSignalHandlers.current[i];
                chartView.current.removeSignalListener(
                    addedSignalHandler.name,
                    addedSignalHandler.handler,
                );
            }

            // Add new handlers.
            for (let i = 0; i < signalHandlers.length; i += 1) {
                const newSignalHandler = signalHandlers[i];
                chartView.current.addSignalListener(
                    newSignalHandler.name,
                    newSignalHandler.handler,
                );
            }

            addedSignalHandlers.current = signalHandlers;
        }
    }, [chartView, signalHandlers]);

    useEffect(() => {
        if (updatableDataDefinitions != null) {
            updatableDataDefinitions.forEach((dataDefinition) => {
                chartView.current?.data(dataDefinition.name, dataDefinition.data);
            });

            runChart();
        }
    }, [updatableDataDefinitions, runChart, chartView]);

    useEffect(() => {
        if (updatableSignalDefinitions != null) {
            updatableSignalDefinitions.forEach((signalDefinition) => {
                chartView.current?.signal(signalDefinition.name, signalDefinition.value);
            });

            runChart();
        }
    }, [updatableSignalDefinitions, runChart, chartView]);

    return (
        <div style={{ width: '100%', height: '100%' }}>
            <div id="view" style={{ width: '100%', height: '100%' }} ref={chartElement} />
        </div>
    );
};

export default forwardRef(Chart);
