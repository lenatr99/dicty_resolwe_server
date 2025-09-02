import React, {
    ReactElement,
    useMemo,
    useRef,
    useEffect,
    useCallback,
    useState,
    forwardRef,
    useImperativeHandle,
} from 'react';
import { connect, ConnectedProps } from 'react-redux';
import { createSelector } from '@reduxjs/toolkit';
import { shallowEqual } from 'react-redux';
import _ from 'lodash';
// Inline canvas renderer to avoid extra files
import { UmapVisualizationContainer } from './umapVisualization.styles';
import { RootState } from '../../../../redux/rootReducer';
import { getSelectedTimeSeries } from '../../../../redux/stores/timeSeries';
import { getSamplesExpressionsById } from '../../../../redux/stores/samplesExpressions';
import { getHighlightedGenesIds, getSelectedGenes } from '../../../../redux/stores/genes';

export interface UmapDataPoint {
    id: string;
    x: number;
    y: number;
    marker?: string;
    time?: string;
}

// Inline Canvas Scatter Plot (viridis palette, smooth zoom)
interface ViewTransform {
    scale: number;
    offsetX: number;
    offsetY: number;
}
export interface CanvasScatterPlotHandle {
    getCanvas: () => HTMLCanvasElement | null;
    exportAsPNG: () => string | null;
}
interface ScatterPlotPoint extends UmapDataPoint {
    screenX: number;
    screenY: number;
    highlighted: boolean;
}
type ColorValues = Record<string, number> | undefined;

type ColorMode = 'genes' | 'time';
// Helpers for color and drawing
const ZERO_COLOR = '#e8e8e8';
const clamp01 = (t: number) => Math.max(0, Math.min(1, t));
const lerp = (a: number, b: number, t: number) => a + (b - a) * t;
const genesColor = (t: number) => {
    const start = { r: 232, g: 232, b: 232 };
    const end = { r: 74, g: 104, b: 141 };
    const tt = clamp01(t);
    return `rgb(${Math.round(lerp(start.r, end.r, tt))}, ${Math.round(lerp(start.g, end.g, tt))}, ${Math.round(lerp(start.b, end.b, tt))})`;
};
const hexToRgb = (hex: string) => {
    const h = hex.replace('#', '');
    return {
        r: parseInt(h.slice(0, 2), 16),
        g: parseInt(h.slice(2, 4), 16),
        b: parseInt(h.slice(4, 6), 16),
    };
};
const spectralStops = [
    '#ea5646',
    '#f46a9c',
    '#ef9b1f',
    '#ede15a',
    '#86bc44',
    '#26adef',
];
const spectralColor = (t: number) => {
    const n = spectralStops.length - 1;
    const x = clamp01(t) * n;
    const i = Math.floor(x);
    const f = Math.min(1, x - i);
    const c0 = hexToRgb(spectralStops[i]);
    const c1 = hexToRgb(spectralStops[Math.min(n, i + 1)]);
    return `rgb(${Math.round(lerp(c0.r, c1.r, f))}, ${Math.round(lerp(c0.g, c1.g, f))}, ${Math.round(lerp(c0.b, c1.b, f))})`;
};
const computeLegendMax = (entries: Array<[string, number]>, useLog1p: boolean): number => {
    let max = 0;
    for (const [, v] of entries) {
        if (typeof v === 'number' && isFinite(v)) {
            const vt = useLog1p ? Math.log1p(v) : v;
            if (vt > max) max = vt;
        }
    }
    return max > 0 ? max : 1;
};
const buildColorCache = (
    entries: Array<[string, number]>,
    useLog1p: boolean,
    colorMode: ColorMode,
    legendMax: number,
): Record<string, string> => {
    const cache: Record<string, string> = {};
    for (const [cellId, rawVal] of entries) {
        const raw = typeof rawVal === 'number' && isFinite(rawVal) ? rawVal : 0;
        const val = useLog1p ? Math.log1p(raw) : raw;
        const t = val / legendMax;
        cache[cellId] =
            colorMode === 'genes' ? (val <= 0 ? ZERO_COLOR : genesColor(t)) : spectralColor(t);
    }
    return cache;
};
const drawCircle = (ctx: CanvasRenderingContext2D, x: number, y: number, r: number) => {
    ctx.beginPath();
    ctx.arc(x, y, r, 0, 2 * Math.PI);
    ctx.fill();
};
const renderGenesPoints = (
    ctx: CanvasRenderingContext2D,
    points: { screenX: number; screenY: number; id: string }[],
    colorCache: Record<string, string>,
    getValue: (id: string) => number,
    radius: number,
    defaultColor: string,
) => {
    const zero: typeof points = [];
    const nonZero: { p: (typeof points)[number]; v: number }[] = [];
    for (const p of points) {
        const c = colorCache[p.id] ?? defaultColor;
        if (c === ZERO_COLOR) zero.push(p);
        else nonZero.push({ p, v: getValue(p.id) });
    }
    nonZero.sort((a, b) => a.v - b.v);
    ctx.fillStyle = ZERO_COLOR;
    for (const p of zero) drawCircle(ctx, p.screenX, p.screenY, radius);
    for (const { p } of nonZero) {
        ctx.fillStyle = colorCache[p.id] ?? defaultColor;
        drawCircle(ctx, p.screenX, p.screenY, radius);
    }
};
const renderTimePoints = (
    ctx: CanvasRenderingContext2D,
    points: { screenX: number; screenY: number; id: string }[],
    colorCache: Record<string, string>,
    radius: number,
    defaultColor: string,
) => {
    for (const p of points) {
        ctx.fillStyle = colorCache[p.id] ?? defaultColor;
        drawCircle(ctx, p.screenX, p.screenY, radius);
    }
};
const CanvasScatterPlot = forwardRef<
    CanvasScatterPlotHandle,
    {
        data: UmapDataPoint[];
        highlightedPoints?: string[];
        colorValues?: ColorValues;
        onPointClick?: (p: UmapDataPoint) => void;
        tooltipValuesByCellId?: Record<string, Record<string, number>>;
        tooltipGeneOrder?: string[];
        colorMode?: ColorMode;
    }
>(
    (
        {
            data,
            highlightedPoints = [],
            colorValues,
            onPointClick,
            tooltipValuesByCellId,
            tooltipGeneOrder = [],
            colorMode = 'genes',
        },
        ref,
    ) => {
        const canvasRef = useRef<HTMLCanvasElement>(null);
        const containerRef = useRef<HTMLDivElement>(null);
        const spatialIndexRef = useRef<any>(null);
        const colorCacheRef = useRef<Record<string, string>>({});
        const legendMaxRef = useRef<number>(1);
        const [isHovering, setIsHovering] = useState(false);
        const [hoveredPoint, setHoveredPoint] = useState<ScatterPlotPoint | null>(null);
        const [dimensions, setDimensions] = useState({ width: 800, height: 600 });
        const [isDragging, setIsDragging] = useState(false);
        const [dragStart, setDragStart] = useState({ x: 0, y: 0 });
        const [lastTransform, setLastTransform] = useState<ViewTransform>({
            scale: 1,
            offsetX: 0,
            offsetY: 0,
        });
        const transformRef = useRef<ViewTransform>({ scale: 1, offsetX: 0, offsetY: 0 });
        const wheelAccumRef = useRef(0);
        const wheelRafRef = useRef<number | null>(null);
        const devicePixelRatio = window.devicePixelRatio || 1;
        const lastScreenPointsRef = useRef<ScatterPlotPoint[]>([]);
        const [zoomDisplay, setZoomDisplay] = useState(100);
        const [useLog1p, setUseLog1p] = useState(false);
        const [legendMax, setLegendMax] = useState(1);
        const [uniqueTimeValues, setUniqueTimeValues] = useState<number[]>([]);
        const highlightedSet = new Set(highlightedPoints);

        useEffect(() => {
            const container = containerRef.current;
            if (!container) return;
            const resizeObserver = new ResizeObserver((entries) => {
                for (const entry of entries) {
                    const { width, height } = entry.contentRect;
                    setDimensions({ width: Math.max(200, width), height: Math.max(200, height) });
                }
            });
            resizeObserver.observe(container);
            const rect = container.getBoundingClientRect();
            setDimensions({ width: Math.max(200, rect.width), height: Math.max(200, rect.height) });
            return () => resizeObserver.disconnect();
        }, []);

        useEffect(() => {
            const container = containerRef.current;
            if (!container) return;
            const handleNativeWheel = (event: WheelEvent) => {
                event.preventDefault();
            };
            container.addEventListener('wheel', handleNativeWheel, {
                passive: false,
                capture: true,
            });
            return () =>
                container.removeEventListener('wheel', handleNativeWheel, { capture: true } as any);
        }, []);

        useImperativeHandle(ref, () => ({
            getCanvas: () => canvasRef.current,
            exportAsPNG: () =>
                canvasRef.current ? canvasRef.current.toDataURL('image/png') : null,
        }));

        const setupCanvas = useCallback(
            (canvas: HTMLCanvasElement, ctx: CanvasRenderingContext2D) => {
                canvas.width = dimensions.width * devicePixelRatio;
                canvas.height = dimensions.height * devicePixelRatio;
                canvas.style.width = dimensions.width + 'px';
                canvas.style.height = dimensions.height + 'px';
                ctx.scale(devicePixelRatio, devicePixelRatio);
                ctx.imageSmoothingEnabled = true;
                ctx.imageSmoothingQuality = 'high';
            },
            [dimensions.width, dimensions.height, devicePixelRatio],
        );

        const calculateScreenCoordinates = useCallback(
            (points: UmapDataPoint[]): ScatterPlotPoint[] => {
                if (points.length === 0) return [];
                let minX = points[0].x,
                    maxX = points[0].x,
                    minY = points[0].y,
                    maxY = points[0].y;
                for (const p of points) {
                    minX = Math.min(minX, p.x);
                    maxX = Math.max(maxX, p.x);
                    minY = Math.min(minY, p.y);
                    maxY = Math.max(maxY, p.y);
                }
                const padding = 40;
                const plotWidth = dimensions.width - 2 * padding;
                const plotHeight = dimensions.height - 2 * padding;
                const xScale = plotWidth / (maxX - minX);
                const yScale = plotHeight / (maxY - minY);
                const mapped = points.map((point) => {
                    const baseX = padding + (point.x - minX) * xScale;
                    const baseY = dimensions.height - padding - (point.y - minY) * yScale;
                    return {
                        ...point,
                        screenX:
                            (baseX - dimensions.width / 2) * transformRef.current.scale +
                            dimensions.width / 2 +
                            transformRef.current.offsetX,
                        screenY:
                            (baseY - dimensions.height / 2) * transformRef.current.scale +
                            dimensions.height / 2 +
                            transformRef.current.offsetY,
                        highlighted: highlightedSet.has(point.id),
                    };
                });
                lastScreenPointsRef.current = mapped;
                setZoomDisplay(Math.round(transformRef.current.scale * 100));
                return mapped;
            },
            [dimensions.width, dimensions.height, highlightedSet],
        );

        useEffect(() => {
            const entries = colorValues ? Object.entries(colorValues) : [];
            if (colorMode === 'genes') {
                const effectiveMax = computeLegendMax(entries, useLog1p);
                legendMaxRef.current = effectiveMax;
                setLegendMax(effectiveMax);
                colorCacheRef.current = buildColorCache(
                    entries,
                    useLog1p,
                    colorMode,
                    effectiveMax,
                );
                setUniqueTimeValues([]);
            } else {
                // color by time: discrete categories
                const values: number[] = [];
                for (const [, v] of entries) {
                    if (typeof v === 'number' && isFinite(v)) values.push(v);
                }
                const uniqueSorted = Array.from(new Set(values)).sort((a, b) => a - b);
                setUniqueTimeValues(uniqueSorted);
                const n = Math.max(1, uniqueSorted.length - 1);
                const colorByValue = new Map<number, string>();
                uniqueSorted.forEach((val, i) => {
                    const t = n === 0 ? 0 : i / n;
                    colorByValue.set(val, spectralColor(t));
                });
                const cache: Record<string, string> = {};
                for (const [cellId, rawVal] of entries) {
                    const val = typeof rawVal === 'number' && isFinite(rawVal) ? rawVal : 0;
                    cache[cellId] = colorByValue.get(val) ?? ZERO_COLOR;
                }
                colorCacheRef.current = cache;
            }
        }, [colorValues, useLog1p, colorMode]);

        const render = useCallback(
            (points: ScatterPlotPoint[]) => {
                const canvas = canvasRef.current;
                if (!canvas) return;
                const ctx = canvas.getContext('2d');
                if (!ctx) return;
                setupCanvas(canvas, ctx);
                ctx.fillStyle = '#ffffff';
                ctx.fillRect(0, 0, dimensions.width, dimensions.height);
                const normalRadius = 1.2;
                const normalOpacity = 0.9;
                const defaultColor = ZERO_COLOR;
                ctx.globalAlpha = normalOpacity;
                if (colorMode === 'genes') {
                    const getValue = (id: string): number => {
                        const v =
                            colorValues && typeof colorValues[id] === 'number'
                                ? (colorValues[id] as number)
                                : 0;
                        return useLog1p ? Math.log1p(v) : v;
                    };
                    renderGenesPoints(
                        ctx,
                        points,
                        colorCacheRef.current,
                        getValue,
                        normalRadius,
                        defaultColor,
                    );
                } else {
                    renderTimePoints(
                        ctx,
                        points,
                        colorCacheRef.current,
                        normalRadius,
                        defaultColor,
                    );
                }
                if (hoveredPoint) {
                    ctx.globalAlpha = 1.0;
                    ctx.strokeStyle = '#000';
                    ctx.lineWidth = 2;
                    ctx.beginPath();
                    ctx.arc(
                        hoveredPoint.screenX,
                        hoveredPoint.screenY,
                        normalRadius + 2,
                        0,
                        2 * Math.PI,
                    );
                    ctx.stroke();
                }
                ctx.globalAlpha = 1.0;
            },
            [
                dimensions.width,
                dimensions.height,
                hoveredPoint,
                setupCanvas,
                colorValues,
                colorMode,
            ],
        );

        const getMouseCoordinates = useCallback(
            (event: React.MouseEvent<HTMLCanvasElement> | React.WheelEvent<HTMLCanvasElement>) => {
                const canvas = canvasRef.current;
                if (!canvas) return { x: 0, y: 0 };
                const rect = canvas.getBoundingClientRect();
                return { x: event.clientX - rect.left, y: event.clientY - rect.top };
            },
            [],
        );

        const handleMouseDown = useCallback(
            (event: React.MouseEvent<HTMLCanvasElement>) => {
                event.preventDefault();
                event.stopPropagation();
                const coords = getMouseCoordinates(event);
                setIsDragging(true);
                setDragStart(coords);
                setLastTransform(transformRef.current);
                document.body.style.userSelect = 'none';
            },
            [getMouseCoordinates],
        );

        const handleMouseMove = useCallback(
            (event: React.MouseEvent<HTMLCanvasElement>) => {
                const coords = getMouseCoordinates(event);
                if (isDragging) {
                    event.preventDefault();
                    event.stopPropagation();
                    const deltaX = coords.x - dragStart.x;
                    const deltaY = coords.y - dragStart.y;
                    transformRef.current = {
                        ...lastTransform,
                        offsetX: lastTransform.offsetX + deltaX,
                        offsetY: lastTransform.offsetY + deltaY,
                    };
                    const screenPoints = calculateScreenCoordinates(data);
                    spatialIndexRef.current = screenPoints;
                    render(screenPoints);
                } else if (spatialIndexRef.current) {
                    // hover detect: find nearest within radius
                    const R = 6;
                    let closest: ScatterPlotPoint | null = null;
                    let best = R * R;
                    const pts: ScatterPlotPoint[] = lastScreenPointsRef.current;
                    for (let i = 0; i < pts.length; i++) {
                        const p = pts[i];
                        const dx = p.screenX - coords.x;
                        const dy = p.screenY - coords.y;
                        const d2 = dx * dx + dy * dy;
                        if (d2 <= best) {
                            best = d2;
                            closest = p;
                        }
                    }
                    setHoveredPoint(closest);
                    setIsHovering(closest != null);
                }
            },
            [
                isDragging,
                dragStart,
                lastTransform,
                getMouseCoordinates,
                calculateScreenCoordinates,
                data,
                render,
            ],
        );

        const handleMouseUp = useCallback(() => {
            setIsDragging(false);
            document.body.style.userSelect = '';
        }, []);
        const handleDoubleClick = useCallback(() => {
            // Reset zoom and pan to starting size
            transformRef.current = { scale: 1, offsetX: 0, offsetY: 0 };
            const screenPoints = calculateScreenCoordinates(data);
            spatialIndexRef.current = screenPoints;
            render(screenPoints);
        }, [calculateScreenCoordinates, data, render]);
        const handleMouseLeave = useCallback(() => {
            setHoveredPoint(null);
            setIsHovering(false);
            setIsDragging(false);
            document.body.style.userSelect = '';
        }, []);
        const handleClick = useCallback(() => {}, []);

        const handleWheel = useCallback(
            (event: React.WheelEvent<HTMLDivElement> | React.WheelEvent<HTMLCanvasElement>) => {
                const coords = getMouseCoordinates(event as any);
                wheelAccumRef.current += event.deltaY;
                if (wheelRafRef.current != null) return;
                wheelRafRef.current = window.requestAnimationFrame(() => {
                    const baseScale = transformRef.current.scale;
                    const k = 0.0015;
                    const ratio = Math.exp(-wheelAccumRef.current * k);
                    const newScale = baseScale * ratio;
                    wheelAccumRef.current = 0;
                    wheelRafRef.current = null;
                    if (newScale !== transformRef.current.scale) {
                        const scaleRatio = newScale / transformRef.current.scale;
                        const centerX = dimensions.width / 2;
                        const centerY = dimensions.height / 2;
                        const newOffsetX =
                            (coords.x - centerX) * (1 - scaleRatio) +
                            transformRef.current.offsetX * scaleRatio;
                        const newOffsetY =
                            (coords.y - centerY) * (1 - scaleRatio) +
                            transformRef.current.offsetY * scaleRatio;
                        transformRef.current = {
                            scale: newScale,
                            offsetX: newOffsetX,
                            offsetY: newOffsetY,
                        };
                        const screenPoints = calculateScreenCoordinates(data);
                        spatialIndexRef.current = screenPoints;
                        render(screenPoints);
                    }
                });
            },
            [
                getMouseCoordinates,
                calculateScreenCoordinates,
                data,
                render,
                dimensions.width,
                dimensions.height,
            ],
        );

        useEffect(() => {
            const screenPoints = calculateScreenCoordinates(data);
            spatialIndexRef.current = screenPoints;
            render(screenPoints);
        }, [data, colorValues, calculateScreenCoordinates, render]);

        return (
            <div
                ref={containerRef}
                style={{
                    position: 'relative',
                    width: '100%',
                    height: '100%',
                    minHeight: '200px',
                    minWidth: '200px',
                    overflow: 'hidden',
                    touchAction: 'none',
                    overscrollBehavior: 'contain',
                }}
                onWheel={handleWheel}
            >
                {/* Controls */}
                <div
                    style={{
                        position: 'absolute',
                        top: 10,
                        right: 10,
                        display: 'flex',
                        flexDirection: 'column',
                        gap: '4px',
                        zIndex: 1000,
                    }}
                >
                    <button
                        onClick={handleDoubleClick}
                        style={{
                            padding: '4px 8px',
                            fontSize: '12px',
                            background: 'rgba(255, 255, 255, 0.9)',
                            border: '1px solid #ccc',
                            borderRadius: '4px',
                            cursor: 'pointer',
                        }}
                        title="Reset zoom and pan (or double-click)"
                    >
                        Reset View
                    </button>
                    <div
                        style={{
                            padding: '4px 8px',
                            fontSize: '12px',
                            background: 'rgba(255, 255, 255, 0.9)',
                            border: '1px solid #ccc',
                            borderRadius: '4px',
                            textAlign: 'center',
                        }}
                    >
                        Zoom: {zoomDisplay}%
                    </div>
                    {colorMode === 'genes' && (
                        <label
                            style={{
                                padding: '4px 8px',
                                fontSize: '12px',
                                background: 'rgba(255, 255, 255, 0.9)',
                                border: '1px solid #ccc',
                                borderRadius: '4px',
                                display: 'flex',
                                alignItems: 'center',
                                gap: '6px',
                            }}
                        >
                            <input
                                type="checkbox"
                                checked={useLog1p}
                                onChange={(e) => setUseLog1p(e.target.checked)}
                            />{' '}
                            log1p
                        </label>
                    )}
                </div>
                <canvas
                    ref={canvasRef}
                    style={{
                        cursor: isDragging ? 'grabbing' : isHovering ? 'pointer' : 'grab',
                        // border: '1px solid #e0e0e0',
                        borderRadius: '4px',
                        boxSizing: 'border-box',
                        width: `${dimensions.width}px`,
                        height: `${dimensions.height}px`,
                    }}
                    onMouseDown={handleMouseDown}
                    onMouseMove={handleMouseMove}
                    onMouseUp={handleMouseUp}
                    onMouseLeave={handleMouseLeave}
                    onClick={handleClick}
                    onDoubleClick={handleDoubleClick}
                    onWheel={handleWheel}
                />
                {/* Legend */}
                {colorValues && Object.keys(colorValues).length > 0 && (
                    <div
                        style={{
                            position: 'absolute',
                            bottom: 10,
                            right: 10,
                            background: 'rgba(255,255,255,0.9)',
                            padding: '6px 8px',
                            borderRadius: '4px',
                            border: '1px solid #ccc',
                            width: colorMode === 'time' ? '70px' : '160px',
                            maxHeight: '40%',
                            overflowY: 'auto',
                        }}
                    >
                        {colorMode === 'genes' ? (
                            <>
                                <div
                                    style={{
                                        height: '10px',
                                        background: 'linear-gradient(to right, #e8e8e8 0%, #4a688d 100%)',
                                        borderRadius: '2px',
                                    }}
                                />
                                <div
                                    style={{
                                        display: 'flex',
                                        justifyContent: 'space-between',
                                        marginTop: '4px',
                                        fontSize: '12px',
                                    }}
                                >
                                    <span>0</span>
                                    <span>{legendMax.toFixed(2)}</span>
                                </div>
                            </>
                        ) : (
                            <div style={{ display: 'flex', flexDirection: 'column', gap: '0px' }}>
                                {uniqueTimeValues.map((t, i) => (
                                    <div
                                        key={`${t}-${i}`}
                                        style={{
                                            display: 'flex',
                                            alignItems: 'center',
                                            gap: '8px',
                                            fontSize: '12px',
                                        }}
                                    >
                                        <div
                                            style={{
                                                width: '8px',
                                                height: '8px',
                                                borderRadius: '12px',
                                                background:
                                                    uniqueTimeValues.length <= 1
                                                        ? spectralColor(0)
                                                        : spectralColor(i / (uniqueTimeValues.length - 1)),
                                            }}
                                        />
                                        <span>{String(t)} hr</span>
                                    </div>
                                ))}
                            </div>
                        )}
                    </div>
                )}
                {hoveredPoint && tooltipValuesByCellId && (
                    <div
                        style={{
                            position: 'absolute',
                            left: Math.min(
                                (hoveredPoint.screenX ?? 0) + 10,
                                dimensions.width - 240,
                            ),
                            top: Math.max((hoveredPoint.screenY ?? 0) - 30, 10),
                            background: 'rgba(0,0,0,0.85)',
                            color: 'white',
                            padding: '6px 8px',
                            borderRadius: '4px',
                            fontSize: '12px',
                            pointerEvents: 'none',
                            zIndex: 1000,
                            maxWidth: '280px',
                        }}
                    >
                        <div style={{ fontWeight: 600, marginBottom: '4px' }}>
                            Cell: {hoveredPoint.id}
                        </div>
                        {tooltipGeneOrder.map((gene) => (
                            <div
                                key={gene}
                                style={{
                                    display: 'flex',
                                    justifyContent: 'space-between',
                                    gap: '12px',
                                }}
                            >
                                <span>{gene}</span>
                                <span>
                                    {(() => {
                                        const raw =
                                            tooltipValuesByCellId[hoveredPoint.id]?.[gene] ?? 0;
                                        return (useLog1p ? Math.log1p(raw) : raw).toFixed(3);
                                    })()}
                                </span>
                            </div>
                        ))}
                    </div>
                )}
            </div>
        );
    },
);

interface UmapVisualizationProps {
    storageId?: number;
    geneExpressionsStorageId?: number;
}

// Optimized selector with enhanced memoization - following differential expressions pattern
const getUmapCoordinates = createSelector(
    (state: RootState) => getSelectedTimeSeries(state.timeSeries),
    (state: RootState) => getSamplesExpressionsById(state.samplesExpressions),
    (selectedTimeSeries, samplesExpressionsById) => {
        if (!selectedTimeSeries || _.isEmpty(samplesExpressionsById)) {
            return [];
        }

        // Find partitions for UMAP coordinates
        const umapXPartitions = selectedTimeSeries.partitions.filter(
            (partition) => partition.label?.toLowerCase() === 'x_uce_umap_0',
        );
        const umapYPartitions = selectedTimeSeries.partitions.filter(
            (partition) => partition.label?.toLowerCase() === 'x_uce_umap_1',
        );

        if (umapXPartitions.length === 0 || umapYPartitions.length === 0) {
            return [];
        }

        const points: UmapDataPoint[] = [];

        // Get all cell IDs from X coordinates
        umapXPartitions.forEach((xPartition) => {
            const xExpressions = samplesExpressionsById[xPartition.entity];
            if (!xExpressions) return;

            umapYPartitions.forEach((yPartition) => {
                const yExpressions = samplesExpressionsById[yPartition.entity];
                if (!yExpressions) return;

                // Combine X and Y coordinates for each cell
                Object.keys(xExpressions).forEach((cellId) => {
                    const x = xExpressions[cellId];
                    const y = yExpressions[cellId];

                    if (x !== undefined && y !== undefined && !isNaN(x) && !isNaN(y)) {
                        points.push({
                            id: cellId,
                            x,
                            y,
                        });
                    }
                });
            });
        });

        return points;
    },
    { memoizeOptions: { resultEqualityCheck: shallowEqual } },
);

// Compute highlighted cell IDs based on selected genes' single-cell storages
const getHighlightedCellIdsAndValues = createSelector(
    (state: RootState) => getSelectedTimeSeries(state.timeSeries),
    (state: RootState) => getSelectedGenes(state.genes),
    (state: RootState) => getSamplesExpressionsById(state.samplesExpressions),
    (state: RootState) => state.singleCellSeries.bySlug as Record<string, any>,
    (selectedTimeSeries, selectedGenes, samplesExpressionsById, singleCellBySlug) => {
        if (!selectedTimeSeries || selectedGenes.length === 0)
            return { ids: [] as string[], values: {} as Record<string, number> };
        const scSlug = `${selectedTimeSeries.slug}_sc`;
        const scRelation = singleCellBySlug?.[scSlug];
        if (!scRelation) return { ids: [] as string[], values: {} as Record<string, number> };

        const selectedGeneNames = new Set(selectedGenes.map((g) => g.name));
        const selectedGeneIds = new Set(selectedGenes.map((g) => g.feature_id));
        const highlighted = new Set<string>();
        const cellSum: Record<string, number> = {};
        const cellCount: Record<string, number> = {};
        scRelation.partitions
            .filter(
                (p: any) =>
                    p.label && (selectedGeneNames.has(p.label) || selectedGeneIds.has(p.label)),
            )
            .forEach((p: any) => {
                const expr = samplesExpressionsById[p.entity];
                if (!expr) return;
                Object.entries(expr).forEach(([cellId, value]) => {
                    if (typeof value === 'number') {
                        cellSum[cellId] = (cellSum[cellId] ?? 0) + value;
                        cellCount[cellId] = (cellCount[cellId] ?? 0) + 1;
                    }
                });
            });
        const valueMap: Record<string, number> = {};
        Object.keys(cellSum).forEach((cellId) => {
            const mean = cellSum[cellId] / (cellCount[cellId] || 1);
            valueMap[cellId] = mean;
            if (mean > 0) highlighted.add(cellId);
        });
        return { ids: Array.from(highlighted), values: valueMap };
    },
);

// Color by time when no genes are selected
const getCellTimeValues = createSelector(
    (state: RootState) => getSelectedTimeSeries(state.timeSeries),
    (state: RootState) => getSamplesExpressionsById(state.samplesExpressions),
    (selectedTimeSeries, samplesExpressionsById) => {
        const timeValues: Record<string, number> = {};
        if (!selectedTimeSeries) return timeValues;
        const timePartitions = selectedTimeSeries.partitions.filter(
            (p) => (p.label?.toLowerCase?.() ?? '') === 'time',
        );
        timePartitions.forEach((p) => {
            const expr = samplesExpressionsById[p.entity];
            if (!expr) return;
            Object.entries(expr).forEach(([cellId, value]) => {
                if (typeof value === 'number') {
                    // assign; if duplicate partitions exist, last one wins
                    timeValues[cellId] = value;
                }
            });
        });
        return timeValues;
    },
);

// Build per-cell per-gene values for tooltip (for selected genes)
const getCellGeneValues = createSelector(
    (state: RootState) => getSelectedTimeSeries(state.timeSeries),
    (state: RootState) => getSelectedGenes(state.genes),
    (state: RootState) => getSamplesExpressionsById(state.samplesExpressions),
    (state: RootState) => state.singleCellSeries.bySlug as Record<string, any>,
    (selectedTimeSeries, selectedGenes, samplesExpressionsById, singleCellBySlug) => {
        const result: Record<string, Record<string, number>> = {};
        if (!selectedTimeSeries || selectedGenes.length === 0) return result;
        const scSlug = `${selectedTimeSeries.slug}_sc`;
        const scRelation = singleCellBySlug?.[scSlug];
        if (!scRelation) return result;

        // Map label -> geneName for quick lookup
        const nameByLabel = new Map<string, string>();
        selectedGenes.forEach((g) => {
            nameByLabel.set(g.name, g.name);
            nameByLabel.set(g.feature_id, g.name);
        });

        scRelation.partitions
            .filter((p: any) => p.label && nameByLabel.has(p.label))
            .forEach((p: any) => {
                const geneName = nameByLabel.get(p.label)!;
                const expr = samplesExpressionsById[p.entity];
                if (!expr) return;
                Object.entries(expr).forEach(([cellId, value]) => {
                    if (typeof value === 'number') {
                        if (!result[cellId]) result[cellId] = {};
                        result[cellId][geneName] = value;
                    }
                });
            });
        return result;
    },
);

const mapStateToProps = (state: RootState) => ({
    selectedTimeSeries: getSelectedTimeSeries(state.timeSeries),
    umapCoordinates: getUmapCoordinates(state),
    highlightedGenesIds: getHighlightedGenesIds(state.genes),
    highlightedCells: getHighlightedCellIdsAndValues(state),
    cellGeneValues: getCellGeneValues(state),
    selectedGeneNames: getSelectedGenes(state.genes).map((g) => g.name),
    timeValues: getCellTimeValues(state),
});

const connector = connect(mapStateToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;
type CombinedProps = UmapVisualizationProps & PropsFromRedux;

const UmapVisualization = ({
    selectedTimeSeries,
    umapCoordinates,
    highlightedGenesIds,
    highlightedCells,
    cellGeneValues,
    selectedGeneNames,
    timeValues,
}: CombinedProps): ReactElement => {
    // Chart reference for potential future interactions
    const chartRef = useRef<CanvasScatterPlotHandle>(null);

    // Memoize cell count display to prevent unnecessary re-renders during resizing
    const cellCountDisplay = useMemo(() => {
        return umapCoordinates.length > 0 ? `${umapCoordinates.length} cells` : '';
    }, [umapCoordinates.length]);

    if (!selectedTimeSeries) {
        return (
            <UmapVisualizationContainer>
                <div>Please select a time series to view UMAP visualization.</div>
            </UmapVisualizationContainer>
        );
    }

    if (umapCoordinates.length === 0) {
        return (
            <UmapVisualizationContainer>
                <div>No UMAP coordinates available in the selected time series.</div>
            </UmapVisualizationContainer>
        );
    }

    return (
        <UmapVisualizationContainer>
            <div
                style={{
                    flex: 1,
                    display: 'flex',
                    width: '100%',
                    height: '100%',
                    minHeight: 0, // allow child to shrink in flex containers
                    maxHeight: '100vh', // never exceed viewport height
                    overflow: 'hidden', // prevent overflow below the window
                }}
            >
                <CanvasScatterPlot
                    data={umapCoordinates}
                    highlightedPoints={selectedGeneNames.length > 0 ? highlightedCells.ids : []}
                    colorValues={
                        selectedGeneNames.length > 0 ? highlightedCells.values : timeValues
                    }
                    // @ts-ignore
                    colorMode={selectedGeneNames.length > 0 ? 'genes' : 'time'}
                    // pass tooltip gene values and order
                    // @ts-ignore - inline component accepts these props
                    tooltipValuesByCellId={cellGeneValues}
                    // @ts-ignore
                    tooltipGeneOrder={selectedGeneNames}
                    ref={chartRef}
                />
            </div>
        </UmapVisualizationContainer>
    );
};

export default connector(UmapVisualization);
