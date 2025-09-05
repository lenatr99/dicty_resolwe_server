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
const MAX_GENES_TO_SHOW = 500;
// Helpers for color and drawing
const ZERO_COLOR = '#e8e8e8';
const clamp01 = (t: number) => Math.max(0, Math.min(1, t));
const lerp = (a: number, b: number, t: number) => a + (b - a) * t;
const genesColor = (t: number) => {
    const tt = clamp01(t);
    // Piecewise interpolation to match legend gradient:
    // #e8e8e8 (0%) -> #4a688d (80%) -> #375068 (100%)
    const start = hexToRgb('#e8e8e8');
    const mid = hexToRgb('#4a688d');
    const end = hexToRgb('#375068');
    if (tt <= 0.8) {
        const k = tt / 0.8;
        return `rgb(${Math.round(lerp(start.r, mid.r, k))}, ${Math.round(lerp(start.g, mid.g, k))}, ${Math.round(lerp(start.b, mid.b, k))})`;
    }
    const k = (tt - 0.8) / 0.2;
    return `rgb(${Math.round(lerp(mid.r, end.r, k))}, ${Math.round(lerp(mid.g, end.g, k))}, ${Math.round(lerp(mid.b, end.b, k))})`;
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
const computeLegendMax = (entries: Array<[string, number]>, transformMode: 'linear' | 'log2' | 'log1p'): number => {
    let max = 0;
    for (const [, v] of entries) {
        if (typeof v === 'number' && isFinite(v)) {
            let vt = v;
            if (transformMode === 'log1p') {
                vt = Math.log1p(v);
            } else if (transformMode === 'log2') {
                vt = v > 0 ? Math.log2(v + 1) : 0;
            }
            if (vt > max) max = vt;
        }
    }
    return max > 0 ? max : 1;
};
const buildColorCache = (
    entries: Array<[string, number]>,
    transformMode: 'linear' | 'log2' | 'log1p',
    colorMode: ColorMode,
    legendMax: number,
): Record<string, string> => {
    const cache: Record<string, string> = {};
    for (const [cellId, rawVal] of entries) {
        const raw = typeof rawVal === 'number' && isFinite(rawVal) ? rawVal : 0;
        let val = raw;
        if (transformMode === 'log1p') {
            val = Math.log1p(raw);
        } else if (transformMode === 'log2') {
            val = raw > 0 ? Math.log2(raw + 1) : 0;
        }
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
        transformMode?: 'linear' | 'log2' | 'log1p';
        aggregationMode?: 'average' | 'sum' | 'min' | 'max';
        onTransformModeChange?: (mode: 'linear' | 'log2' | 'log1p') => void;
        onAggregationModeChange?: (mode: 'average' | 'sum' | 'min' | 'max') => void;
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
            transformMode = 'log2',
            aggregationMode = 'average',
            onTransformModeChange,
            onAggregationModeChange,
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
                const effectiveMax = computeLegendMax(entries, transformMode);
                legendMaxRef.current = effectiveMax;
                setLegendMax(effectiveMax);
                colorCacheRef.current = buildColorCache(
                    entries,
                    transformMode,
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
        }, [colorValues, transformMode, colorMode]);

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
                        if (transformMode === 'log1p') {
                            return Math.log1p(v);
                        } else if (transformMode === 'log2') {
                            return v > 0 ? Math.log2(v + 1) : 0;
                        }
                        return v;
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
                        <>
                            <div
                                style={{
                                    display: 'flex',
                                    flexDirection: 'column',
                                    gap: '2px',
                                }}
                            >
                                <label
                                    style={{
                                        fontSize: '10px',
                                        fontWeight: 'bold',
                                        color: '#666',
                                        marginBottom: '2px',
                                    }}
                                >
                                    Transform:
                                </label>
                                <select
                                    value={transformMode}
                                    onChange={(e) => onTransformModeChange?.(e.target.value as 'linear' | 'log2' | 'log1p')}
                                    style={{
                                        padding: '4px 6px',
                                        fontSize: '12px',
                                        background: 'rgba(255, 255, 255, 0.9)',
                                        border: '1px solid #ccc',
                                        borderRadius: '4px',
                                        minWidth: '80px',
                                    }}
                                >
                                    <option value="linear">Linear</option>
                                    <option value="log2">Log2</option>
                                    <option value="log1p">Log1p</option>
                                </select>
                            </div>
                            {tooltipGeneOrder.length > 1 && (
                                <div
                                    style={{
                                        display: 'flex',
                                        flexDirection: 'column',
                                        gap: '2px',
                                    }}
                                >
                                    <label
                                        style={{
                                            fontSize: '10px',
                                            fontWeight: 'bold',
                                            color: '#666',
                                            marginBottom: '2px',
                                        }}
                                    >
                                        Aggregation:
                                    </label>
                                    <select
                                        value={aggregationMode}
                                        onChange={(e) => onAggregationModeChange?.(e.target.value as 'average' | 'sum' | 'min' | 'max')}
                                        style={{
                                            padding: '4px 6px',
                                            fontSize: '12px',
                                            background: 'rgba(255, 255, 255, 0.9)',
                                            border: '1px solid #ccc',
                                            borderRadius: '4px',
                                            minWidth: '80px',
                                        }}
                                    >
                                        <option value="average">Average</option>
                                        <option value="sum">Sum</option>
                                        <option value="min">Min</option>
                                        <option value="max">Max</option>
                                    </select>
                                </div>
                            )}
                        </>
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
                                        background: 'linear-gradient(to right, #e8e8e8 0%, #4a688d 80%, #375068 100%)',
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
                {/* UMAP Axis Indicator */}
                <div
                    style={{
                        position: 'absolute',
                        bottom: 10,
                        left: 10,
                        pointerEvents: 'none',
                        color: '#666',
                    }}
                >
                    <svg width="52" height="52" viewBox="0 0 52 52">
                        {/* Horizontal axis (UMAP1) */}
                        <line x1="12" y1="40" x2="50" y2="40" stroke="#666" strokeWidth="1.5" />
                        <polygon points="46,36 52,40 46,44" fill="#666" />
                        <text x="32" y="50" fontSize="9" textAnchor="middle" fill="#666" fontWeight="600">UMAP1</text>

                        {/* Vertical axis (UMAP2) */}
                        <line x1="12" y1="40" x2="12" y2="2" stroke="#666" strokeWidth="1.5" />
                        <polygon points="8,6 12,0 16,6" fill="#666" />
                        <text x="-10" y="4" fontSize="9" fill="#666" fontWeight="600" transform="translate(4,24) rotate(-90)">UMAP2</text>
                    </svg>
                </div>
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
                        {(() => {
                            const cellValues = tooltipValuesByCellId[hoveredPoint.id] || {};
                            const scored = tooltipGeneOrder.map((gene) => {
                                const raw = cellValues[gene] ?? 0;
                                let transformed = raw;
                                if (transformMode === 'log1p') {
                                    transformed = Math.log1p(raw);
                                } else if (transformMode === 'log2') {
                                    transformed = raw > 0 ? Math.log2(raw + 1) : 0;
                                }
                                return { gene, value: transformed };
                            });
                            scored.sort((a, b) => b.value - a.value);
                            const shown = scored.slice(0, 10);
                            const remaining = Math.max(0, scored.length - shown.length);
                            return (
                                <>
                                    {shown.map(({ gene, value }) => (
                                        <div
                                            key={gene}
                                            style={{
                                                display: 'flex',
                                                justifyContent: 'space-between',
                                                gap: '12px',
                                            }}
                                        >
                                            <span>{gene}</span>
                                            <span>{value.toFixed(3)}</span>
                                        </div>
                                    ))}
                                    {remaining > 0 && (
                                        <div style={{ marginTop: '4px', opacity: 0.8 }}>
                                            … and {remaining} more
                                        </div>
                                    )}
                                </>
                            );
                        })()}
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
    (state: RootState, aggregationMode: 'average' | 'sum') => aggregationMode,
    (selectedTimeSeries, selectedGenes, samplesExpressionsById, singleCellBySlug, aggregationMode) => {
        if (!selectedTimeSeries || selectedGenes.length === 0)
            return { ids: [] as string[], values: {} as Record<string, number> };
        const scSlugDash = `${selectedTimeSeries.slug}-sc`;
        const scRelation = singleCellBySlug?.[scSlugDash];
        if (!scRelation) return { ids: [] as string[], values: {} as Record<string, number> };

        // Cap to first MAX_GENES_TO_SHOW genes to avoid performance issues
        const limitedGenes = selectedGenes.slice(0, MAX_GENES_TO_SHOW);
        const selectedGeneNames = new Set(limitedGenes.map((g) => g.name));
        const selectedGeneIds = new Set(limitedGenes.map((g) => g.feature_id));
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
            const aggregatedValue = aggregationMode === 'sum' 
                ? cellSum[cellId] 
                : cellSum[cellId] / (cellCount[cellId] || 1);
            valueMap[cellId] = aggregatedValue;
            if (aggregatedValue > 0) highlighted.add(cellId);
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
        const scSlugDash = `${selectedTimeSeries.slug}-sc`;
        const scRelation = singleCellBySlug?.[scSlugDash];
        if (!scRelation) return result;

        // Map label -> geneName for quick lookup
        const nameByLabel = new Map<string, string>();
        const limitedGenes = selectedGenes.slice(0, MAX_GENES_TO_SHOW);
        limitedGenes.forEach((g) => {
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
    selectedGenes: getSelectedGenes(state.genes),
    samplesExpressionsById: getSamplesExpressionsById(state.samplesExpressions),
    singleCellBySlug: state.singleCellSeries.bySlug as Record<string, any>,
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
    selectedGenes,
    samplesExpressionsById,
    singleCellBySlug,
    cellGeneValues,
    selectedGeneNames,
    timeValues,
}: CombinedProps): ReactElement => {
    // Determine if capping is in effect for the current selection
    const isCapped = selectedGenes.length > MAX_GENES_TO_SHOW;
    const displayedGeneNames = useMemo(
        () => (isCapped ? selectedGenes.slice(0, MAX_GENES_TO_SHOW).map((g) => g.name) : selectedGenes.map((g) => g.name)),
        [selectedGenes, isCapped],
    );
    // Chart reference for potential future interactions
    const chartRef = useRef<CanvasScatterPlotHandle>(null);
    
    // Aggregation mode state
    const [aggregationMode, setAggregationMode] = useState<'average' | 'sum' | 'min' | 'max'>('average');
    const [transformMode, setTransformMode] = useState<'linear' | 'log2' | 'log1p'>('log2');

    // Compute highlighted cells based on aggregation mode
    const highlightedCells = useMemo(() => {
        if (!selectedTimeSeries || selectedGenes.length === 0)
            return { ids: [] as string[], values: {} as Record<string, number> };
        const scSlugDash = `${selectedTimeSeries.slug}-sc`;
        const scRelation = singleCellBySlug?.[scSlugDash];
        if (!scRelation) return { ids: [] as string[], values: {} as Record<string, number> };

        // Build fast lookup of allowed labels from capped genes
        const allowedLabels: Record<string, 1> = Object.create(null);
        const limitedGenes = selectedGenes.slice(0, MAX_GENES_TO_SHOW);
        for (let i = 0; i < limitedGenes.length; i++) {
            const g = limitedGenes[i];
            if (g.name) allowedLabels[g.name] = 1;
            if (g.feature_id) allowedLabels[g.feature_id] = 1;
        }

        const valueMap: Record<string, number> = Object.create(null);
        if (aggregationMode === 'min' || aggregationMode === 'max') {
            // Track per-cell min/max without allocating arrays per cell
            const hasVal: Record<string, 1> = Object.create(null);
            const isMin = aggregationMode === 'min';
            const parts = scRelation.partitions as any[];
            for (let pi = 0; pi < parts.length; pi++) {
                const p = parts[pi];
                const label = p.label;
                if (!label || allowedLabels[label] !== 1) continue;
                const expr = samplesExpressionsById[p.entity];
                if (!expr) continue;
                for (const cellId in expr) {
                    const v = (expr as any)[cellId];
                    if (typeof v !== 'number') continue;
                    if (!hasVal[cellId]) {
                        valueMap[cellId] = v;
                        hasVal[cellId] = 1;
                    } else if (isMin) {
                        if (v < valueMap[cellId]) valueMap[cellId] = v;
                    } else {
                        if (v > valueMap[cellId]) valueMap[cellId] = v;
                    }
                }
            }
            const ids: string[] = [];
            for (const cellId in valueMap) if (valueMap[cellId] > 0) ids.push(cellId);
            return { ids, values: valueMap };
        }

        // Sum/Average optimized path
        const sum: Record<string, number> = Object.create(null);
        const cnt: Record<string, number> = Object.create(null);
        const parts = scRelation.partitions as any[];
        for (let pi = 0; pi < parts.length; pi++) {
            const p = parts[pi];
            const label = p.label;
            if (!label || allowedLabels[label] !== 1) continue;
            const expr = samplesExpressionsById[p.entity];
            if (!expr) continue;
            for (const cellId in expr) {
                const v = (expr as any)[cellId];
                if (typeof v !== 'number') continue;
                sum[cellId] = (sum[cellId] || 0) + v;
                cnt[cellId] = (cnt[cellId] || 0) + 1;
            }
        }
        const ids: string[] = [];
        for (const cellId in sum) {
            const aggregated = aggregationMode === 'sum' ? sum[cellId] : sum[cellId] / (cnt[cellId] || 1);
            valueMap[cellId] = aggregated;
            if (aggregated > 0) ids.push(cellId);
        }
        return { ids, values: valueMap };
    }, [selectedTimeSeries, selectedGenes, samplesExpressionsById, singleCellBySlug, aggregationMode]);

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
                <div>No single cell data available in the selected time series.</div>
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
                    colorMode={selectedGeneNames.length > 0 ? 'genes' : 'time'}
                    tooltipValuesByCellId={cellGeneValues}
                    tooltipGeneOrder={displayedGeneNames}
                    transformMode={transformMode}
                    aggregationMode={aggregationMode}
                    onTransformModeChange={setTransformMode}
                    onAggregationModeChange={setAggregationMode}
                    ref={chartRef}
                />
                {isCapped && (
                    <div
                        style={{
                            position: 'absolute',
                            top: 50,
                            left: 10,
                            background: 'rgba(255, 243, 205, 0.95)',
                            color: '#664d03',
                            border: '1px solid #ffecb5',
                            borderRadius: '4px',
                            padding: '6px 8px',
                            fontSize: '12px',
                            boxShadow: '0 1px 2px rgba(0,0,0,0.05)',
                            zIndex: 1100,
                        }}
                    >
                        Showing first {MAX_GENES_TO_SHOW} genes of {selectedGenes.length} selected.
                    </div>
                )}
            </div>
        </UmapVisualizationContainer>
    );
};

export default connector(UmapVisualization);
