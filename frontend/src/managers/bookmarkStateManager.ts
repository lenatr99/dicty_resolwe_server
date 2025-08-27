import { v4 as uuidv4 } from 'uuid';
import _ from 'lodash';
import { BookmarkReduxState, RootState } from 'redux/rootReducer';
import {
    comparisonTimeSeriesChanged,
    getBasketId,
    getComparisonTimeSeries,
    getSelectedTimeSeries,
    setTimeSeriesSelected,
} from 'redux/stores/timeSeries';
import { createAppState, getAppState } from 'api/appState';
import { makeBasketReadOnly } from 'api/basketApi';
import { AppDispatch } from 'redux/appStore';
import { genesHighlighted, getHighlightedGenesIds, getSelectedGenes } from 'redux/stores/genes';
import { fetchAndSelectPredefinedGenes } from 'redux/epics/epicsActions';
import { BookmarkComponentsState } from 'redux/models/internal';
import { getPValueThreshold, pValueThresholdChanged } from 'redux/stores/gOEnrichment';
import {
    genesSimilaritiesDistanceMeasureChanged,
    genesSimilaritiesQueryGeneSet,
    getGenesSimilaritiesDistanceMeasure,
    getGenesSimilaritiesQueryGeneId,
} from 'redux/stores/genesSimilarities';
import {
    differentialExpressionSelected,
    getSelectedDifferentialExpression,
} from 'redux/stores/differentialExpressions';
import { BookmarkStatePath } from 'components/genexpress/common/constants';

type BookmarkableState<T> = {
    bookmarkStatePath: BookmarkStatePath;
    get: () => T;
    set: (value: T) => void;
};

const registeredBookmarkableStates: {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    [id: string]: BookmarkableState<any>;
} = {};

export const register = <T>(bookmarkableState: BookmarkableState<T>): (() => void) => {
    const id = uuidv4();

    registeredBookmarkableStates[id] = bookmarkableState;

    return (): void => {
        delete registeredBookmarkableStates[id];
    };
};

const getBookmarkReduxState = (state: RootState): BookmarkReduxState => {
    const selectedGenes = getSelectedGenes(state.genes);
    return {
        timeSeries: {
            selectedId: getSelectedTimeSeries(state.timeSeries)?.id ?? null,
            comparisonIds: getComparisonTimeSeries(state.timeSeries).map(
                (timeSeries) => timeSeries.id,
            ),
        },
        genes: {
            selectedGenesIds: selectedGenes.map((gene) => gene.feature_id),
            highlightedGenesIds: getHighlightedGenesIds(state.genes),
            source: selectedGenes[0]?.source,
            species: selectedGenes[0]?.species,
        },
        gOEnrichment: {
            pValueThreshold: getPValueThreshold(state.gOEnrichment),
        },
        genesSimilarities: {
            distanceMeasure: getGenesSimilaritiesDistanceMeasure(state.genesSimilarities),
            queryGeneId: getGenesSimilaritiesQueryGeneId(state.genesSimilarities),
        },
        differentialExpressions: {
            selectedId:
                getSelectedDifferentialExpression(state.differentialExpressions)?.id ?? null,
        },
    };
};

const getComponentsStates = (): BookmarkComponentsState => {
    const componentsStates: BookmarkComponentsState = {};

    _.each(registeredBookmarkableStates, (bookmarkableState) => {
        _.set(componentsStates, bookmarkableState.bookmarkStatePath, bookmarkableState.get());
    });

    return componentsStates;
};

export const saveBookmarkState = async (state: RootState): Promise<string> => {
    const appStateId = uuidv4();

    const basketId = getBasketId(state.timeSeries);
    if (basketId != null) {
        await makeBasketReadOnly(basketId);
    }

    const bookmarkedState = { ...getBookmarkReduxState(state), ...getComponentsStates() };
    await createAppState({ contributor: 1, uuid: appStateId, state: bookmarkedState });

    return appStateId;
};

export const loadBookmarkedState = async (
    appStateId: string,
    dispatch: AppDispatch,
): Promise<void> => {
    const bookmarkedState = await getAppState(appStateId);

    if (bookmarkedState.timeSeries.selectedId != null) {
        dispatch(setTimeSeriesSelected(bookmarkedState.timeSeries.selectedId));
    }

    if (
        bookmarkedState.timeSeries.comparisonIds != null &&
        bookmarkedState.timeSeries.comparisonIds.length > 0
    ) {
        dispatch(comparisonTimeSeriesChanged(bookmarkedState.timeSeries.comparisonIds));
    }

    if (bookmarkedState.genes.selectedGenesIds.length > 0) {
        dispatch(
            fetchAndSelectPredefinedGenes({
                geneIds: bookmarkedState.genes.selectedGenesIds,
                source: bookmarkedState.genes.source,
                species: bookmarkedState.genes.species,
            }),
        );

        dispatch(genesHighlighted(bookmarkedState.genes.highlightedGenesIds));
    }

    if (bookmarkedState.differentialExpressions.selectedId != null) {
        dispatch(
            differentialExpressionSelected(bookmarkedState.differentialExpressions.selectedId),
        );
    }

    dispatch(pValueThresholdChanged(bookmarkedState.gOEnrichment.pValueThreshold));

    dispatch(
        genesSimilaritiesDistanceMeasureChanged(bookmarkedState.genesSimilarities.distanceMeasure),
    );
    dispatch(genesSimilaritiesQueryGeneSet(bookmarkedState.genesSimilarities.queryGeneId));

    _.map(registeredBookmarkableStates, (bookmarkableState) => {
        const bookmarkedStateValue = _.get(bookmarkedState, bookmarkableState.bookmarkStatePath);
        if (bookmarkedStateValue != null) {
            bookmarkableState.set(_.get(bookmarkedState, bookmarkableState.bookmarkStatePath));
        }
    });
};
