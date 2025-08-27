import { Action } from '@reduxjs/toolkit';
import { combineEpics, Epic } from 'redux-observable';
import { Observable } from 'rxjs';
import { catchError } from 'rxjs/operators';
import authenticationEpics from './authenticationEpics';
import connectToServerEpics from './connectToServerEpics';
import genesEpics from './genesEpics';
import timeSeriesEpics from './timeSeriesEpics';
import layoutsEpics from './layoutsEpics';
import gOEnrichmentEpics from './gOEnrichmentEpics';
import gafEpics from './gafEpics';
import findSimilarGenesEpics from './findSimilarGenesEpics';
import samplesExpressionsEpics from './samplesExpressionsEpics';
import differentialExpressionsEpics from './differentialExpressionsEpics';
import ontologyOboEpics from './ontologyOboEpics';
import { logError } from 'utils/errorUtils';
import { RootState } from 'redux/rootReducer';

const rootEpic: Epic<Action, Action, RootState> = (
    action$,
    store$,
    dependencies,
): Observable<Action> => {
    return combineEpics(
        authenticationEpics,
        connectToServerEpics,
        layoutsEpics,
        timeSeriesEpics,
        gafEpics,
        genesEpics,
        gOEnrichmentEpics,
        findSimilarGenesEpics,
        samplesExpressionsEpics,
        differentialExpressionsEpics,
        ontologyOboEpics,
    )(action$, store$, dependencies).pipe(
        catchError((error, source) => {
            logError('RootEpic global error', error);
            return source;
        }),
    );
};

export default rootEpic;
