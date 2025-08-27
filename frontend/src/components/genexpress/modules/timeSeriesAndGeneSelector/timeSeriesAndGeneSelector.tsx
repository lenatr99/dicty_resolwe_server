import React, { ReactElement } from 'react';
import { connect, ConnectedProps } from 'react-redux';
import { Relation } from '@genialis/resolwe/dist/api/types/rest';
import GeneSelector from './geneSelector/geneSelector/geneSelector';
import TimeSeriesSelector from './timeSeriesSelector/timeSeriesSelector';
import {
    TimeSeriesAndGeneSelectorContainer,
    TimeSeriesSelectorWrapper,
} from './timeSeriesAndGeneSelector.styles';
import { RootState } from 'redux/rootReducer';
import {
    getTimeSeries,
    getTimeSeriesIsFetching,
    getSelectedTimeSeries,
    timeSeriesSelected,
} from 'redux/stores/timeSeries';
import useReport from 'components/genexpress/common/reportBuilder/useReport';
import { objectsArrayToTsv } from 'utils/reportUtils';

export const moduleKey = 'timeSeriesAndGeneSelector';

const mapStateToProps = (state: RootState) => {
    return {
        timeSeries: getTimeSeries(state.timeSeries),
        selectedTimeSeries: getSelectedTimeSeries(state.timeSeries),
        isFetching: getTimeSeriesIsFetching(state.timeSeries),
    };
};

const connector = connect(mapStateToProps, {
    connectedTimeSeriesSelected: timeSeriesSelected,
});

type PropsFromRedux = ConnectedProps<typeof connector>;

const TimeSeriesAndGeneSelector = ({
    timeSeries,
    selectedTimeSeries,
    isFetching,
    connectedTimeSeriesSelected,
}: PropsFromRedux): ReactElement => {
    useReport(
        (processFile) => {
            processFile(
                'Collection/selectedCollection.tsv',
                selectedTimeSeries != null
                    ? objectsArrayToTsv([
                          {
                              id: selectedTimeSeries.id,
                              created: selectedTimeSeries.created,
                              modified: selectedTimeSeries.modified,
                              name: selectedTimeSeries.collection.name,
                              contributor_username: selectedTimeSeries.contributor.username,
                          },
                      ])
                    : '',
                false,
            );
        },
        [selectedTimeSeries],
    );

    const onRowSelectedHandler = (newSelectedTimeSeries: Relation): void => {
        if (selectedTimeSeries == null || newSelectedTimeSeries.id !== selectedTimeSeries.id) {
            connectedTimeSeriesSelected(newSelectedTimeSeries.id);
        }
    };

    return (
        <TimeSeriesAndGeneSelectorContainer>
            <TimeSeriesSelectorWrapper>
                <TimeSeriesSelector
                    timeSeries={timeSeries}
                    selectedTimeSeries={selectedTimeSeries != null ? [selectedTimeSeries] : []}
                    onRowSelected={onRowSelectedHandler}
                    isFetching={isFetching}
                />
            </TimeSeriesSelectorWrapper>
            <GeneSelector />
        </TimeSeriesAndGeneSelectorContainer>
    );
};

export default connector(TimeSeriesAndGeneSelector);
