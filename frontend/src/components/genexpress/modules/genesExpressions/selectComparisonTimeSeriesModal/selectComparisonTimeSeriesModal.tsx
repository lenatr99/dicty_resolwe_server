import React, { ReactElement } from 'react';
import { Button } from '@mui/material';
import { connect, ConnectedProps } from 'react-redux';
import _ from 'lodash';
import TimeSeriesSelector from '../../timeSeriesAndGeneSelector/timeSeriesSelector/timeSeriesSelector';
import { TimeSeriesGridWrapper } from './selectComparisonTimeSeriesModal.style';
import {
    CenteredModal,
    ModalBody,
    ModalContainer,
    ModalFooter,
    ModalHeader,
} from 'components/genexpress/common/dictyModal/dictyModal.styles';
import {
    comparisonTimeSeriesChanged,
    getComparisonTimeSeries,
    getSelectedTimeSeries,
    getTimeSeries,
} from 'redux/stores/timeSeries';
import { RootState } from 'redux/rootReducer';
import useStateWithEffect from 'components/genexpress/common/useStateWithEffect';

const mapStateToProps = (state: RootState) => {
    return {
        timeSeries: getTimeSeries(state.timeSeries),
        selectedTimeSeries: getSelectedTimeSeries(state.timeSeries),
        comparisonTimeSeries: getComparisonTimeSeries(state.timeSeries),
    };
};

const connector = connect(mapStateToProps, {
    connectedComparisonTimeSeriesChanged: comparisonTimeSeriesChanged,
});

type PropsFromRedux = ConnectedProps<typeof connector>;

type SelectComparisonTimeSeriesModalProps = {
    handleOnClose: () => void;
} & PropsFromRedux;

const SelectComparisonTimeSeriesModal = ({
    timeSeries,
    selectedTimeSeries,
    comparisonTimeSeries,
    connectedComparisonTimeSeriesChanged,
    handleOnClose,
}: SelectComparisonTimeSeriesModalProps): ReactElement => {
    const filteredTimeSeries = useStateWithEffect(() => {
        return timeSeries.filter(
            (singleTimeSeries) => singleTimeSeries.id !== selectedTimeSeries?.id,
        );
    }, [selectedTimeSeries, timeSeries]);

    return (
        <CenteredModal
            open
            aria-labelledby="modalTitle"
            aria-describedby="modalDescription"
            onClose={handleOnClose}
        >
            <ModalContainer>
                <ModalHeader id="modalTitle">Select Time series To Compare</ModalHeader>
                <ModalBody>
                    <TimeSeriesGridWrapper>
                        <TimeSeriesSelector
                            timeSeries={filteredTimeSeries}
                            selectedTimeSeries={comparisonTimeSeries}
                            selectionMode="multiple"
                            onSelectionChanged={(selected): void => {
                                if (_.xor(selected, comparisonTimeSeries).length > 0) {
                                    connectedComparisonTimeSeriesChanged(
                                        selected.map((singleTimeSeries) => singleTimeSeries.id),
                                    );
                                }
                            }}
                        />
                    </TimeSeriesGridWrapper>
                </ModalBody>
                <ModalFooter>
                    <Button onClick={handleOnClose}>Close</Button>
                </ModalFooter>
            </ModalContainer>
        </CenteredModal>
    );
};

export default connector(SelectComparisonTimeSeriesModal);
