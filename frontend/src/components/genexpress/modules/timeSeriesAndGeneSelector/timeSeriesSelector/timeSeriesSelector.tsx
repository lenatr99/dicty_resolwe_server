import React, { ReactElement, useEffect, useState } from 'react';
import { DescriptorSchema, Relation } from '@genialis/resolwe/dist/api/types/rest';
import { useDispatch } from 'react-redux';
import { ColDef } from 'ag-grid-community';
import _ from 'lodash';
import { TimeSeriesSelectorContainer } from './timeSeriesSelector.styles';
import CitationCell from './citationCell/citationCell';
import { getDictyDescriptorSchema } from 'api/descriptorSchemaApi';
import { handleError } from 'utils/errorUtils';
import useStateWithEffect from 'components/genexpress/common/useStateWithEffect';
import DictyGrid, { DictyGridProps } from 'components/genexpress/common/dictyGrid/dictyGrid';

const selectionCell = {
    headerCheckboxSelection: true,
    checkboxSelection: true,
    width: 25,
    field: 'selection',
};

type TimeSeriesSelectorProps = {
    timeSeries: Relation[];
    selectedTimeSeries: Relation[];
    selectionMode?: DictyGridProps<Relation>['selectionMode'];
    onRowSelected?: DictyGridProps<Relation>['onRowSelected'];
    onSelectionChanged?: DictyGridProps<Relation>['onSelectionChanged'];
    isFetching?: boolean;
};

const descriptorFieldsOrder = [
    'selection',
    'descriptor.project',
    'collection.name',
    'descriptor.details',
    'descriptor.strain',
    'descriptor.treatment',
    'descriptor.growth',
    'descriptor.citation',
];

let descriptorSchemaCache: DescriptorSchema | undefined;

const TimeSeriesSelector = ({
    timeSeries,
    selectedTimeSeries,
    selectionMode = 'single',
    onRowSelected,
    onSelectionChanged,
    isFetching,
}: TimeSeriesSelectorProps): ReactElement => {
    const [descriptorSchema, setDescriptorSchema] = useState<DescriptorSchema | undefined>(
        descriptorSchemaCache,
    );
    const [isFetchingDescriptorSchema, setIsFetchingDescriptorSchema] = useState(false);
    const dispatch = useDispatch();

    useEffect(() => {
        const fetchAndSetDescriptorSchema = async () => {
            try {
                setIsFetchingDescriptorSchema(true);

                descriptorSchemaCache = await getDictyDescriptorSchema();

                if (descriptorSchemaCache != null) {
                    setDescriptorSchema(descriptorSchemaCache);
                }
            } catch (error) {
                dispatch(
                    handleError(
                        'Error fetching time series descriptor schema, limited metadata is displayed in Time series Selection.',
                        error,
                    ),
                );
            } finally {
                setIsFetchingDescriptorSchema(false);
            }
        };

        if (descriptorSchema == null) {
            void fetchAndSetDescriptorSchema();
        }
    }, [descriptorSchema, dispatch]);

    const columnDefs = useStateWithEffect(() => {
        if (descriptorSchema == null) {
            return [
                ...(selectionMode === 'multiple' ? [selectionCell] : []),
                { field: 'collection.name', headerName: 'Name' },
            ];
        }

        return [
            ...(selectionMode === 'multiple' ? [selectionCell] : []),
            ...descriptorSchema.schema.map((fieldSchema, index) => {
                if (fieldSchema.type.includes('url')) {
                    return {
                        field: `descriptor.${fieldSchema.name}`,
                        headerName: fieldSchema.label,
                        cellStyle: { padding: 0 },
                        cellRenderer: CitationCell,
                        valueFormatter: ({ data }) => data.descriptor?.citation?.name,
                        getQuickFilterText: (params) => {
                            return params.data.descriptor?.citation?.name;
                        },
                    } as ColDef;
                }

                return {
                    field: `descriptor.${fieldSchema.name}`,
                    headerName: fieldSchema.label,
                    ...(fieldSchema.name === 'details' && {
                        sort: 'asc',
                        sortIndex: 1,
                    }),
                    ...(index === 0 && {
                        sort: 'asc',
                        sortIndex: 0,
                        valueGetter: ({ data }) =>
                            _.get(data, `descriptor.${fieldSchema.name}`) ??
                            _.get(data, 'collection.name'),
                    }),
                } as ColDef;
            }),
        ].sort((a: ColDef, b: ColDef) => {
            if (descriptorFieldsOrder.indexOf(a.field ?? '') < 0) {
                return 1;
            }
            if (descriptorFieldsOrder.indexOf(b.field ?? '') < 0) {
                return -1;
            }
            return (
                descriptorFieldsOrder.indexOf(a.field ?? '') -
                descriptorFieldsOrder.indexOf(b.field ?? '')
            );
        });
    }, [descriptorSchema, selectionMode]);

    return (
        <TimeSeriesSelectorContainer>
            <DictyGrid
                isFetching={isFetching || isFetchingDescriptorSchema}
                data={timeSeries}
                selectionMode={selectionMode}
                filterLabel="Filter time series"
                disableSizeColumnsToFit
                columnDefs={columnDefs}
                getRowId={(data): string => data.id.toString()}
                onRowSelected={onRowSelected}
                onSelectionChanged={onSelectionChanged}
                selectedData={selectedTimeSeries.length > 0 ? selectedTimeSeries : undefined}
            />
        </TimeSeriesSelectorContainer>
    );
};

export default TimeSeriesSelector;
