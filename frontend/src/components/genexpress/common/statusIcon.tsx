import type * as TIcons from '@mui/icons-material';
import * as Icons from '@mui/icons-material';
import {
    DataStatus,
    DIRTY_DATA_STATUS,
    DONE_DATA_STATUS,
    ERROR_DATA_STATUS,
    PREPARING_DATA_STATUS,
    PROCESSING_DATA_STATUS,
    RESOLVING_DATA_STATUS,
    UPLOADING_DATA_STATUS,
    WAITING_DATA_STATUS,
} from '@genialis/resolwe/dist/api/types/rest';
import { Tooltip } from '@mui/material';

export const StatusIconMap = {
    [UPLOADING_DATA_STATUS]: 'UploadFile',
    [RESOLVING_DATA_STATUS]: 'Schedule',
    [PREPARING_DATA_STATUS]: 'Archive',
    [WAITING_DATA_STATUS]: 'Schedule',
    [PROCESSING_DATA_STATUS]: 'DataUsage',
    [DONE_DATA_STATUS]: 'DoneAll',
    [ERROR_DATA_STATUS]: 'Error',
    [DIRTY_DATA_STATUS]: 'BugReport',
};

export const StatusLabelMap = {
    [UPLOADING_DATA_STATUS]: 'Uploading',
    [RESOLVING_DATA_STATUS]: 'Waiting in queue',
    [PREPARING_DATA_STATUS]: 'Preparing',
    [WAITING_DATA_STATUS]: 'Waiting in queue',
    [PROCESSING_DATA_STATUS]: 'Processing',
    [DONE_DATA_STATUS]: 'Done',
    [ERROR_DATA_STATUS]: 'Error',
    [DIRTY_DATA_STATUS]: 'Dirty',
};

export const StatusIcon = ({ status }: { status: DataStatus }) => {
    const IconComponentName = StatusIconMap[status];
    const Icon = Icons[IconComponentName as keyof typeof TIcons];

    if (status === DONE_DATA_STATUS) {
        return null;
    }

    return (
        <Tooltip title={StatusLabelMap[status]}>
            <Icon fontSize="small" />
        </Tooltip>
    );
};
