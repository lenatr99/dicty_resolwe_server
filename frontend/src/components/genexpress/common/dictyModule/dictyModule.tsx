import { DataStatus } from '@genialis/resolwe/dist/api/types/rest';
import React, { ReactElement, ReactNode } from 'react';
import { StatusIcon } from '../statusIcon';
import { ModuleHeader, ModuleContent, ModuleContainer, LoadingBar } from './dictyModule.styles';

type SelectedGenesProps = {
    children: ReactNode;
    title: string;
    isLoading: boolean;
    status?: DataStatus | null;
};

const DictyModule = ({ children, title, isLoading, status }: SelectedGenesProps): ReactElement => {
    return (
        <ModuleContainer>
            <ModuleHeader className="dragHandle">
                {title}
                {status != null && <StatusIcon status={status} />}
                {isLoading && <LoadingBar color="secondary" />}
            </ModuleHeader>
            <ModuleContent>{children}</ModuleContent>
        </ModuleContainer>
    );
};

export default DictyModule;
