import { useEffect } from 'react';
import * as reportBuilder from 'components/genexpress/common/reportBuilder/reportBuilder';
import { GetComponentReport } from 'components/genexpress/common/reportBuilder/reportBuilder';

const useReport = (getComponentReport: GetComponentReport, dependencies: unknown[]): void => {
    useEffect(() => {
        return reportBuilder.register(getComponentReport);
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, dependencies);
};

export default useReport;
