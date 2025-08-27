import JSZip from 'jszip';
import _ from 'lodash';
import { v4 as uuidv4 } from 'uuid';
import { saveAs } from 'utils/documentHelpers';
import { logError } from 'utils/errorUtils';

export type ExportFile = {
    path: string;
    content: string | Blob;
    base64: boolean;
};
export type ProcessExportFile = (
    path: string,
    content: string | Blob | null,
    base64: boolean,
) => void;

export type GetComponentReport = (processExportFile: ProcessExportFile) => Promise<void> | void;

const registeredComponents: {
    [id: string]: GetComponentReport;
} = {};

export const reportName = 'Report.zip';

export const register = (getComponentReport: GetComponentReport): (() => void) => {
    const id = uuidv4();
    registeredComponents[id] = getComponentReport;

    return (): void => {
        delete registeredComponents[id];
    };
};

const saveJsZip = async (zip: JSZip, prefix: string): Promise<void> => {
    // No compression.
    const zipContent = await zip.generateAsync({ type: 'blob' });
    saveAs(zipContent, prefix !== '' ? `${prefix}_${reportName}` : reportName, zipContent.type);
};

export const getRegisteredComponentsExportFiles = async (): Promise<ExportFile[]> => {
    const componentsReportsPromises = _.map(registeredComponents, async (getReport) => {
        const files: ExportFile[] = [];
        const processExportFile = (
            path: string,
            content: string | Blob | null,
            base64: boolean,
        ): void => {
            if (!path) {
                throw new Error('Missing file path');
            }
            if (content != null) {
                files.push({ path, content, base64 });
            }
        };

        await getReport(processExportFile);
        return files;
    });

    const arraysOfFiles = await Promise.all(componentsReportsPromises);
    return _.flatten(arraysOfFiles);
};

const createJsZip = (uniqueFiles: ExportFile[]): JSZip => {
    const zip = new JSZip();
    _.each(uniqueFiles, (file) => {
        zip.file(file.path, file.content, { base64: file.base64 });
    });

    return zip;
};

export const exportToZip = async (prefix: string): Promise<void> => {
    const files = await getRegisteredComponentsExportFiles();

    const uniqueFiles = _.uniqBy(files, 'path');
    if (uniqueFiles.length < files.length) {
        logError('All exported files must have unique paths.');
    }

    return saveJsZip(createJsZip(uniqueFiles), prefix);
};
