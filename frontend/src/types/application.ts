import { DescriptorSchemaSlug } from 'components/genexpress/common/constants';

export type UrlDescriptor = {
    name: string;
    url: string;
};
export type Descriptor = {
    [DescriptorSchemaSlug.DictyTimeSeries]: {
        project?: string;
        citation?: UrlDescriptor;
        details?: string;
        strain?: string;
        growth?: string;
        treatment?: string;
    };
};
