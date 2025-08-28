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

export interface UmapPoint {
    cellId: string;
    umap_1: number;
    umap_2: number;
    expressionValue?: number; // For gene-based coloring
}

export interface SingleCellPoint {
    cellId: string;
    umap_1: number;
    umap_2: number;
    expressions: Record<string, number>; // geneId -> expression value
}

export interface UmapStorageData {
    id: number;
    name: string;
    slug: string;
    created: string;
    modified: string;
    contributor: {
        id: number;
        username: string;
        first_name: string;
        last_name: string;
    };
    data: number[];
    json: {
        cells: Record<string, { umap_1: number; umap_2: number }>;
        metadata: {
            n_cells: number;
            umap_key: string;
        };
    };
}

export interface SingleCellStorageData {
    id: number;
    name: string;
    slug: string;
    created: string;
    modified: string;
    contributor: {
        id: number;
        username: string;
        first_name: string;
        last_name: string;
    };
    data: number[];
    json: {
        cells: Record<string, {
            umap_1: number;
            umap_2: number;
        }>;
        genes: Record<string, string>; // geneId -> geneName
        metadata: {
            n_cells: number;
            n_genes: number;
            umap_key: string;
            storage_pattern: string;
        };
    };
}

export interface GeneExpressionsStorageData {
    id: number;
    name: string;
    slug: string;
    created: string;
    modified: string;
    contributor: {
        id: number;
        username: string;
        first_name: string;
        last_name: string;
    };
    data: number[];
    json: Record<string, { // geneId as key
        cells: Record<string, number>; // cellId -> expression value
        metadata: {
            gene_name: string;
            n_cells_with_expression: number;
            max_expression: number;
            min_expression: number;
        };
    }>;
}
