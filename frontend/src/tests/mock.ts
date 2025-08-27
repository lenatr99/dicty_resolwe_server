import configureMockStore, { MockStoreEnhanced } from 'redux-mock-store';
import {
    Relation,
    Contributor,
    ItemPermissionsOf,
    CollectionPermissions,
    Collection,
    DescriptorSchema,
    RelationPartition,
    User,
    DataGafAnnotation,
    DONE_DATA_STATUS,
    Data,
    Process,
    Storage,
} from '@genialis/resolwe/dist/api/types/rest';
import { v4 as uuidv4 } from 'uuid';
import _ from 'lodash';
import { createEpicMiddleware } from 'redux-observable';
import {
    RelationsById,
    GenesById,
    BasketInfo,
    SnackbarNotifications,
    SnackbarNotification,
    SamplesGenesExpressionsById,
    GeneSet,
    Gene,
    VolcanoPoint,
    DifferentialExpression,
    DifferentialExpressionsById,
    GOEnrichmentRow,
    EnhancedGOEnrichmentJson,
    GenesExpressionById,
    BasketExpression,
    GeneSimilarity,
    BookmarkComponentsState,
} from '../redux/models/internal';
import { BookmarkReduxState, RootState } from '../redux/rootReducer';
import { BackendAppState, BasketAddSamplesResponse } from 'redux/models/rest';
import { AppDispatch } from 'redux/appStore';
import { generateRandomString, generateRandomStrings } from 'utils/stringUtils';
import { generateRandomNumbers } from 'utils/numberUtils';
import { flattenGoEnrichmentTree } from 'utils/gOEnrichmentUtils';
import { DistanceMeasure, DictyUrlQueryParameter } from 'components/genexpress/common/constants';
import { pValueThresholdsOptions } from 'redux/stores/gOEnrichment';

const getDateISOString = (): string => new Date().toISOString();

export const generateContributor = (id: number): Contributor => ({
    id,
    first_name: 'Janez',
    last_name: 'Novak',
    username: 'janeznovak',
});

export const generateDescriptorSchema = (id: number): DescriptorSchema => ({
    id,
    created: getDateISOString(),
    modified: getDateISOString(),
    current_user_permissions: [
        {
            type: 'user',
            permissions: ['view'],
        },
        {
            type: 'group',
            permissions: ['edit'],
        },
    ],
    slug: `testDescriptorSchema-${id}`,
    name: `Test descriptor schema ${id}`,
    version: '1.0',
    schema: [],
    contributor: generateContributor(1),
});

export const generateCollection = (id: number): Collection => ({
    id,
    created: getDateISOString(),
    modified: getDateISOString(),
    duplicated: 'No',
    slug: `testCollection-${id}`,
    name: `Test collection ${id}`,
    description: 'Collection description',
    descriptor: {},
    descriptor_schema: generateDescriptorSchema(3),
    contributor: generateContributor(1),
    current_user_permissions: [] as ItemPermissionsOf<CollectionPermissions>[],
    tags: ['community:bcm'],
    settings: {},
    entity_count: 2,
    data_count: 3,
    status: 'OK',
});

export const generateSingleTimeSeries = (id: number): Relation => ({
    id,
    slug: `relation-${id}`,
    created: getDateISOString(),
    modified: getDateISOString(),
    contributor: generateContributor(1),
    type: 'series',
    collection: generateCollection(id),
    category: 'Time series',
    partitions: [{ id: 1, position: 1, entity: 1, label: 'Test label' }],
    unit: null,
    descriptor: {
        citation: {
            name: 'Rosengarten et. al.',
            url: 'http://www.biomedcentral.com/1471-2164/16/294',
        },
        details: 'Filter development',
        growth: 'HL5',
        project: '2. Filter Development vs. cAMP Pulsing; Frequent Sampling',
        strain: 'AX4',
        treatment: 'Filter development',
    },
    descriptor_schema: generateDescriptorSchema(3),
});

export const generateGene = (
    id: number | string,
    source = 'ENSEMBL',
    species = 'Homo sapiens',
): Gene => ({
    feature_id: id.toString(),
    aliases: [`TestGene${id}_alias1`, `TestGene${id}_alias2`],
    name: `TestGene${id}`,
    full_name: `Test gene ${id}`,
    description: `This is a test gene ${id}`,
    source,
    species,
});

export const generateGenes = (genesIds: string[], source?: string, species?: string): Gene[] => {
    return genesIds.map((geneId) => generateGene(geneId, source, species));
};

export const generateNotification = (key: string | number): SnackbarNotification => ({
    key,
    message: 'Test message',
    variant: 'info',
});

export const generateBasketAddSamplesResponse = (id: string): BasketAddSamplesResponse => ({
    id,
    modified: getDateISOString(),
    conflict_organisms: [],
    conflict_sources: [],
    ignored: [],
    permitted_organisms: ['Homo sapiens'],
    permitted_sources: ['UCSC'],
});

export const generateBasketInfo = (id: string): BasketInfo => ({
    id,
    source: 'DICTYBASE',
    species: 'Dictyostelium purpureum',
    type: 'gene',
});

export const generateGenesExpressionsById = (genesIds: string[]): GenesExpressionById => {
    return genesIds.reduce(
        (byId, geneId) => ({
            ...byId,
            [geneId]: generateRandomNumbers(1, () => Math.random() * 10)[0],
        }),
        {} as GenesExpressionById,
    );
};

const generateDifferentialExpressionFile = (): { file: string; size: number } => ({
    file: `${generateRandomString(5)}.tab.gz`,
    size: 1234,
});

export const generateRelationPartition = (id: number, entityId: number): RelationPartition => ({
    id,
    entity: entityId,
    position: 0,
    label: `${id}Hr`,
});

export const generateRelationPartitions = (entityIds: number[]): RelationPartition[] => {
    return _.times(entityIds.length, (i) => generateRelationPartition(i, entityIds[i]));
};

const generateProcess = (id: number): Process => ({
    id: 1,
    slug: 'clustering-hierarchical-samples',
    name: 'Hierarchical clustering of samples',
    created: getDateISOString(),
    modified: getDateISOString(),
    version: '2.0.0',
    type: 'data:clustering:hierarchical:sample:',
    category: 'Analyses:',
    persistence: 'TMP',
    description: 'Hierarchical clustering of samples.',
    input_schema: [
        {
            name: 'processing',
            group: [
                {
                    default: 'spearman',
                    label: 'Distance metric',
                    type: 'basic:string:',
                    name: 'distance_metric',
                    choices: [
                        {
                            value: 'spearman',
                            label: 'spearman',
                        },
                        {
                            value: 'pearson',
                            label: 'pearson',
                        },
                        {
                            value: 'euclidean',
                            label: 'euclidean',
                        },
                    ],
                },
                {
                    default: 'average',
                    label: 'Linkage method',
                    type: 'basic:string:',
                    name: 'linkage_method',
                    choices: [
                        {
                            value: 'average',
                            label: 'average',
                        },
                        {
                            value: 'single',
                            label: 'single',
                        },
                        {
                            value: 'complete',
                            label: 'complete',
                        },
                    ],
                },
            ],
        },
    ],
    output_schema: null,
    run: null,
    contributor: generateContributor(id),
    current_user_permissions: [],
    is_active: true,
    data_name: '',
    entity_descriptor_schema: undefined,
    entity_input: '',
    entity_type: undefined,
    scheduling_class: 'BA',
    entity_always_create: true,
});

export const generateData = (id: number): Data => ({
    id,
    created: getDateISOString(),
    modified: getDateISOString(),
    contributor: {
        id: 1,
        first_name: '',
        last_name: '',
        username: '',
    },
    current_user_permissions: [],
    started: getDateISOString(),
    finished: getDateISOString(),
    scheduled: getDateISOString(),
    duplicated: undefined,
    checksum: '',
    status: DONE_DATA_STATUS,
    process_progress: 100,
    process_rc: 0,
    size: 123,
    process_memory: 123,
    process_cores: 2,
    process_info: [],
    process_warning: [],
    process_error: [],
    process: generateProcess(id),
    output: {},
    slug: `test-data-${id}`,
    name: `Test Data ${id}`,
    descriptor_schema: {
        id: 1,
        slug: 'test-schema',
        name: 'Test Schema',
        version: '1.0.0',
        contributor: generateContributor(id),
        schema: [],
        created: getDateISOString(),
        modified: getDateISOString(),
    },
    descriptor: {},
    input: {},
    tags: [],
    entity: undefined,
    collection: undefined,
});

export const generateGeneExpression = (geneId: number): GenesExpressionById => ({
    [geneId.toString()]: new Date().getTime(),
});

export const generateDifferentialExpressionJson = (
    numberOfGenes = 10,
): {
    fdr: number[];
    gene_id: string[];
    logfc: number[];
    pvalue: number[];
    stat: number[];
} => ({
    fdr: generateRandomNumbers(numberOfGenes, () => Math.random() * 10),
    gene_id: generateRandomStrings(numberOfGenes, 8),
    logfc: generateRandomNumbers(numberOfGenes, () => Math.random() * 10),
    pvalue: generateRandomNumbers(numberOfGenes, () => Math.random() * 10),
    stat: generateRandomNumbers(numberOfGenes, () => Math.random() * 10),
});

export const generateDifferentialExpression = (id: number): DifferentialExpression => ({
    ...generateData(id),
    name: `TestDifferentialExpression${id}`,
    status: DONE_DATA_STATUS,
    logfc_threshold: 1,
    // Has to be the same as property in differentialExpressionJson.
    prob_field: 'fdr',
    prob_threshold: 0.05,
    up_regulated: Math.random() * 10,
    down_regulated: Math.random() * 10,
    output: {
        de_json: Math.random() * 10,
        de_file: generateDifferentialExpressionFile(),
        raw: generateDifferentialExpressionFile(),
        source: 'DICTYBASE',
        species: 'Dictyostelium discoideum',
    },
    slug: `testDifferentialExpression${id}`,
});

export const generateVolcanoPoint = (geneId: number): VolcanoPoint => ({
    geneId: geneId.toString(),
    geneName: geneId.toString(),
    logFcValue: Math.random(),
    logProbFiniteValue: Math.random(),
    logProbValue: Math.random(),
    probValue: Math.random(),
});

export const generateVolcanoPoints = (n: number): VolcanoPoint[] => {
    return _.times(n, (i) => {
        return generateVolcanoPoint(i + 1);
    });
};

export const generateGeneSet = (): GeneSet => ({
    dateTime: new Date(Math.random() * 1000000000000),
    genesNames: [generateRandomString(3), generateRandomString(3), generateRandomString(3)],
});

export const generateGeneSets = (n: number): GeneSet[] => {
    return _.times(n, () => {
        return generateGeneSet();
    });
};

export const generatePartition = (id: number, entityId: number): RelationPartition => ({
    id,
    entity: entityId,
    position: 0,
    label: `${id}Hr`,
});

export const generateUser = (id: number): User => ({
    id,
    username: 'testUser',
    email: 'test@mail.com',
    first_name: 'Test',
    last_name: 'User',
    job_title: 'Developer',
    company: 'Genialis',
    department: 'Software development',
    location: 'Ljubljana',
    lab: 'Bio',
    phone_number: '031123456',
    last_login: getDateISOString(),
    date_joined: getDateISOString(),
});

export const generateStorage = <T>(id: number, json: T): Storage => ({
    id,
    slug: 'test-storage',
    name: 'Test Storage',
    data: 1,
    json,
    contributor: generateContributor(id),
    created: getDateISOString(),
    modified: getDateISOString(),
});

export const generateGOEnrichmentRow = (
    id: number,
    genesIds = [] as string[],
): GOEnrichmentRow => ({
    gene_ids: genesIds,
    term_name: `Parent-${id}-${generateRandomString(5)}`,
    term_id: `Parent-${id}-${generateRandomString(5)}`,
    pval: 0.05,
    score: 500,
    total: 250,
    matched: 125,
    path: [],
    score_percentage: 0,
    gene_associations: [],
    children: [
        {
            gene_ids: [],
            term_name: `Child-${id}-${generateRandomString(5)}`,
            term_id: `Child-${id}-${generateRandomString(5)}`,
            pval: 0.01,
            score: 250,
            total: 250,
            matched: 125,
            path: [],
            score_percentage: 0,
            gene_associations: [],
            children: [],
        },
    ],
});

export const generateGeneOntologyStorageJson = (genesIds: string[]): EnhancedGOEnrichmentJson => {
    const json: EnhancedGOEnrichmentJson = {
        total_genes: 500,
        gene_associations: {},
        tree: {
            BP: [
                generateGOEnrichmentRow(0, genesIds),
                generateGOEnrichmentRow(1, genesIds),
                generateGOEnrichmentRow(2, genesIds),
                generateGOEnrichmentRow(3, genesIds),
                generateGOEnrichmentRow(4, genesIds),
            ],
            CC: [
                generateGOEnrichmentRow(5, genesIds),
                generateGOEnrichmentRow(6, genesIds),
                generateGOEnrichmentRow(7, genesIds),
            ],
            MF: [
                generateGOEnrichmentRow(8, genesIds),
                generateGOEnrichmentRow(9, genesIds),
                generateGOEnrichmentRow(10, genesIds),
            ],
        },
    };
    const allRows = [
        ...flattenGoEnrichmentTree(json.tree.BP, []),
        ...flattenGoEnrichmentTree(json.tree.CC, []),
        ...flattenGoEnrichmentTree(json.tree.MF, []),
    ];

    allRows.forEach((row): void => {
        json.gene_associations[row.term_id] = genesIds;
    });

    return json;
};

export const generateGaf = (
    id: number,
): {
    humanGaf: DataGafAnnotation;
    mouseMGIGaf: DataGafAnnotation;
    mouseUCSCGaf: DataGafAnnotation;
} => {
    const humanGaf = {
        ...generateData(id),
        name: 'Human gaf',
        slug: 'human-gaf',
        output: {
            source: 'UniProtKB',
            species: 'Homo sapiens',
            gaf: {
                file: 'goa_human.gaf.txt.gz',
                size: 1231231,
            },
            gaf_obj: {
                file: 'gaf_obj',
                size: 1231231,
            },
        },
    };

    const mouseMGIGaf = {
        ...generateData(id),
        name: 'Mouse MGI gaf',
        slug: 'mouse-mgi-gaf',
        output: {
            source: 'MGI',
            species: 'Mus musculus',
            gaf: {
                file: 'gene_association.mgi.txt.gz',
                size: 1231231,
            },
            gaf_obj: {
                file: 'gaf_obj',
                size: 1231231,
            },
        },
    };

    const mouseUCSCGaf = {
        ...generateData(id),
        name: 'Mouse UCSC gaf',
        slug: 'mouse-ucsc-gaf',
        output: {
            source: 'UCSC',
            species: 'Mus musculus',
            gaf: {
                file: 'gene_association.mgi.txt.gz',
                size: 1231231,
            },
            gaf_obj: {
                file: 'gaf_obj',
                size: 1231231,
            },
        },
    };

    return { humanGaf, mouseMGIGaf, mouseUCSCGaf };
};

export const generateBasketExpression = (): BasketExpression => ({
    id: generateRandomNumbers(1, () => Math.random() * 10000)[0],
    exp_type: 'polyA',
});

export const generateGeneSimilarity = (geneId: string): GeneSimilarity => ({
    gene: geneId,
    distance: generateRandomNumbers(1, () => Math.random() * 10)[0],
});

export const generateBackendBookmark = (
    selectedTimeSeriesId?: number,
    selectedGenesIds: string[] = [],
): BackendAppState<BookmarkReduxState & BookmarkComponentsState> => ({
    contributor: 1,
    uuid: uuidv4(),
    state: {
        timeSeries: { selectedId: selectedTimeSeriesId ?? null, comparisonIds: [] },
        genes: {
            selectedGenesIds,
            highlightedGenesIds: [],
            source: 'DICTYBASE',
            species: 'Dictyostelium purpureum',
        },
        differentialExpressions: { selectedId: null },
        gOEnrichment: { pValueThreshold: 0.01 },
        genesSimilarities: { queryGeneId: null, distanceMeasure: DistanceMeasure.pearson },
    },
});

export const generateBookmarkQueryParameter = (): string =>
    `?${DictyUrlQueryParameter.appState}=${uuidv4()}`;
export const generateGenesQueryParameter = (genesIds: string[]): string =>
    `?${DictyUrlQueryParameter.genes}=${genesIds.join(',')}`;

/**
 * Helper function that generates mock instances of objects for use in unit tests.
 *
 * @param n Number of mock instances to generate.
 * @returns Instances in object, accessible by id.
 */
type InstanceById<T> = {
    [_ in string | number]: T;
};

let lastInstanceId = 1;
export const generateInstances = <T>(
    n: number,
    getId: (instance: T, instanceIndex: number) => string | number,
    generateSingleInstance: (id: number) => T,
): InstanceById<T> => {
    const instancesById = {} as InstanceById<T>;

    for (let i = lastInstanceId; i < lastInstanceId + n; i += 1) {
        const instance = generateSingleInstance(i);
        instancesById[getId(instance, i)] = instance;
    }

    lastInstanceId += n;

    return instancesById;
};

export const generateTimeSeriesById = (n: number): RelationsById => {
    return generateInstances<Relation>(n, (instance) => instance.id, generateSingleTimeSeries);
};

export const generateGenesById = (n: number): GenesById => {
    return generateInstances<Gene>(n, (instance) => instance.feature_id, generateGene);
};

export const generateGenesByIdPredefinedIds = (genesIds: string[]): GenesById => {
    return _.keyBy(
        genesIds.map((geneId) => generateGene(geneId)),
        'feature_id',
    );
};

export const generateRelationsById = (n: number): RelationsById => {
    return generateInstances<Relation>(n, (instance) => instance.id, generateSingleTimeSeries);
};

export const generateSamplesExpressionsById = (
    n: number,
    genesIds: string[],
): SamplesGenesExpressionsById => {
    return generateInstances<GenesExpressionById>(
        n,
        (_instance, instanceIndex) => instanceIndex,
        () => generateGenesExpressionsById(genesIds),
    );
};

export const generateDifferentialExpressionsById = (n: number): DifferentialExpressionsById => {
    return generateInstances<DifferentialExpression>(
        n,
        (instance) => instance.id,
        generateDifferentialExpression,
    );
};

export const testState = (): RootState => {
    const timeSeriesById = generateTimeSeriesById(2);

    return {
        layouts: {},
        authentication: {
            user: generateUser(2),
            isLoggedIn: true,
            isLoggingIn: false,
            isLoggingOut: false,
            isFetchingUser: false,
        },
        timeSeries: {
            byId: timeSeriesById,
            selectedId: _.flatMap(timeSeriesById)[0].id,
            comparisonIds: [],
            isFetching: false,
            isAddingToBasket: false,
            isFetchingGenesMappings: false,
            basketInfo: generateBasketInfo('1'),
            basketExpressionsIds: [],
        },
        genes: {
            byId: generateGenesById(2),
            selectedGenesIds: [],
            highlightedGenesIds: [],
            isFetchingDifferentialExpressionGenes: false,
            isFetchingAssociationsGenes: false,
            isFetchingSimilarGenes: false,
            isFetchingBookmarkedGenes: false,
        },
        genesSimilarities: {
            status: null,
            data: [],
            queryGeneId: null,
            distanceMeasure: DistanceMeasure.spearman,
            isFetchingGenesSimilarities: false,
        },
        samplesExpressions: {
            byId: {},
            isFetchingSamplesExpressions: false,
        },
        differentialExpressions: {
            byId: {},
            isFetchingDifferentialExpressions: false,
            isFetchingDifferentialExpressionsData: false,
            selectedId: 0,
        },
        gOEnrichment: {
            status: null,
            json: null,
            gaf: generateGaf(1).humanGaf,
            source: '',
            species: '',
            pValueThreshold: pValueThresholdsOptions[0],
            isFetchingJson: false,
            ontologyObo: { id: 24330 } as Data,
        },
        notifications: { notifications: [] as SnackbarNotifications },
    };
};

export const mockStore = (
    initialTestState: RootState,
): MockStoreEnhanced<RootState, AppDispatch> => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const epicMiddleware = createEpicMiddleware<any, any, RootState, any>();

    const createMockStore = configureMockStore<RootState, AppDispatch>([epicMiddleware]);
    const store = createMockStore(initialTestState);
    store.clearActions();

    return store;
};
