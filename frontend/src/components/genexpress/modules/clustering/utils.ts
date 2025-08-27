import { groupBy, reverse, round, sum } from 'lodash';
import { agnes, indices, Cluster } from './cluster';
import { ClusteringLinkageFunction, DistanceMeasure } from 'components/genexpress/common/constants';
import { ClusterNode, GeneExpression } from 'redux/models/internal';

export const euclidean = (v1: number[], v2: number[]) => {
    let powSum = 0;
    for (let index = 0; index < v1.length; index++) {
        const difference = v1[index] - v2[index];
        powSum += difference * difference;
    }
    return Math.sqrt(powSum);
};

export const pearson = (v1: number[], v2: number[]) => {
    const sum1 = sum(v1);
    const sum2 = sum(v2);
    const dat1 = v1.map((dataPoint) => Math.pow(dataPoint, 2));
    const dat2 = v2.map((dataPoint) => Math.pow(dataPoint, 2));
    const sum1sq = sum(dat1);
    const sum2sq = sum(dat2);

    const pData = v1.map((dataPoint, index) => dataPoint * v2[index]);

    const pSum = sum(pData);

    const num = pSum - (sum1 * sum2) / v1.length;
    const den = Math.sqrt(
        (sum1sq - Math.pow(sum1, 2) / v1.length) * (sum2sq - Math.pow(sum2, 2) / v1.length),
    );

    if (den === 0) {
        return 0;
    }

    return round(1.0 - num / den, 5);
};

const getRankArray = (arr: number[]) => {
    const sortedArr = [...arr].sort((a, b) => a - b);
    const rankMap = new Map<number, number>();

    for (let i = 0; i < sortedArr.length; i++) {
        const num = sortedArr[i];
        if (!rankMap.has(num)) {
            rankMap.set(num, i + 1);
        }
    }

    // eslint-disable-next-line @typescript-eslint/no-non-null-assertion
    return arr.map((num) => rankMap.get(num)!);
};

const spearman = (arr1: number[], arr2: number[]) => {
    const n = arr1.length;

    const sortedIndexArr1 = getRankArray(arr1);
    const sortedIndexArr2 = getRankArray(arr2);

    let sumOfSquaredRankDifferences = 0;
    for (let i = 0; i < n; i++) {
        const rankDiff = sortedIndexArr1[i] - sortedIndexArr2[i];
        sumOfSquaredRankDifferences += rankDiff * rankDiff;
    }

    const spearmanDistance = (6 * sumOfSquaredRankDifferences) / (n * (n * n - 1));

    return spearmanDistance;
};

const distanceMeasureFunctionMap = {
    [DistanceMeasure.euclidean]: euclidean,
    [DistanceMeasure.pearson]: pearson,
    [DistanceMeasure.spearman]: spearman,
};

export const clusterByGenes = (
    sourceGeneExpressions: GeneExpression[],
    distanceMeasure: DistanceMeasure,
    linkageFunction: ClusteringLinkageFunction,
) => {
    const dataGroupedByGeneId = groupBy(sourceGeneExpressions, 'geneId');

    const genesWithoutValues: string[] = [];
    const genesWithAllConstantValues: string[] = [];
    const genesWithValues: string[] = [];
    const geneExpressions: number[][] = [];

    for (let index = 0; index < Object.keys(dataGroupedByGeneId).length; index++) {
        const geneId = Object.keys(dataGroupedByGeneId)[index];
        const items = dataGroupedByGeneId[geneId];

        if (items.some((item) => !isFinite(item.value)) || items.length === 0) {
            genesWithoutValues.push(geneId);
        } else if (items.every((item) => item.value === items[0].value)) {
            genesWithAllConstantValues.push(geneId);
        } else {
            genesWithValues.push(geneId);
            geneExpressions.push(items.map((item) => item.value));
        }
    }

    if (geneExpressions.length === 0) {
        return {
            withoutValues: genesWithoutValues,
            allConstantValues: genesWithAllConstantValues,
        };
    }

    const tree = agnes(geneExpressions, {
        method: linkageFunction,
        distanceFunction: distanceMeasureFunctionMap[distanceMeasure],
    });
    const intermediateTable: ClusterNode[] = [];

    // First element has order=0, that's why order value needs to be inverted
    // so that first element is shown at the top.
    const leafNodesOrder = reverse(indices(tree));
    // Get max distance from all nodes.
    const max = tree.height * 100;

    // Each clusterNode has to have a nodeIndex in order to track which gene is highlighted.
    let nonLeafNodeIndex = leafNodesOrder.length;

    const processNode = (node: Cluster, parent?: ClusterNode) => {
        const geneId = Object.keys(dataGroupedByGeneId)[node.index];
        if (node.isLeaf) {
            const yCoordinate = leafNodesOrder.indexOf(node.index);
            intermediateTable.push({
                nodeIndex: node.index,
                geneName: dataGroupedByGeneId[geneId][0].geneName,
                geneId: dataGroupedByGeneId[geneId][0].geneId,
                expressions: dataGroupedByGeneId[geneId],
                x: max,
                y: yCoordinate,
                parent,
            });
            return yCoordinate;
        } else if (node.children.length > 0) {
            const parentNode: ClusterNode = {
                x: max - node.height * 100,
                y: 0,
                nodeIndex: nonLeafNodeIndex++,
                parent,
            };
            intermediateTable.push(parentNode);

            let middle = 0;
            node.children.forEach((child) => {
                middle += processNode(child, parentNode);
            });

            middle = middle / 2;

            parentNode.y = middle;

            return middle;
        }

        return 0;
    };

    processNode(tree);

    const clusterNodes = intermediateTable.filter((node) => node.parent != null);
    const order = [
        ...leafNodesOrder.map((index) => genesWithValues[index]),
        ...genesWithoutValues,
        ...genesWithAllConstantValues,
    ];

    return {
        clusterNodes,
        order,
        withoutValues: genesWithoutValues,
        allConstantValues: genesWithAllConstantValues,
    };
};
