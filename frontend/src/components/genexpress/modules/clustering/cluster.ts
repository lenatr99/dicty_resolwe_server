import getDistanceMatrix from 'ml-distance-matrix';
import { Matrix } from 'ml-matrix';
import { pearson } from './utils';
import { ClusteringLinkageFunction } from 'components/genexpress/common/constants';

// Inspired by https://github.com/mljs/hclust.
export type Cluster = {
    children: Cluster[];
    height: number;
    size: number;
    index: number;
    isLeaf: boolean;
};

/**
 * Traverses the tree depth-first and calls the provided callback with each individual node.
 * @param {Cluster} cluster - Cluster to traverse.
 * @param {function} callback - The callback to be called on each node encounter.
 */
export const traverse = (cluster: Cluster, callback: (node: Cluster) => void) => {
    callback(cluster);
    if (cluster.children) {
        for (const child of cluster.children) {
            traverse(child, callback);
        }
    }
};

/**
 * Returns a list of indices for all the leaves of this cluster.
 * The list is ordered in such a way that a dendrogram could be drawn without crossing branches.
 * @param {Cluster} cluster - Root cluster.
 * @returns {Array<number>}
 */
export const indices = (cluster: Cluster) => {
    const result: number[] = [];
    traverse(cluster, (node) => {
        if (node.isLeaf) {
            result.push(node.index);
        }
    });
    return result;
};

const completeLink = (dKI: number, dKJ: number) => {
    return Math.max(dKI, dKJ);
};

const averageLink = (dKI: number, dKJ: number, dIJ: number, ni: number, nj: number) => {
    const ai = ni / (ni + nj);
    const aj = nj / (ni + nj);
    return ai * dKI + aj * dKJ;
};

const singleLink = (dKI: number, dKJ: number) => {
    return Math.min(dKI, dKJ);
};

const getSmallestDistance = (distance: Matrix) => {
    let smallest = Infinity;
    let smallestI = 0;
    let smallestJ = 0;
    for (let i = 1; i < distance.rows; i++) {
        for (let j = 0; j < i; j++) {
            if (distance.get(i, j) < smallest) {
                smallest = distance.get(i, j);
                smallestI = i;
                smallestJ = j;
            }
        }
    }
    return [smallestI, smallestJ, smallest];
};

const getPreviousIndex = (newIndex: number, prev1: number, prev2: number) => {
    let updatedIndex = newIndex - 1;
    if (updatedIndex >= prev1) {
        updatedIndex++;
    }
    if (updatedIndex >= prev2) {
        updatedIndex++;
    }
    return updatedIndex;
};

/**
 * Continuously merge nodes that have the least dissimilarity.
 * @param {Array<Array<number>>} data - Array of points to be clustered.
 * @param {object} [options]
 * @param {Function} [options.distanceFunction] - Distance function between two number arrays.
 * @param {string} [options.method] - Default: `'average'`.
 */
export const agnes = (
    data: number[][],
    {
        distanceFunction = pearson,
        method = ClusteringLinkageFunction.average,
    }: {
        distanceFunction?: (v1: number[], v2: number[]) => number;
        method: ClusteringLinkageFunction;
    },
) => {
    let updateFunc: typeof averageLink | typeof completeLink = averageLink;
    let distanceMatrix = new Matrix(getDistanceMatrix(data, distanceFunction));
    const numLeaves = distanceMatrix.rows;

    if (typeof method === 'string') {
        switch (method.toLowerCase()) {
            case ClusteringLinkageFunction.complete.toString():
                updateFunc = completeLink;
                break;
            case ClusteringLinkageFunction.average.toString():
                updateFunc = averageLink;
                break;
            case ClusteringLinkageFunction.single.toString():
                updateFunc = singleLink;
                break;
            default:
                throw new RangeError(`unknown clustering method: ${method}`);
        }
    }

    let clusters: Cluster[] = [];
    for (let i = 0; i < numLeaves; i++) {
        clusters.push({ isLeaf: true, index: i, height: 0, children: [], size: 1 });
    }

    for (let n = 0; n < numLeaves - 1; n++) {
        const [row, column, distance] = getSmallestDistance(distanceMatrix);
        const cluster1 = clusters[row];
        const cluster2 = clusters[column];
        const newCluster: Cluster = {
            children: [cluster1, cluster2],
            height: distance,
            isLeaf: false,
            index: -1,
            size: cluster1.size + cluster2.size,
        };

        const newClusters = [newCluster];
        const newDistanceMatrix = new Matrix(distanceMatrix.rows - 1, distanceMatrix.rows - 1);
        const previous = (newIndex: number) =>
            getPreviousIndex(newIndex, Math.min(row, column), Math.max(row, column));

        for (let i = 1; i < newDistanceMatrix.rows; i++) {
            const prevI = previous(i);
            newClusters.push(clusters[prevI]);
            for (let j = 0; j < i; j++) {
                if (j === 0) {
                    const dKI = distanceMatrix.get(row, prevI);
                    const dKJ = distanceMatrix.get(prevI, column);
                    const val = updateFunc(dKI, dKJ, distance, cluster1.size, cluster2.size);
                    newDistanceMatrix.set(i, j, val);
                    newDistanceMatrix.set(j, i, val);
                } else {
                    const val = distanceMatrix.get(prevI, previous(j));
                    newDistanceMatrix.set(i, j, val);
                    newDistanceMatrix.set(j, i, val);
                }
            }
        }

        clusters = newClusters;
        distanceMatrix = newDistanceMatrix;
    }

    return clusters[0];
};
