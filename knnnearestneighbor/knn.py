"""
KNN Nearest Neigbor Algorithm - Derived from edX-Python for Research Course

@ Author: Res. Assist. Ömer Akın
@ Instituiton: Istanbul Technical University Geomatics Engineering Departmant
@ e-mail: akinom@itu.edu.tr

"""


import random
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt


def distance(p1, p2):
    """
    Distance between two points

    Parameters
    ----------
    p1 : numpy.ndarray
        first point
    p2 : numpy.ndarray
        second point

    Returns
    -------
    dist : numpy.float64
        distance between two points

    """

    dist = np.sqrt(np.sum(np.power(p2 - p1, 2)))
    return dist


def majorityVote(votes):
    """
    Find the mod of data

    """

    vote_counts = {}
    for vote in votes:
        if vote in vote_counts:
            vote_counts[vote] += 1
        else:
            vote_counts[vote] = 1

    winners = []
    max_value = max(vote_counts.keys())
    max_count = max(vote_counts.values())
    for vote, count in vote_counts.items():
        if count == max_count:
            winners.append(vote)
    winner = random.choice(winners)
    return winner


def majorityVoteFast(votes):
    """
    Faster version of majorityVote function with scipy library

    """

    mode, count = ss.mstats.mode(votes)
    return mode


def findNearestNeighbors(p, points, k=5):
    """
    Find k nearest neighbors of point p in given points

    Parameters
    ----------
    p : numpy array
        points that will be classified
    points : numpy array
        point list
    k : int, optional
        number of neighbors. The default value is 5.

    """

    distances = np.zeros(points.shape[0])
    for i in range(len(distances)):
        distances[i] = distance(p, points[i])
    ind = np.argsort(distances)
    return ind[:k]


def knnPredict(p, points, outcomes, k=5):
    # Find k nearest neighbors
    ind = findNearestNeighbors(p, points, k)
    # Predict the class of P based on majority vote
    return majorityVote(outcomes[ind])


def generateSyntheticData(n=50):
    """
    Create two sets of points from bivariate normal distributions

    """

    points = np.concatenate((ss.norm(0, 1).rvs((n, 2)),
                             ss.norm(1, 1).rvs((n, 2))), axis=0)
    outcomes = np.concatenate((np.repeat(0, n), np.repeat(1, n)))

    plt.figure()
    plt.plot(points[:n, 0], points[:n, 1], "ro")
    plt.plot(points[n:, 0], points[n:, 1], "go")
    plt.savefig("initial_data.png")
    plt.show()
    return (points, outcomes)


def makePredictionGrid(predictors, outcomes, limits, h, k):
    """
    Classify each point on the prediction grid

    """

    (x_min, x_max, y_min, y_max) = limits
    xs = np.arange(x_min, x_max, h)
    ys = np.arange(y_min, y_max, h)
    xx, yy = np.meshgrid(xs, ys)

    prediction_grid = np.zeros(xx.shape, dtype=int)
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            p = np.array([x, y])
            prediction_grid[j, i] = knnPredict(p, predictors, outcomes, k)

    return (xx, yy, prediction_grid)


def plotPredictionGrid(xx, yy, predictors, outcomes, predictionGrid, filename):
    """
    Plot KNN predictions for every point on the grid

    """

    from matplotlib.colors import ListedColormap
    background_colormap = ListedColormap(
        ["hotpink", "lightskyblue", "yellowgreen"]
        )
    observation_colormap = ListedColormap(["red", "blue", "green"])
    plt.figure(figsize=(12, 10))
    plt.pcolormesh(
        xx, yy,
        predictionGrid,
        cmap=background_colormap,
        alpha=0.5)
    plt.scatter(predictors[:, 0],
                predictors[:, 1],
                c=outcomes,
                cmap=observation_colormap,
                s=50)
    plt.title("Classification Result")
    plt.xlabel('Variable 1')
    plt.ylabel('Variable 2')
    plt.xticks(())
    plt.yticks(())
    plt.xlim(np.min(xx), np.max(xx))
    plt.ylim(np.min(yy), np.max(yy))
    plt.savefig(filename)


def main():
    from sklearn import datasets
    iris = datasets.load_iris()
    predictors = iris["data"][:, 0:2]
    outcomes = iris["target"]

    plt.plot(
        predictors[outcomes == 0][:, 0], predictors[outcomes == 0][:, 1], "ro"
        )
    plt.plot(
        predictors[outcomes == 1][:, 0], predictors[outcomes == 1][:, 1], "bo"
        )
    plt.plot(
        predictors[outcomes == 2][:, 0], predictors[outcomes == 2][:, 1], "go"
        )
    plt.title("Data that will be classified")
    plt.savefig("iris.png")

    k = 5
    filename = "iris_grid.pdf"
    limits = (4, 8, 1.5, 4.5)
    h = 0.1
    (xx, yy, prediction_grid) = makePredictionGrid(
        predictors, outcomes, limits, h, k
        )
    plotPredictionGrid(xx, yy, predictors, outcomes, prediction_grid, filename)
    my_predictions = np.array(
        [knnPredict(p, predictors, outcomes, 5) for p in predictors]
        )


if __name__ == "__main__":
    main()
