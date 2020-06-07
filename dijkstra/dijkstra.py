"""
Dijkstra Algorithm

@ Author: Res. Assist. Ömer Akın
@ Instituiton: Istanbul Technical University Geomatics Engineering Departmant
@ e-mail: akinom@itu.edu.tr

"""


import os, sys, math
import matplotlib.pyplot as plt


def getInput():
    """
    Generate points with corresponding coordinates and graph elements
(nodes, edges, weights) from the input project files

    Returns
    -------
    points : dict
        points with their coordinates
    nodes : set
        set of nodes
    edges : dict
        nodes and their connections
    weights : dict
        distances between connected points
    """

    global file
    while True:
        file = input("\nPress q & Enter anywhere on the program to exit\n\
Please write the project name without extension: ")
        # Exit Conditions
        if file == "q":
            sys.exit()
        # Check file existence
        if (os.path.exists(file + ".xyz") and os.path.exists(file + ".nno")):
            print("Project exists in the directory")
            try:
                pointfile = open(file + ".xyz")
                polylinefile = open(file + ".nno")
                pointlines = pointfile.readlines()[2:]
                polylines = polylinefile.readlines()[2:]
                pointfile.close()
                polylinefile.close()
                pointIDs = [i.split()[0] for i in pointlines]
                x = [float(i.split()[1]) for i in pointlines]
                y = [float(i.split()[2]) for i in pointlines]
                startIDs = [i.split()[0] for i in polylines]
                endIDs = [i.split()[1] for i in polylines]
                points = {point: (i, j) for point, i, j in zip(
                    pointIDs, x, y
                    )}
                nodes = set(points.keys())
                edges = {}
                for node in nodes:
                    endIDlist = []
                    for startID, endID in zip(startIDs, endIDs):
                        if startID == node:
                            endIDlist.append(endID)
                        if endID == node:
                            endIDlist.append(startID)
                        edges[node] = endIDlist
                weights = {
                    (i, j): math.sqrt(
                        (points[i][0]-points[j][0])**2 +
                        (points[i][1]-points[j][1])**2
                        )
                    for i, j in zip(startIDs, endIDs)
                    }
                return points, nodes, edges, weights

            except:
                print("File is not properly designed.")

        else:
            print("File does not exist in the directory")


def origdestInput(nodes):
    """
    Get origin and destination point numbers from the user

    Parameters
    ----------
    nodes : set
        set of nodes

    Returns
    -------
    origin : str
        user inputted origin point
    destination : str
        user inputted destination point
    """

    while True:
        origin = input("Please enter the origin point ID: ")
        if origin == "q":
            print("Exiting.")
            sys.exit()
        elif origin not in nodes:
            print("Selected ID is not found in given file.")
        else:
            break
    while True:
        destination = input("Please enter the destination point ID: ")
        if destination == "q":
            print("Exiting.")
            sys.exit()
        elif destination not in nodes:
            print("Selected ID is not found in given file.")
        else:
            break
    return origin, destination


def dijkstra(nodes, edges, weights, origin):
    """
    Dijkstra algorithm for finding the shortest weights of the accessible
points and path (for an undirected graph)

    Parameters
    ----------
    nodes : set
        set of nodes
    edges : dict
        nodes and their connections
    weights : dict
        distances between connected points
    origin : str
        user inputted origin point

    Returns
    -------
    visited : dict (pointID: weight)
        nodes that are accessible from the origin point and their weights
    path : dict (childnode: parentnode)
        child & parent relationship of accessible nodes
    """

    # Initialize dijkstra algorithm
    visited = {origin: 0}  # {point ID: weight}
    path = {}  # {child: parent}

    unvisited = set(nodes)
    while unvisited:
        nearest_node = None
        for node in unvisited:
            if node in visited:
                if nearest_node is None:
                    nearest_node = node
                # Finding the edge that has minimum weight
                elif visited[node] < visited[nearest_node]:
                    nearest_node = node
        if nearest_node is None:
            break
        unvisited.remove(nearest_node)

        weight_ = visited[nearest_node]
        for edge in edges[nearest_node]:

            # Undirected Structure
            try:
                weight = weight_ + weights[(edge, nearest_node)]
            except:
                pass

            try:
                weight = weight_ + weights[(nearest_node, edge)]
            except:
                pass

            # Update the tables
            if edge not in visited or weight < visited[edge]:
                visited[edge] = weight
                path[edge] = nearest_node
    return visited, path


def shortestPath(path, origin, destination):
    """
    Find the shortest path between origin and destination

    Parameters
    ----------
    path : dict
        child & parent relationship of accessible nodes
    origin : str
        user inputted origin point
    destination : str
        user inputted destination point

    Returns
    -------
    shortest_path : list
        list contains the visited nodes between origin and destination points
    """

    # Check the destination is accessible or not
    if destination not in path:
        print("There is no path between point {} and {}".format(
            origin, destination)
            )
        return None

    shortest_path = []
    origcheck = destination
    while origcheck != origin:
        shortest_path.append(path[origcheck])
        origcheck = path[origcheck]

    shortest_path.reverse()
    shortest_path.append(destination)
    return shortest_path


def plotShortestPath(points, edges, shortestPath):
    """
    Plot the points, lines and shortest path between user entered origin and
destination

    Parameters
    ----------
    points : dict
        points with their coordinates
    edges : dict
        nodes and their connections
    shortestPath : list
        list contains the visited nodes between origin and destination points

    Returns
    -------
    None.
    """

    fig, ax = plt.subplots(figsize=(12, 8))
    for key, value in edges.items():
        for v in value:
            ax.plot([points[key][1], points[v][1]],
                    [points[key][0], points[v][0]], "black", alpha=0.5)

    # Write PointIDs
    for ID, (x, y) in points.items():
        ax.annotate(ID, (y, x),
                    verticalalignment="bottom",
                    horizontalalignment="right")

    for i, point in enumerate(shortestPath):
        if i == len(shortestPath) - 1:
            break

        ax.plot([points[point][1], points[shortestPath[i+1]][1]],
                [points[point][0], points[shortestPath[i+1]][0]], "r-")

    plt.title("Dijkstra Algorithm")
    plt.xlabel("Y")
    plt.ylabel("X")
    plt.show()
    return


def writeOutput(shortestPath, visited):
    """
    Write the result in an output report file in the directory

    Parameters
    ----------
    shortestPath : list
        list contains the visited nodes between origin and destination points
    visited : dict
        nodes that are accessible from the origin point and their weights

    Returns
    -------
    None.
    """

    outputfile = open(file + ".out", "w")
    outputfile.write("Starting Point: {}\nEnd Point: {}\n".format(
        shortestPath[0], shortestPath[-1])
        )
    for point in shortestPath:
        outputfile.write("Point {}: {} Total Distance From Start: {:.2f} m\n\
".format(shortestPath.index(point) + 1, point, visited[point]))

    print("Output report is generated as {}.out in the directory".format(file))
    outputfile.close()
    return


def main():
    points, nodes, edges, weights = getInput()
    origin, destination = origdestInput(nodes)
    visited_points, path = dijkstra(nodes, edges, weights, origin)
    shortest_path = shortestPath(path, origin, destination)
    if shortest_path:
        print("Path between point {} and {} is:\n\
{}".format(origin, destination, "------".join(shortest_path)))
        plotShortestPath(points, edges, shortest_path)
        writeOutput(shortest_path, visited_points)


main()
