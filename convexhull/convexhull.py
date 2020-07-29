"""
Convex Hull with Graham Scan Algorithm

@ Author: Res. Assist. Ömer Akın
@ Institution: Istanbul Technical University Geomatics Engineering Departmant
@ e-mail: akinom@itu.edu.tr

"""


import os, sys, math, random
import matplotlib.pyplot as plt


def getInput():
    """
    Generate points randomly or from the input project files

    Returns
    -------
    pointIDs : list
        point IDs in order
    points : list
        contains coordinate tuples

    """

    global file
    while True:
        selection = input("Welcome to Convex Hull Finder\nPress q to exit.\n\
To find the convex hull of random points, please press 1\nTo find the convex \
hull of given point file, please press 2: \n")
        if selection == "q":
            sys.exit()
        elif selection == "1":
            pointnumber = input("How many points you want to generate: ")
            if pointnumber == "q":
                sys.exit()
            elif pointnumber.isdigit():
                pointIDs = [i + 1 for i in range(int(pointnumber))]
                x = [random.uniform(20, 60) for _ in range(int(pointnumber))]
                y = [random.uniform(20, 60) for _ in range(int(pointnumber))]
                points = [(i, j) for i, j in zip(x, y)]
                file = "randompoints"
                return pointIDs, points
            else:
                print("Please enter a valid number")
        elif selection == "2":
            file = input("Please enter the file name of points without \
its extension: ")
            # Exit Conditions
            if file == "q":
                sys.exit()
            # Check file existence
            if os.path.exists(file + ".xyz"):
                print("File exists in the directory")
                try:
                    # To add points
                    pointfile = open(file + ".xyz")
                    p_lines = pointfile.readlines()[2:]
                    pointfile.close()
                    pointIDs = [i.split()[0] for i in p_lines]
                    x = [float(i.split()[1]) for i in p_lines]
                    y = [float(i.split()[2]) for i in p_lines]
                    points = [(i, j) for i, j in zip(x, y)]
                    return pointIDs, points

                except:
                    print("However it is not supported or properly designed")

            else:
                print("File does not exist in the directory")

        else:
            print("Please enter a valid selection")


def findAnchor(IDs, points):
    """
    Find the most southwest point to be used as initial start point
    for searching

    Parameters
    ----------
    IDs : list
        point IDs in order
    points : list
        contains coordinate tuples

    Returns
    -------
    anchor : tuple
        Most soutwestern point of given point list

    """

    min_x = None
    for i, (x, y) in zip(IDs, points):
        if min_x is None or x < points[IDs.index(min_x)][0]:
            min_x = i
        if (
                x == points[IDs.index(min_x)][0] and
                y < points[IDs.index(min_x)][1]
                ):
            min_x = i
    anchor = points[IDs.index(min_x)]
    return anchor


def calcAzimuth(p0, p1):
    """
    Calculate azimuthal angle between two points

    Parameters
    ----------
    p0 : tuple
        first point
    p1 : tuple
        second point

    Returns
    -------
    azimuth : float
        azimuthal angle between input points

    """

    delta_y = p1[1]-p0[1]
    delta_x = p1[0]-p0[0]
    azimuth = math.atan2(delta_y, delta_x)
    return azimuth


def sortPoints(anchor, points):
    """
    Sort all the points based on azimuthal (polar) angle they make with
    the anchor to start searching in counter-clockwise rotation

    Parameters
    ----------
    anchor : tuple
        Most soutwestern point of given point list
    points : list
        contains coordinate tuples

    Returns
    -------
    sorted_points : list
        coordinate tuples sorted by their azimuth angles between them
        in ascending order

    """

    azimuthlist = {
        calcAzimuth(anchor, point): (point[0], point[1])
        for point in points
                   }
    azimuths = list(azimuthlist.keys())
    azimuths.remove(0)

    sort_azi = sorted(azimuths, reverse=True)
    sorted_points = [azimuthlist[azimuth] for azimuth in sort_azi]
    return sorted_points


def ccw(p0, p1, p2):
    """
    Find traversing to a point from the previous two points makes a clockwise
    or a counter-clockwise direction

    Parameters
    ----------
    p0 : tuple
        first point
    p1 : tuple
        second point
    p3 : tuple
        third point

    Returns
    -------
    rotation : float
        Difference between the slopes to identify the third point is on left
        or right

    """

    rotation = (p1[0]- p0[0]) * (p2[1] - p1[1]) - (p2[0] - p1[0]) * (p1[1] - p0[1])
    return rotation


def convexHull(anchor, pointIDs, points, sortedPoints):
    """
    Find convex hull of given points

    Parameters
    ----------
    anchor : tuple
        Most soutwestern point of given point list
    pointIDs : list
        point IDs in order
    points : list
        contains coordinate tuples
    sortedPoints : list
        coordinate tuples sorted by their azimuth angles between them
        in ascending order

    Returns
    -------
    convex_IDs : list
        point IDs of hull points
    convex_hull : list
        coordinates of hull points

    """

    convex_hull = [anchor, sortedPoints[0]]
    for point in sortedPoints[1:]:
        while ccw(convex_hull[-2], convex_hull[-1], point) > 0:
            del convex_hull[-1]
        convex_hull.append(point)
    convex_IDs = [
        pointIDs[points.index(convex_point)]
        for convex_point in convex_hull
        ]
    hulls = {ID: hull for ID, hull in zip(convex_IDs, convex_hull)}
    return convex_IDs, convex_hull


def PlotConvexHull(pIDs, points, polygon):
    """
    Plot the points and convex hull

    Parameters
    ----------
    pIDs : list
        point IDs in order
    points : list
        contains coordinate tuples.
    polygon : list
        coordinates of hull points

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots(figsize=(12, 8))
    x, y = list(zip(*points))[0], list(zip(*points))[1]
    ax.scatter(y, x)
    polydraw = list(polygon)
    polydraw.append(polygon[0])
    c_x, c_y = list(zip(*polydraw))[0], list(zip(*polydraw))[1]
    ax.plot(c_y, c_x, color="red")
    plt.title("Convex Hull with Graham Scan Algorithm")
    plt.xlabel("Y")
    plt.ylabel("X")
    for i, txt in enumerate(pIDs):
        ax.annotate(txt, (y[i], x[i]),
                    verticalalignment="bottom",
                    horizontalalignment="right")
    return


def writeOutput(pointIDs, points, convexIDs):
    """
    Write the results in output report files in the directory

    Parameters
    ----------
    pointIDs : list
        point IDs in order
    points : list
        contains coordinate tuples.
    convexIDs : list
        point IDs of hull points

    Returns
    -------
    None.

    """

    outputfile = open(file + ".out", "w")
    outputfile.write("Point ID 	x [m] 	y [m]\n\
---------------------------\n")
    for ID in convexIDs:
        outputfile.write("{} {:.3f} {:.3f}\n".format(
            ID, points[pointIDs.index(ID)][0], points[pointIDs.index(ID)][1]
            ))
    print("Output report is generated as {}.out in the directory".format(file))
    outputfile.close()

    if file == "randompoints":
        outputfile = open("randompoints.xyz", "w")
        outputfile.write("{}\n".format(str(len(points))))
        for ID, point in zip(pointIDs, points):
            outputfile.write(
                "{} {:.3f} {:.3f}\n".format(ID, point[0], point[1])
                )
        print("Generated random points' coordinates are generated as \
{}.xyz in the directory".format(file))
        outputfile.close()
    return


def main():
    pointIDs, points = getInput()
    anchor = findAnchor(pointIDs, points)
    sorted_points = sortPoints(anchor, points)
    convex_IDs, convex_hull = convexHull(anchor, pointIDs,
                                         points, sorted_points)
    PlotConvexHull(pointIDs, points, convex_hull)
    writeOutput(pointIDs, points, convex_IDs)


if __name__ == "__main__":
    main()
