"""
Point in Polygon Ray Casting Method

@ Author: Res. Assist. Ömer Akın
@ Instituiton: Istanbul Technical University Geomatics Engineering Departmant
@ e-mail: akinom@itu.edu.tr

"""

import os, sys
import matplotlib.pyplot as plt


def getInput():
    """
    Generate points, polygon and checkpoints from the input project files

    Returns
    -------
    pointIDs : list
        IDs of points
    points : list
        point coordinates
    polyIDs : list
        vertex point IDs of polygon
    polygon : list
        vertex points of polygon
    checkIDs : list
        point IDs of points that will be checked
    checkpoints : list
        point coordinates that will be checked
    """

    global file
    while True:
        file = input("\nPress q & Enter to exit\n\
Please write the project name without extension: ")
        # Exit Conditions
        if file == "q":
            sys.exit()
        # Check file existence
        if (
                os.path.exists(file + ".set") and
                os.path.exists(file + ".pol") and
                os.path.exists(file + ".in")
                ):
            try:
                print("File exists in the directory")
                pointfile = open(file + ".set")
                polyfile = open(file + ".pol")
                checkfile = open(file + ".in")
                p_lines = pointfile.readlines()[2:]
                poly_lines = polyfile.read()
                check_lines = checkfile.read()
                pointfile.close()
                polyfile.close()
                checkfile.close()
                checkIDs = check_lines.split()
                polyIDs = poly_lines.split()
                pointIDs = [i.split()[0] for i in p_lines]
                x = [float(i.split()[1]) for i in p_lines]
                y = [float(i.split()[2]) for i in p_lines]
                points = [(i, j) for i, j in zip(x, y)]
                checkpoints = [points[pointIDs.index(pID)] for pID in checkIDs]
                polygon = [points[pointIDs.index(pID)] for pID in polyIDs]
                return pointIDs, points, polyIDs, polygon, checkIDs, checkpoints

            except ValueError:
                print("Some points are not found in given project.")
            except:
                print("File is not properly designed.")
            

        else:
            print("File does not exist in the directory")


def pointInPoly(checkpoints, poly):
    """
    Check points whether they are inside of the polygon or not

    Parameters
    ----------
    checkpoints : list
        point coordinates that will be checked
    poly : list
        vertex points of polygon

    Returns
    -------
    checkedpoints : list
        Boolean list with the same order of checkpoints that shows in/out
        situation (True for in, False for out)
    borderpoints : list
        Boolean list with the same order of checkpoints that shows whether the
        point on border or not (True for border, False for not)
    """

    polypoints = len(poly)
    j = polypoints - 1
    checkedpoints = []
    borderpoints = []
    for point in checkpoints:
        checklist = []
        borderlist = []
        for i in range(polypoints):
            check = False
            border = False
            if ((poly[i][0] > point[0]) != (poly[j][0] > point[0])) and \
                (point[1] < poly[i][1] + (poly[j][1] - poly[i][1]) *
                (point[0] - poly[i][0]) / (poly[j][0] - poly[i][0])):
                check = True
            elif ((poly[i][0] > point[0]) != (poly[j][0] > point[0])) and \
                (point[1] == poly[i][1] + (poly[j][1] - poly[i][1]) *
                (point[0] - poly[i][0]) / (poly[j][0] - poly[i][0])):
                border = True
            checklist.append(check)
            borderlist.append(border)
            countcheck = lambda x: True if x.count(True) % 2 != 0 else False
            bordercheck = lambda x: True if any(borderlist) else False
            j = i
        checkedpoints.append(countcheck(checklist))
        borderpoints.append(bordercheck(borderlist))

    return checkedpoints, borderpoints


def getIDs(checkedpoints, borderpoints, checkIDs, checkpoints):
    """
    Get Point IDs and coordinates by using boolean lists

    Parameters
    ----------
    checkedpoints : list
        Boolean list with the same order of checkpoints that shows in/out
        situation (True for in, False for out)
    borderpoints : list
        Boolean list with the same order of checkpoints that shows whether the
        point on border or not (True for border, False for not)
    checkIDs : list
        point IDs of points that will be checked
    checkpoints : list
        point coordinates that will be checked

    Returns
    -------
    inIDs : list
        point IDs of points that are inside of the polygon
    incoords : list
        coordinates of points that are inside of the polygon
    borderIDs : list
        point IDs of points that are on border of the polygon
    bordercoords : list
        coordinates of points that are on border of the polygon
    """

    inIDs = []
    for i in range(len(checkedpoints)):
        if checkedpoints[i]:
            inIDs.append(checkIDs[i])
    incoords = [checkpoints[checkIDs.index(ID)] for ID in inIDs]
    borderIDs = []
    for i in range(len(borderpoints)):
        if borderpoints[i]:
            borderIDs.append(checkIDs[i])
    bordercoords = [checkpoints[checkIDs.index(ID)] for ID in borderIDs]

    return inIDs, incoords, borderIDs, bordercoords


def pltPlot(pointIDs, points, polygon, inCoords, borderCoords, checkpoints):
    """
    Plot the points, polygon and checked points in corresponding colors

    Parameters
    ----------
    pointIDs : list
        IDs of points
    points : list
        point coordinates
    polygon : list
        vertex points of polygon

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
    ax.plot(c_y, c_x, color="black")
    check_x, check_y = list(zip(*checkpoints))[0], list(zip(*checkpoints))[1]
    ax.scatter(check_y, check_x, color="red", label="Checked Points")
    if any(inCoords):
        in_x, in_y = list(zip(*inCoords))[0], list(zip(*inCoords))[1]
        ax.scatter(in_y, in_x, color="green", label="In Points")
    if any(borderCoords):
        b_x, b_y = list(zip(*borderCoords))[0], list(zip(*borderCoords))[1]
        ax.scatter(b_y, b_x, color="yellow", label="Border Points")
    for ID, (x, y) in zip(pointIDs, points):
        ax.annotate(ID, (y, x),
                    verticalalignment="bottom",
                    horizontalalignment="right")

    plt.title("Point in Polygon (Ray Casting)")
    plt.xlabel("Y")
    plt.ylabel("X")
    plt.legend(loc="upper left")
    plt.show()

    return


def writeOutput(checkIDs, polyIDs, checkedpoints, borderIDs):
    """
    Write the result in an output report file in the directory

    Parameters
    ----------
    checkIDs : list
        point IDs of points that will be checked
    polyIDs : list
        vertex point IDs of polygon
    checkedpoints : list
        Boolean list with the same order of checkpoints that shows in/out
        situation (True for in, False for out)
    borderIDs : list
        point IDs of points that are on border of the polygon

    Returns
    -------
    None.
    """

    inOut = ["inside" if boo else "outside" for boo in checkedpoints]
    outputfile = open(file + ".out", "w")
    poly_string = ','.join(polyIDs)
    for ID, inout in zip(checkIDs, inOut):
        if ID in borderIDs:
            outputfile.write(
                "Point {} is on the border of polygon of {}\n".format(
                    ID, poly_string
                    )
                )
        else:
            outputfile.write(
                "Point {} is {} the polygon of {}\n".format(
                    ID, inout, poly_string
                    )
                )
    print("Output report is generated as {}.out in the directory".format(file))
    outputfile.close()

    return


def main():
    pointIDs, points, polyIDs, polygon, checkIDs, checkpoints = getInput()
    inPoints, border_points = pointInPoly(checkpoints, polygon)
    inIDs, incoords, borderIDs, border_coords = getIDs(
        inPoints, border_points, checkIDs, checkpoints
        )
    pltPlot(pointIDs, points, polygon, incoords, border_coords, checkpoints)
    writeOutput(checkIDs, polyIDs, inPoints, borderIDs)

main()
