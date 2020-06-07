"""
Open, Closed and Linked Traverse Calculator

@ Author: Res. Assist. Ömer Akın
@ Instituiton: Istanbul Technical University Geomatics Engineering Departmant
@ e-mail: akinom@itu.edu.tr

"""


import os, sys, math


def showHelp():
    print("\n\n###################----HELP MENU----###################\n\
User should give the project name that holds the following conditions.\n\
Please note that, there must be three seperate files having same name\
to run the script. \n\
A file with .xyz extension that holds the reference points' coordinates. \n\
A file with .ang extension that holds the traverse angle measurements. \n\
A file with .leg extension that holds the distance measurements.\n\n")


def getInputsfromFile():
    """
    Get Reference points coordinates, traverse angles, distances between
    traverse points and point IDs

    Returns
    -------
    referencepoints : list
        known points' coordinates
    traverseangles : list
        traverse angles
    distances : list
        distances between points
    pointIDs : list
        known points' IDs
    traverseIDs : list
        traverse points' IDs
    to_pID[-1] :
        DESCRIPTION.

    """

    global projectfile
    while True:
        projectfile = input("Please enter the name of files: ")
        try:
            assert (os.path.exists(projectfile + ".xyz") and
                    os.path.exists(projectfile + ".ang") and
                    os.path.exists(projectfile + ".leg"))

            # Get Reference Points
            pointfile = open(projectfile + ".xyz")
            p_lines = pointfile.readlines()[2:]
            pointIDs = [i.split()[0] for i in p_lines]
            xlist = [i.split()[1] for i in p_lines]
            ylist = [i.split()[2] for i in p_lines]
            pointfile.close()
            referencepoints = [
                (float(i), float(j)) for i, j in zip(xlist, ylist)
                ]

            # Get Traverse Angles
            anglefile = open(projectfile + ".ang")
            a_lines = anglefile.readlines()[2:]
            traverseIDs = [i.split()[0] for i in a_lines]
            traverseangles = [float(i.split()[1]) for i in a_lines]
            anglefile.close()

            # Get Distances
            distancefile = open(projectfile + ".leg")
            d_lines = distancefile.readlines()[2:]
            #from_pID = [i.split()[0] for i in d_lines]
            to_pID = [i.split()[1] for i in d_lines]
            distances = [float(i.split()[2]) for i in d_lines]
            distancefile.close()
            return (referencepoints, traverseangles, distances,
                    pointIDs, traverseIDs, to_pID[-1])

        except AssertionError:
            help_input = input("Files are not found or \
file names are not properly designed. \nPlease press 'h' for instructions, \
q to quit, any key to try again: ")
            if help_input == "h":
                showHelp()
            elif help_input == "q":
                sys.exit()


def findAzimuth(x1, y1, x2, y2):
    """
    Calculate the azimuth angle between two points

    Parameters
    ----------
    x1 : float
        x coordinate of the first point
    y1 : float
        y coordinate of the first point
    x2 : float
        x coordinate of the second point
    y2 : float
        y coordinate of the second point

    Returns
    -------
    azimuth : float
        azimuth angle between given two points

    """

    delta_y = y2 - y1
    delta_x = x2 - x1
    try:
        azimuth = math.atan(abs(delta_y)/abs(delta_x))*(200/math.pi)
        if delta_y > 0 and delta_x < 0:
            azimuth = 200 - azimuth
        elif delta_y < 0 and delta_x < 0:
            azimuth = 200 + azimuth
        elif delta_y < 0 and delta_x > 0:
            azimuth = 400 - azimuth
    except ZeroDivisionError:
        if delta_y > 0:
            azimuth = 100
        elif delta_y < 0:
            azimuth = 300
    azimuth = round(azimuth, 4)
    return azimuth


def calcAzimuths(initial_azimuth, traverseangles):
    """
    Calculate azimuth angles between the points within the given traverse

    Parameters
    ----------
    initial_azimuth : float
        azimuth angle between first point and the second point
    traverseangles : list
        traverse angles

    Returns
    -------
    azimuthList : list
        azimuths between all possible points

    """

    t = initial_azimuth
    azimuthList = [t]
    for angle in traverseangles:
        tx = t + angle
        if tx > 200 and tx < 600:
            t = tx - 200
        elif tx > 600:
            t = tx - 600
        elif tx < 200:
            t = tx + 200
        azimuthList.append(round(t, 4))
    return azimuthList


def openTraverse(referencepoints, traverseangles, distances):
    """
    Perform open traverse calculations

    Parameters
    ----------
    referencepoints : list
        known points' coordinates
    traverseangles : list
        traverse angles
    distances : list
        distances between points

    Returns
    -------
    traversecoords : list
        coordinates of traverse points
    azimuths : list
        azimuths between all possible points
    deltaX : list
        differences between X coordinates
    deltaY : list
        differences between X coordinates

    """
    # Open Traverse Calculations
    t_initial = findAzimuth(referencepoints[0][0], referencepoints[0][1],
                            referencepoints[1][0], referencepoints[1][1])

    azimuths = calcAzimuths(t_initial, traverseangles)

    deltaX = []
    deltaY = []
    for i in range(len(distances)):
        deltay = distances[i] * math.sin(
            azimuths[i + 1]*(math.pi / 200)
            )
        deltax = distances[i] * math.cos(
            azimuths[i + 1]*(math.pi / 200)
            )
        deltaY.append(round(deltay, 3))
        deltaX.append(round(deltax, 3))

    x, y = referencepoints[1][0], referencepoints[1][1]
    traversecoords = []
    for i, j in zip(deltaX, deltaY):
        x += i
        y += j
        traversecoords.append((round(x, 3), round(y, 3)))

    return traversecoords, azimuths, deltaX, deltaY


def closedTraverse(referencepoints, traverseangles, distances):
    """
    Perform closed traverse calculations

    Parameters
    ----------
    referencepoints : list
        known points' coordinates
    traverseangles : list
        traverse angles
    distances : list
        distances between points

    Returns
    -------
    traversecoords : list
        coordinates of traverse points
    correctedazimuths : list
        corrected azimuths between all possible points
    deltaXcorrected : list
        corrected differences between X coordinates
    deltaYcorrected : list
        corrected differences between X coordinates
    angle_error : float
        angular mis-closure error
    angle_tolerance : float
        angular mis-closure error tolerance
    fs : float
        linear mis-closure error
    Fs : float
        linear mis-closure error tolerance

    """

    t_initial = findAzimuth(referencepoints[0][0], referencepoints[0][1],
                            referencepoints[1][0], referencepoints[1][1])

    # Angular mis-closure error check
    angle_tolerance = (1.5 * math.sqrt(len(traverseangles))) / 100
    if (sum(traverseangles) % 200) > 190:
        angle_error = (sum(traverseangles) % 200) - 200
    else:
        angle_error = sum(traverseangles) % 200
    if abs(angle_error) > angle_tolerance:
        print("fb: {}\nFb: {}\nAngular mis-closure error is bigger than \
the tolerance value".format(round(angle_error, 4), round(angle_tolerance, 4)))
        sys.exit()
    else:
        a_correction = -angle_error / len(traverseangles)
        correctedangles = [
            angle + round(a_correction, 4) for angle in traverseangles
            ]
        correctedazimuths = calcAzimuths(t_initial, correctedangles)

    # Linear mis-closure error check
    deltaXlist = []
    deltaYlist = []
    for i in range(len(distances)):
        deltay = distances[i] * math.sin(
            correctedazimuths[i + 1]*(math.pi / 200)
            )
        deltax = distances[i] * math.cos(
            correctedazimuths[i + 1]*(math.pi / 200)
            )
        deltaYlist.append(round(deltay, 3))
        deltaXlist.append(round(deltax, 3))

    totaldeltaY, totaldeltaX = sum(deltaYlist), sum(deltaXlist)
    fx = round(-totaldeltaX, 3)
    fy = round(-totaldeltaY, 3)
    fs = math.sqrt(fx**2 + fy**2)
    Fs = 0.004*math.sqrt(sum(distances)) + 0.0003 * sum(distances) + 0.02
    if abs(fs) > Fs:
        print("fs: {} \nFs: {}\nLinear mis-closure error is bigger than the \
tolerance value.".format(round(fs, 3), round(Fs, 3)))
        sys.exit()
    else:
        deltaXcorrected = []
        deltaYcorrected = []
        for i in range(len(distances)):
            y_correction = ((fy / sum(distances)) * distances[i])
            deltaYcorrected.append(deltaYlist[i] + round(y_correction, 3))
            x_correction = ((fx / sum(distances)) * distances[i])
            deltaXcorrected.append(deltaXlist[i] + round(x_correction, 3))

    x, y = referencepoints[1][0], referencepoints[1][1]
    traversecoords = []
    for i, j in zip(deltaXcorrected, deltaYcorrected):
        x += i
        y += j
        traversecoords.append((round(x, 3), round(y, 3)))

    print("Calculations are done successfully")

    return (traversecoords, correctedazimuths, deltaXcorrected,
            deltaYcorrected, angle_error, angle_tolerance, fs, Fs)


def linkedTraverse(referencepoints, traverseangles, distances):
    """
    Perform linked traverse calculations

    Parameters
    ----------
    referencepoints : list
        known points' coordinates
    traverseangles : list
        traverse angles
    distances : list
        distances between points

    Returns
    -------
    traversecoords : list
        coordinates of traverse points
    correctedazimuths : list
        corrected azimuths between all possible points
    deltaXcorrected : list
        corrected differences between X coordinates
    deltaYcorrected : list
        corrected differences between X coordinates
    angle_error : float
        angular mis-closure error
    angle_tolerance : float
        angular mis-closure error tolerance
    fl : float
        linear mis-closure error
    fq : float
        linear mis-closure error
    Fl : float
        linear mis-closure error tolerance
    Fq : float
        linear mis-closure error tolerance

    """

    t_initial = findAzimuth(referencepoints[0][0], referencepoints[0][1],
                            referencepoints[1][0], referencepoints[1][1])
    t_final = findAzimuth(referencepoints[-2][0], referencepoints[-2][1],
                          referencepoints[-1][0], referencepoints[-1][1])

    first_azimuths = calcAzimuths(t_initial, traverseangles)

    # Angular mis-closure error check
    angle_tolerance = (1.5 * math.sqrt(len(traverseangles))) / 100
    angle_error = first_azimuths[-1] - t_final
    if abs(angle_error) > angle_tolerance:
        print("fb: {}\nFb: {}\nAngular mis-closure error is bigger than \
the tolerance value".format(round(angle_error, 4), round(angle_tolerance, 4)))
        sys.exit()
    else:
        a_correction = -angle_error / len(traverseangles)
        correctedangles = [
            angle + round(a_correction, 4) for angle in traverseangles
            ]
        correctedazimuths = calcAzimuths(t_initial, correctedangles)

    # Linear mis-closure error check
    deltaXlist = []
    deltaYlist = []
    for i in range(len(distances)):
        deltay = distances[i] * math.sin(
            correctedazimuths[i + 1]*(math.pi / 200)
            )
        deltax = distances[i] * math.cos(
            correctedazimuths[i + 1]*(math.pi / 200)
            )
        deltaYlist.append(round(deltay, 3))
        deltaXlist.append(round(deltax, 3))

    totaldeltaY, totaldeltaX = sum(deltaYlist), sum(deltaXlist)
    fx = round((referencepoints[-2][0] - referencepoints[1][0]
               - totaldeltaX), 3)
    fy = round((referencepoints[-2][1] - referencepoints[1][1]
               - totaldeltaY), 3)
    S = math.sqrt(totaldeltaY**2 + totaldeltaX**2)
    fq = (1 / S) * (fy * totaldeltaX - fx * totaldeltaY)
    fl = (1 / S) * (fy * totaldeltaY + fx * totaldeltaX)
    Fq = 0.05 + 0.15 * math.sqrt(sum(distances)*0.001)
    Fl = 0.05 + 0.04 * math.sqrt(len(traverseangles) - 1)
    if abs(fq) > Fq or abs(fl) > Fl:
        print("fq: {}, fl: {} \nFq: {}, Fl: {} \nLinear mis-closure error \
is bigger than the tolerance value.".format(round(fq, 3), round(fl, 3),
                                            round(Fq, 3), round(Fl, 3)))
        sys.exit()
    else:
        deltaXcorrected = []
        deltaYcorrected = []
        for i in range(len(distances)):
            y_correction = ((fy / sum(distances)) * distances[i])
            deltaYcorrected.append(deltaYlist[i] + round(y_correction, 3))
            x_correction = ((fx / sum(distances)) * distances[i])
            deltaXcorrected.append(deltaXlist[i] + round(x_correction, 3))

    x, y = referencepoints[1][0], referencepoints[1][1]
    traversecoords = []
    for i, j in zip(deltaXcorrected, deltaYcorrected):
        x += i
        y += j
        traversecoords.append((round(x, 3), round(y, 3)))

    print("Calculations are done successfully")

    return (traversecoords, correctedazimuths, deltaXcorrected,
            deltaYcorrected, angle_error, angle_tolerance, fl, fq, Fl, Fq)


def writeOutput(coordinates, lp, tIDs, traversetype="Linked"):
    """
    Write the traverse points' coordinates to a report file in the directory

    Parameters
    ----------
    coordinates : TYPE
        DESCRIPTION.
    lp : TYPE
        DESCRIPTION.
    tIDs : TYPE
        DESCRIPTION.
    traversetype : TYPE, optional
        DESCRIPTION. The default is "Linked".

    Returns
    -------
    None.

    """

    coordfile = open("traverseCoords.xyz", "w")
    coordfile.write("Point ID   X(m) Y(m) Z(m)\n\
--------   ---  ---  ---\n")
    if traversetype == "Open":
        for i in range(len(coordinates)):
            if i == len(coordinates) - 1:
                coordfile.write(
                    "{} {:.3f} {:.3f} 0.000\n".format(lp,
                                                      coordinates[i][0],
                                                      coordinates[i][1]))
            else:
                coordfile.write(
                    "{} {:.3f} {:.3f} 0.000\n".format(tIDs[i+1],
                                                      coordinates[i][0],
                                                      coordinates[i][1]))
    else:
        for i in range(len(coordinates) - 1):
            coordfile.write(
                "{} {:.3f} {:.3f} 0.000\n".format(tIDs[i+1],
                                                  coordinates[i][0],
                                                  coordinates[i][1]))
    coordfile.close()
    print("Coordinates of traverse points are written on 'traverseCoords.xyz'")
    return


def writeReport(pIDs, tIDs, b, t, d, x, y, coords, points, lp, ttype):
    """
    Write output report table

    """
    try:
        outfl = open(projectfile + ".txt", "w")
        outfl.write("|{}|{}|{}|{}|{}|{}|{}|{}|\n\
".format("-"*10, "-"*10, "-"*10, "-"*10, "-"*10, "-"*10, "-"*15, "-"*15))
        outfl.write("|Station ID|  B(Grad) |  t(Grad) |   S(m)   |   DY(m)  \
|   DX(m)  |      Y(m)     |      X(m)     |\n")
        outfl.write("|{}|{}|{}|{}|{}|{}|{}|{}|\n\
".format("-"*10, "-"*10, "-"*10, "-"*10, "-"*10, "-"*10, "-"*15, "-"*15))
        outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}|  \
{:<13}|\n".format(pIDs[0], " ", " ", " ",
                  " ", " ", points[0][1], points[0][0]))
        outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}|  \
{:<13}|\n".format(" ", " ", t[0], " ", " ", " ", " ", " "))
        outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}|  \
{:<13}|\n".format(pIDs[1], b[0], " ", " ",
                  " ", " ", points[1][1], points[1][0]))
        outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}|  \
{:<13}|\n".format(" ", " ", t[1], d[0],
                  round(y[0], 3), round(x[0], 3), " ", " "))
        outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}|  \
{:<13}|\n".format(
                  tIDs[1], b[1], " ", " ", " ",
                  " ", coords[0][1], coords[0][0])
            )
        outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}|  \
{:<13}|\n".format(
                  " ", " ", t[2], d[1], round(y[1], 3),
                  round(x[1], 3), " ", " ")
            )
        outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}|  \
{:<13}|\n".format(tIDs[2], b[2], " ", " ",
                  " ", " ", coords[1][1], coords[1][0]))
        outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}|  \
{:<13}|\n".format(" ", " ", t[3], d[2],
                  round(y[2], 3), round(x[2], 3), " ", " "))
        if ttype == "Linked":
            outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}\
|  {:<13}|\n".format(
                tIDs[3], b[3], " ", " ", " ", " ", coords[2][1], coords[2][0])
                )
            outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}\
|  {:<13}|\n".format(
                " ", " ", t[4], d[3], round(y[3], 3), round(x[3], 3), " ", " ")
                )
            if len(tIDs) > 5:
                outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  \
{:<13}|  {:<13}|\n".format(tIDs[4], b[4], " ", " ", " ",
                           " ", coords[3][1], coords[3][0]))
                outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  \
{:<13}|  {:<13}|\n".format(" ", " ", t[5], d[4],
                           round(y[4], 3), round(x[4], 3), " ", " "))
            outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}\
|  {:<13}|\n".format(pIDs[-2], b[-1], " ", " ", " ", " ",
                     points[-2][1], points[-2][0]))
            outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}\
|  {:<13}|\n".format(" ", " ", t[-1], " ", " ", " ", " ", " "))
            outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}\
|  {:<13}|\n".format(pIDs[-1], " ", " ", " ", " ", " ",
                     points[-1][1], points[-1][0]))
        elif ttype == "Closed":
            outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}\
|  {:<13}|\n".format(tIDs[3], b[3], " ", " ", " ", " ",
                     coords[2][1], coords[2][0]))
            outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}\
|  {:<13}|\n".format(" ", " ", t[4], d[3],
                     round(y[3], 3), round(x[3], 3), " ", " "))
            if len(tIDs) > 5:
                outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  \
{:<13}|  {:<13}|\n".format(tIDs[4], b[4], " ", " ", " ", " ",
                           coords[3][1], coords[3][0]))
                outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  \
{:<13}|  {:<13}|\n".format(" ", " ", t[5], d[4],
                           round(y[4], 3), round(x[4], 3), " ", " "))
            outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}\
|  {:<13}|\n".format(pIDs[-1], b[-1], " ", " ", " ", " ",
                     points[-1][1], points[-1][0]))
        elif ttype == "Open":
            outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}\
|  {:<13}|\n".format(tIDs[3], b[3], " ", " ", " ", " ",
                     coords[2][1], coords[2][0]))
            outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}\
|  {:<13}|\n".format(" ", " ", t[4], d[3],
                     round(y[3], 3), round(x[3], 3), " ", " "))
            if len(tIDs) > 4:
                outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  \
{:<13}|  {:<13}|\n".format(tIDs[4], b[4], " ", " ", " ", " ",
                           coords[3][1], coords[3][0]))
                outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  \
{:<13}|  {:<13}|\n".format(" ", " ", t[5], d[4],
                           round(y[4], 3), round(x[4], 3), " ", " "))
            outfl.write("|{:<10}| {:<9}| {:<9}| {:<9}| {:<9}| {:<9}|  {:<13}\
|  {:<13}|\n".format(lp, b[-1], " ", " ", " ", " ",
                     coords[-1][1], coords[-1][0]))
        outfl.write("|{}|{}|{}|{}|{}|{}|{}|{}|\n".format(
            "-"*10, "-"*10, "-"*10, "-"*10, "-"*10, "-"*10, "-"*15, "-"*15)
            )
        print("Output report is generated as {}.txt in the directory".format(
            projectfile)
            )
        outfl.close()
    except:
        print("Report could not be generated.")
    return


def main():
    while True:
        inputtype = input("Welcome to traverse calculator. \n\
To exit, please press 'q' \nTo get help, please press 'h' \n\n\
Please select your traverse type: \n\
For Open, please press 1 \n\
For Closed, please press 2 \n\
For Linked, please press 3: \n")
        if inputtype == "q":
            sys.exit()
        elif inputtype == "h":
            showHelp()
        elif inputtype == "1":
            (points, angles, distances, pIDs, tIDs, lp) = getInputsfromFile()
            coords, t, deltax, deltay = openTraverse(
                points, angles, distances
                )
            writeOutput(coords, lp, tIDs, "Open")
            writeReport(pIDs, tIDs, angles, t, distances,
                        deltax, deltay, coords, points, lp, "Open")
            break

        elif inputtype == "2":
            (points, angles, distances, pIDs, tIDs, lp) = getInputsfromFile()
            (coords, t, deltax, deltay, fb, Fb, fs, Fs) = closedTraverse(
                points, angles, distances
                )
            writeOutput(coords, lp, tIDs)
            writeReport(pIDs, tIDs, angles, t, distances,
                        deltax, deltay, coords, points, lp, "Closed")
            print("fb = {}g\t\t\t\t\tfs = {}m\nFb = \
{}g\t\t\t\tFs = {}m\n\nfb < Fb\t\t\t\t\t\tfs < Fs".format(
                round(fb, 4), round(fs, 3),
                round(Fb, 4), round(Fs, 3)))
            break

        elif inputtype == "3":
            (points, angles, distances, pIDs, tIDs, lp) = getInputsfromFile()
            (coords, t, deltax, deltay, fb,
             Fb, fl, fq, Fl, Fq) = linkedTraverse(
                points, angles, distances
                )
            writeOutput(coords, lp, tIDs)
            writeReport(pIDs, tIDs, angles, t, distances,
                        deltax, deltay, coords, points, lp, "Linked")
            print("fb = {}g\t\t\t\t\tfl = {}m\tfq = {}m\nFb = \
{}g\t\t\t\t\tFl = {}m\tFq = {}m\n\nfb < Fb\t\t\t\t\t\t\t\
fl < Fl & fq < Fq".format(
                        round(fb, 4), round(fl, 3), round(fq, 3),
                        round(Fb, 4), round(Fl, 3), round(Fq, 3)))
            break


if __name__ == "__main__":
    main()
