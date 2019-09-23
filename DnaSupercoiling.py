import numpy as np
import math
from mayavi import mlab
import random
import scipy
from scipy import integrate

bases = [(0,0,1), (1,1,0), (1,0,0), (0,1,0)]
height = 100
twists = 2
writhes = 2

#DNA backbone points for values of t from 0 to height
def backbonePoints (t, shift = 0):

    if writhes > 0:
        (x0, y0, z0) = multipleWrithes (t * 2*math.pi/(height - 1))

        x = 1.185 * math.sin(t * 2 * math.pi / height * twists - shift) + x0
        y = 1.185 * math.cos(t * 2 * math.pi / height * twists - shift) + y0

        z = z0

    else:
        x = 1.185 * math.sin(t * 2 * math.pi / height * twists - shift)
        y = 1.185 * math.cos(t * 2 * math.pi / height * twists - shift)

        z = t

    print("t: " + str(t * 2*math.pi/(height - 1)))
    print("x: " + str(x))
    print("y: " + str(y))
    print("z: " + str(z))

    return (x, y, z)

#helper function to round decimals
def _roundUp(n, decimals=0):
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier

#creates the figure 8 like shape
#the function is called only when writhes > 1
def multipleWrithes (t, recursive = False):

    if (t > 2*math.pi):
        t = 2*math.pi

    smallCircles = writhes - 1
    circleRadius = height / ((4 + smallCircles) * math.pi)
    smallCircleRadius = circleRadius / 2
    numCircles = writhes + 1
    intervalSize = math.pi/numCircles

    if (t >= 0 and t <= intervalSize):

        if recursive:
            t1 = 2*math.pi - t

            x = circleRadius * math.sin(numCircles*t1)
            y = math.sin(t1)

        else:
            x = circleRadius * math.sin(numCircles*t)
            y = math.sin(t)

        z = circleRadius * math.cos(numCircles*t)

    elif (t > intervalSize and t <= math.pi - intervalSize):
        circleNumber = int (t/intervalSize)
        shift = 0
        if (circleNumber % 2 == 1):
            shift = math.pi

        if recursive:
            t1 = 2*math.pi - t
            x = smallCircleRadius * math.sin(numCircles*t1)
            y = math.sin(t1)
        else:
            x = smallCircleRadius * math.sin(numCircles*t)
            y = math.sin(t)

        z = smallCircleRadius * (math.cos(numCircles*t - shift)) - (3/2 + circleNumber - 1) * circleRadius

    elif (t > math.pi - intervalSize and t <= math.pi):
        shift = 0
        if (writhes % 2 == 1):
            shift = math.pi

        if recursive:
            t1 = 2*math.pi - t
            x = circleRadius * math.sin(numCircles*t1)
            y = math.sin(t1)
        else:
            x = circleRadius * math.sin(numCircles*t)
            y = math.sin(t)
        z = circleRadius * (math.cos(numCircles*t - shift)) - circleRadius * numCircles

    else:
        t = 2*math.pi - t
        (x, y, z) = multipleWrithes (_roundUp(t, 3), recursive = True)


    return (x, y, z)

#calculates arclength given t from 0 to 2pi
def _arcLength(t, recursive = False):

    smallCircles = writhes - 1
    circleRadius = height / ((4 + smallCircles) * math.pi)
    smallCircleRadius = circleRadius / 2
    numCircles = writhes + 1
    intervalSize = math.pi/numCircles
    if (t >= 0 and t <= intervalSize):

        if recursive:
            t1 = 2*math.pi - t
            arcLength1 = _arcLength(2 * math.pi - intervalSize)
            arcLength = scipy.integrate.quad (lambda x: math.sqrt((circleRadius * numCircles)**2 + math.cos(x) **2 ), 2 * math.pi - intervalSize, t1)[0]
            arcLength = arcLength + arcLength1
        else:
            arcLength = scipy.integrate.quad (lambda x: math.sqrt((circleRadius * numCircles)**2 + math.cos(x) **2 ), 0, t)[0]

    elif (t > intervalSize and t <= math.pi - intervalSize):

        if recursive:
            t1 = 2*math.pi - t
            arcLength1 = _arcLength(math.pi + intervalSize)
            arcLength = _smallArcLength(math.pi + intervalSize, t1)[0]
            arcLength = arcLength + arcLength1
        else:
            arcLength1 = _arcLength(intervalSize)
            arcLength = arcLength1 + _smallArcLength(intervalSize, t)[0]

    elif (t > math.pi - intervalSize and t <= math.pi):

        if recursive:
            t1 = 2*math.pi - t
            arcLength1 = _arcLength(math.pi)
            arcLength = scipy.integrate.quad (lambda x: math.sqrt((circleRadius * numCircles)**2 + math.cos(x) **2 ), math.pi, t1)[0]
            arcLength = arcLength + arcLength1
        else:

            arcLength1 = _arcLength(math.pi - intervalSize)
            arcLength = scipy.integrate.quad (lambda x: math.sqrt((circleRadius * numCircles)**2 + math.cos(x) **2 ), math.pi - intervalSize, t)[0]

            arcLength = arcLength + arcLength1


    else:
        t = 2*math.pi - t

        return(_arcLength (_roundUp(t, 3), recursive = True))

    return arcLength

#helper function to get arclength of the small circles created by multiple writhes
def _smallArcLength(a, b):

    smallCircles = writhes - 1
    numCircles = writhes + 1
    circleRadius = height / ((4 + smallCircles) * math.pi)
    smallCircleRadius = circleRadius / 2
    intervalSize = math.pi/numCircles


    circleNumber = int (b/intervalSize)

    shift = 0
    if (circleNumber % 2 == 1):
        shift = math.pi

    func_X = lambda x: smallCircleRadius * numCircles * math.cos(numCircles*x)
    func_Y = lambda x: math.cos(x)
    func_Z = lambda x: smallCircleRadius * numCircles * (- math.sin(numCircles * x - shift))

    arcLength = scipy.integrate.quad (lambda x: math.sqrt(func_X(x)**2 + func_Y(x)**2 + func_Z(x)**2 ), a, b)

    return arcLength

#plots DNA backbones
def plotBackbones():

    x = np.empty((height, 1))
    y = np.empty((height, 1))
    z = np.empty((height, 1))
    for t in range (height):

        (x1, y1, z1) = backbonePoints (t, 0)
        x[t][0] = x1
        y[t][0] = y1
        z[t][0] = z1

    mlab.plot3d(x, y, z, color = (0,0,0))

    x = np.empty((height, 1))
    y = np.empty((height, 1))
    z = np.empty((height, 1))
    for t in range (height):

        (x1, y1, z1) = backbonePoints (t, math.pi)
        x[t][0] = x1
        y[t][0] = y1
        z[t][0] = z1

    mlab.plot3d(x, y, z, color = (0,0,0))

#helper function to generate DNA color pairs
def _generateBasePairs(numRungs):

    basePairs = np.empty((numRungs, 2))
    for i in range (numRungs):
        num = random.randint(0,3)
        basePairs[i][0] = num
        cases = {
            0: 1,
            1: 0,
            2: 3,
            3: 2
        }
        basePairs[i][1] = cases.get(num)

    return basePairs

#plots DNA rungs
def plotRungs():
    totalArcLength = 0
    numRungs = int(height/0.332) - 2
    basePairs = _generateBasePairs(numRungs)
    counter = 0
    for t in range(int(2*math.pi * 100)):
        arcLength = _arcLength(t/100)
        if arcLength - totalArcLength > 0.332:

            x = np.empty((2, 1))
            y = np.empty((2, 1))
            z = np.empty((2, 1))
            (x1, y1, z1) = backbonePoints (t/(2*math.pi), math.pi)
            (x2, y2, z2) = backbonePoints (t/(2*math.pi))
            average_x = (x1 + x2)/2
            average_y = (y1 + y2)/2

            x[0][0] = x1
            x[1][0] = average_x
            y[0][0] = y1
            y[1][0] = average_y
            z[0][0] = z1
            z[1][0] = z1
            totalArcLength = arcLength

            baseColor = bases[int(basePairs[counter][0])]
            mlab.plot3d(x, y, z, color = baseColor)
            counter += 1

    totalArcLength = 0
    counter = 0
    for t in range(int(2*math.pi * 100)):
        arcLength = _arcLength(t/100)
        if arcLength - totalArcLength > 0.332:

            x = np.empty((2, 1))
            y = np.empty((2, 1))
            z = np.empty((2, 1))
            (x1, y1, z1) = backbonePoints (t/(2*math.pi), math.pi)
            (x2, y2, z2) = backbonePoints (t/(2*math.pi))
            average_x = (x1 + x2)/2
            average_y = (y1 + y2)/2

            x[0][0] = average_x
            x[1][0] = x2
            y[0][0] = average_y
            y[1][0] = y2
            z[0][0] = z1
            z[1][0] = z1
            totalArcLength = arcLength

            baseColor = bases[int(basePairs[counter][1])]
            mlab.plot3d(x, y, z, color = baseColor)
            counter += 1


plotBackbones()
plotRungs()
mlab.show()
