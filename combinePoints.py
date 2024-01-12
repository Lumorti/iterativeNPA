import os
from scipy.spatial import ConvexHull
import numpy as np

# Write points to a file
def writePoints(points, filename):
    print("Writing ", len(points), " points to file " + filename + "...")
    with open(filename, "w") as f:
        for p in points:
            for i in range(len(p)):
                f.write(str(p[i]))
                if i < len(p)-1:
                    f.write(",")
            f.write("\n")

# Load points from a file
def readPoints(filename):
    print("Loading points from file " + filename + "...")
    with open(filename, "r") as f:
        points = []
        for line in f:
            x = [float(i) for i in line.split(",")]
            points.append(x)
        print("Loaded ", len(points), " points")
        return points

# Here I want to combine every csv in this folder, 
# adding a column based on the filename numbers
# Valid filenames are points1.2.csv, points-1.3.csv, etc.
combinedPoints = []
for filename in os.listdir("data"):
    if filename.startswith("points") and filename.endswith(".csv"):
        try:
            num = float(filename[6:-4])
            print(num)
            newPoints = readPoints("data/" + filename)
            for p in newPoints:
                p.append(num)
            combinedPoints.extend(newPoints)
        except:
            pass
combinedPoints = np.array(combinedPoints)
hull = ConvexHull(combinedPoints)
combinedPoints = combinedPoints[hull.vertices,:]
writePoints(combinedPoints, "data/pointsAll.csv")



