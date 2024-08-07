#!/usr/env/bin python3

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Load the data file
data = []
with open('rxx.dat') as f:
    for line in f:
        splitLine = line.split()
        for i in range(len(splitLine)):
            if splitLine[i] != "|":
                splitLine[i] = float(splitLine[i])
        data.append(splitLine)

# Add the size of the SDP matrix to the data
for i in range(len(data)):
    sizeSDP = data[i][1]*2 + 1
    data[i].insert(2, sizeSDP)
    if (data[i][13] > 0):
        optimError = 100.0 * abs((data[i][7] - data[i][13]) / data[i][13])
        data[i].append(optimError)
    else:
        data[i].append(0)

# The names of each column
columnNames = ["Attempt number", "Number of inputs", "Matrix size", "Time to generate (ms)", 
"|", "Peak memory Optim (kB)", "Time taken Optim (ms)", "Value Optim", "Linear error Optim", "PSD error Optim", 
"|", "Peak memory MOSEK (kB)", "Time taken MOSEK (ms)", "Value MOSEK", "Linear error MOSEK", "PSD error MOSEK",
"|", "Peak memory SCS (kB)", "Time taken SCS (ms)", "Value SCS", "Linear error SCS", "PSD error SCS", 
"Percentage error in solution", "Standard deviation of percentage error", "Min and max percentage error"]

for i in range(len(columnNames)):
    print(i, columnNames[i])

# Average various quantities between the attempts
for i in range(0, len(data), 100):
    totalMemoryOptim = 0
    totalTimeOptim = 0
    totalMemoryMOSEK = 0
    totalTimeMOSEK = 0
    totalMemorySCS = 0
    totalTimeSCS = 0
    totalError = 0
    sd = 0
    maxError = 0
    minError = 100
    for j in range(0, 99):
        totalMemoryOptim += data[i+j][5]
        totalTimeOptim += data[i+j][6]
        totalMemoryMOSEK += data[i+j][11]
        totalTimeMOSEK += data[i+j][12]
        totalMemorySCS += data[i+j][17]
        totalTimeSCS += data[i+j][18]
        totalError += data[i+j][22]
        maxError = max(maxError, data[i+j][22])
        minError = min(minError, data[i+j][22])
    totalMemoryOptim /= 100
    totalTimeOptim /= 100
    totalMemoryMOSEK /= 100
    totalTimeMOSEK /= 100
    totalMemorySCS /= 100
    totalTimeSCS /= 100
    totalError /= 100
    for j in range(0, 99):
        sd += (data[i+j][22] - totalError)**2
    sd = (sd / 100)**0.5
    data[i][5] = totalMemoryOptim
    data[i][6] = totalTimeOptim
    data[i][11] = totalMemoryMOSEK
    data[i][12] = totalTimeMOSEK
    data[i][17] = totalMemorySCS
    data[i][18] = totalTimeSCS
    data[i][22] = totalError
    data[i].append(sd)
    data[i].append([minError, maxError])

# Increase the font size
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)

# Memory usage
xCol = 2
yCol1 = 5
yCol2 = 11
yCol3 = 17
xData1 = []
xData2 = []
yData1 = []
yData2 = []
yData3 = []
for i in range(0, len(data), 100):
    xData1.append(data[i][xCol])
    yData1.append(data[i][yCol1])
    if data[i][yCol2] > 0:
        xData2.append(data[i][xCol])
        yData2.append(data[i][yCol2])
    if data[i][yCol3] > 0:
        yData3.append(data[i][yCol3])
plt.clf()
plt.xlabel(columnNames[xCol])
yAxis = columnNames[yCol1]
if ("Optim" in yAxis):
    yAxis = yAxis.replace("Optim ", "")
plt.ylabel(yAxis)
plt.plot(xData1, yData1, 'bx', label="Ours")
plt.plot(xData2, yData2, 'ro', label="MOSEK")
plt.plot(xData2, yData3, 'g^', label="SCS")
plt.grid()
plt.yscale('log')
plt.legend()
plt.tight_layout()
# plt.show()
plt.savefig("memory.pdf")

# Time taken
xCol = 2
yCol1 = 6
yCol2 = 12
yCol3 = 18
xData1 = []
xData2 = []
yData1 = []
yData2 = []
yData3 = []
for i in range(0, len(data), 100):
    xData1.append(data[i][xCol])
    yData1.append(data[i][yCol1])
    if data[i][yCol2] > 0:
        yData2.append(data[i][yCol2])
        xData2.append(data[i][xCol])
    if data[i][yCol3] > 0:
        yData3.append(data[i][yCol3])
plt.clf()
plt.xlabel(columnNames[xCol])
yAxis = columnNames[yCol1]
if ("Optim" in yAxis):
    yAxis = yAxis.replace("Optim ", "")
plt.ylabel(yAxis)
plt.plot(xData1, yData1, 'bx', label="Ours")
plt.plot(xData2, yData2, 'ro', label="MOSEK")
plt.plot(xData2, yData3, 'g^', label="SCS")
plt.grid()
plt.yscale('log')
plt.legend()
plt.tight_layout()
# plt.show()
plt.savefig("time.pdf")

# Optim error
xCol = 2
yCol1 = 22
xData1 = []
xData2 = []
yData1 = []
yError1 = []
for i in range(0, len(data), 100):
    if data[i][22] > 0:
        xData1.append(data[i][xCol])
        yData1.append(data[i][yCol1])
        yError1.append(data[i][23])
yError1 = np.array(yError1).T
plt.clf()
plt.xlabel(columnNames[xCol])
yAxis = columnNames[yCol1]
if ("Optim" in yAxis):
    yAxis = yAxis.replace("Optim ", "")
plt.ylabel(yAxis)
plt.errorbar(xData1, yData1, yerr=yError1, marker='x', fmt=' ', capsize=5)
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("error.pdf")
plt.show()




