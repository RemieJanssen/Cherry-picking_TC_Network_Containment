from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import csv
import seaborn as sns
import pandas

def LinearFit(csvFileLocation):
    allData = list()
    allTimes = list()
    allSum = list()
    trueData = list()
    trueTimes = list()
    falseData = list()
    falseTimes = list()

    with open(csvFileLocation, 'r') as csvFile:
        reader = csv.reader(csvFile, delimiter=';')
        reps = -1
        for i,line in enumerate(reader):
            if i!=0:
                if int(line[3])>reps:
                    reps = int(line[3])
                else:
                    break
        reps = reps+1
        reader = csv.reader(csvFile, delimiter=';')
        for i,line in enumerate(reader):
            if i>0:
                allData.append((float(line[0]),float(line[1]),float(line[2])))
                allSum.append((float(line[0])+float(line[1])+float(line[2])))
                allTimes.append(float(line[5]))
                if line[4] == " True ":
                    trueData.append((float(line[0]),float(line[1]),float(line[2])))
                    trueTimes.append(float(line[5]))
                    if int(line[3]) <reps/2:
                        print("oh no, this should not be a subnetwork!", i) 
                else:
                    falseData.append((float(line[0]),float(line[1]),float(line[2])))
                    falseTimes.append(float(line[5]))
                    if int(line[3])>=reps/2:
                        print("oh no, this should be a subnetwork!", i) 

    allModel = LinearRegression(fit_intercept = False).fit(allData, allTimes)
    all_r_sq = allModel.score(allData, allTimes)
    print('all data')
    print('coefficient of determination:', all_r_sq)
    print('intercept:', allModel.intercept_)
    print('slope:', allModel.coef_)


    
    trueModel = LinearRegression(fit_intercept = False).fit(trueData, trueTimes)
    true_r_sq = trueModel.score(trueData, trueTimes)
    print("a subnetwork")
    print('coefficient of determination:', true_r_sq)
    print('intercept:', trueModel.intercept_)
    print('slope:', trueModel.coef_)

    falseModel = LinearRegression(fit_intercept = False).fit(falseData, falseTimes)
    false_r_sq = falseModel.score(falseData, falseTimes)
    print("not a subnetwork")
    print('coefficient of determination:', false_r_sq)
    print('intercept:', falseModel.intercept_)
    print('slope:', falseModel.coef_)



def OpenDataSet(csvFileLocation):
    with open(csvFileLocation) as inputfile:
        datafile = pandas.read_csv(inputfile, sep=';', header=0)
        return datafile
    return Null

def FigureRunningTimeLeaves(dataSet,r1,r2):
    fixedRetics1 = r1
    fixedRetics2 = r2
    dataSetFixed = dataSet[(dataSet["reticulations"]==fixedRetics1) & (dataSet["reticulations_subnetwork"]==fixedRetics2)]
    sns_plot = sns.lmplot(x="leaves", y="running_time", hue="subnetwork", data=dataSetFixed, markers=["o", "x"])
    axes = sns_plot.axes.flatten()
    axes[0].set_title("r="+str(fixedRetics1)+", r'="+str(fixedRetics2))
    sns_plot.set_axis_labels("leaves","time (s)")
    sns_plot._legend.set_title("Subnetwork")
    sns_plot._legend.texts[0].set_text("No")
    sns_plot._legend.texts[1].set_text("Yes")
    sns_plot.savefig("LeavesRunningTime_r1_"+str(fixedRetics1)+"_ r2_"+str(fixedRetics2)+".pdf")        



def FigureNegligableReticulationSubnetwork(dataSet,n,r1):
    leaves = n
    fixedRetics1 = r1
    dataSetFixed = dataSet[(dataSet["reticulations"]==fixedRetics1) & (dataSet["leaves"]==n)]
    sns_plot = sns.lmplot(x="reticulations_subnetwork", y="running_time", hue="subnetwork", data=dataSetFixed, markers=["o", "x"])
    axes = sns_plot.axes.flatten()
    axes[0].set_title("n="+str(leaves)+", r="+str(fixedRetics1))
    sns_plot.set_axis_labels("r'","time (s)")
    sns_plot._legend.set_title("Subnetwork")
    sns_plot._legend.texts[0].set_text("No")
    sns_plot._legend.texts[1].set_text("Yes")
    sns_plot.savefig("ReticulationSubnetworkRunningTime_n_"+str(leaves)+"_ r1_"+str(fixedRetics1)+".pdf")        



csvFileLocation = "./data.txt"
LinearFit(csvFileLocation)
dataSet = OpenDataSet(csvFileLocation)
if not dataSet.empty:
    FigureRunningTimeLeaves(dataSet,1000,1000)
    FigureNegligableReticulationSubnetwork(dataSet,1000,1000)













