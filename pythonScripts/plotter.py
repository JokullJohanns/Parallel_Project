import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate



def plot(data, title, legend_pos):
    plt.cla()
    fig = plt.figure(1)
    ax = plt.gca()
    ax.yaxis.grid(True) 
    
    #plt.title("k-means on Tegner")
    fig.suptitle(title, fontsize=20)
    plt.xlabel('#Processors')
    plt.ylabel('Speedup')

    plt.plot(labels, data[0],label='32')
    plt.plot(labels, data[1],label='16')
    plt.plot(labels, data[2],label='8')
    plt.plot(labels, data[3],label='4')
    plt.plot(labels, data[4],label='2')
    #plt.plot(t_X, data[2],'b',label='Predicted function')
    plt.xlim([3,24])
    plt.legend(loc=legend_pos, title="#dimensions")
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize='medium')
    fig.patch.set_facecolor('white')
    plt.show()


timed = [np.array([ 66.84,  35.47,  17.63,   9.56,   6.06]), 
        np.array([ 34.2 ,  18.23,   9.21,   4.82,   3.32]), 
        np.array([ 17.75,   9.5 ,   5.03,   2.58,   1.72]), 
        np.array([ 8.21,  4.36,  2.59,  1.43,  0.87])]
        #np.array([ 163.34,   88.97,   50.56,   27.68,   17.23])]

speedup =      np.array([np.array([2.44374626, 2.508316888, 2.867838911, 2.89539749, 2.843234323]),
                np.array([4.776023392, 4.880416895, 5.489685125, 5.742738589, 5.189759036]),
                np.array([9.202253521, 9.365263158, 10.05168986, 10.72868217, 10.01744186]),
                np.array([19.8952497, 20.4059633, 19.52123552, 19.35664336, 19.8045977])])

labels = [3, 6, 12, 24]
xAxis = [2,4,8,16,32]
timed = np.array(timed)


#plot(timed.T, "k-means on Tegner", "upper right")
plot(speedup.T, "Speedups","upper left")
#print tabulate(speedup,tablefmt="latex")