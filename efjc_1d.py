#!/usr/bin/python
import os
import math
import time
import multiprocessing
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from multiprocessing import Pool
from libmath import *

plt.rc('text',usetex=True)
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})

ResultsDir = 'results'

#Default Values

link_length=1
extendable_length=0.1


def ConstantParameter(n_):

    return 1/float(2*factorial(n_ - 1)) * pow(1/float(2*extendable_length) ,n_)
    

def Eta(n_, k1_, k2_, R_):

    return link_length*(n_ - 2*k1_) + (extendable_length/2)*(n_ - 2*k2_) - R_


def F(n_, k1_, k2_):

    return link_length*(n_ - 2*k1_) + (extendable_length/2)*(n_ - 2*k2_)


def TopHat(n_, k1_, k2_, R_):

    return heaviside_step(R_ - F(n_, k1_, n_)) * heaviside_step(F(n_, k1_, 0) - R_)


def Kernal(n_, k1_, k2_, R_):

    return pow(-1,k2_) * double_binomial(n_, k1_, k2_) * pow(Eta(n_, k1_, k2_, R_),n_-1) * sgn(Eta(n_, k1_, k2_, R_))


def KernalModified(n_, k1_, k2_, R_):

    return pow(-1,k2_) * double_binomial(n_, k1_, k2_) * pow(Eta(n_, k1_, k2_, R_),n_-1) * sgn(Eta(n_, k1_, k2_, R_)) * TopHat(n_, k1_, k2_, R_)


def PartitionFunction(n_, R_):

    SumPartitionFunction = 0
    for i in range(0,int(n_+1)):
        for j in range(0,int(n_+1)):
            SumPartitionFunction = SumPartitionFunction + Kernal(n_,i,j,R_)
    return ConstantParameter(n_) * SumPartitionFunction

def PartitionFunctionModified(n_, R_):

    SumPartitionFunction = 0
    for i in range(0,int(n_+1)):
        for j in range(0,int(n_+1)):
            SumPartitionFunction = SumPartitionFunction + KernalModified(n_,i,j,R_)
    return ConstantParameter(n_) * SumPartitionFunction


def Z0_data():

    nlistmax = 20

    Nlist = np.linspace(2,nlistmax,10)
    Z0list = []

    for n in Nlist:
        Z0list.append(PartitionFunction(n,0))
    
    plt.plot(Nlist, Z0list, marker='o', markeredgecolor='none', linestyle='--', color='#0e84b5', linewidth=2.0, zorder=100)
    plt.grid(color='0.85', linestyle='--')
    plt.xlabel(r'$N$',fontsize=12)
    plt.ylabel(r'$Z(R)/A^{N}$', fontsize=12)
    
    plt.xlim([0,nlistmax+2])
    plt.xticks(np.arange(min(Nlist), max(Nlist)+1, 4.0))
    plt.ylim([0,5.4])

    plt.title(r'Z as a function of $N$ for $R=0$', fontsize=14)
    plt.savefig(os.path.join(ResultsDir,'Z0_data.pdf'), bbox_inches='tight')
    plt.figure()


def ZR_data(n_):

    nrange = n_ + 1

    output_filename1 = 'n'+str(n_)+'_ZR_data.txt'
    out1 = open(os.path.join(ResultsDir,output_filename1),'w')

    output_filename2 = 'n'+str(n_)+'_ZR_mod_data.txt'
    out2 = open(os.path.join(ResultsDir,output_filename2),'w')

    nx = np.linspace(-nrange,nrange,(2*nrange+1)*100)
    zr = []
    zr_modified = []

    for i in range(len(nx)):
        zr.append(PartitionFunction(n_,nx[i]))
        zr_modified.append(PartitionFunctionModified(n_,nx[i]))

        out1.write(str(nx[int(i)]) + ',' + str(zr[int(i)]) + '\n')
        out2.write(str(nx[int(i)]) + ',' + str(zr_modified[int(i)]) + '\n')


    output_filename_plot_PDF = os.path.join(ResultsDir,'n'+str(n_)+'_ZR.pdf')
    PlotData(nx,zr,'R/a','Z(R)/A^{N}',output_filename_plot_PDF,'#05467C','#E6ECF2')

    output_filename_plot_PDF_mod = os.path.join(ResultsDir,'n'+str(n_)+'_ZR_mod.pdf')
    PlotData(nx,zr_modified,'R/a','Z(R)/A^{N}',output_filename_plot_PDF_mod,'#05467C','#E6ECF2')



def PlotData(xData_, yData_, xlabel_, ylabel_, outfilenamepdf_, plotcolor_, fillcolor_):

    plt.plot(xData_, yData_, color=plotcolor_, linewidth=2.0, zorder=100)
    plt.fill(xData_, yData_, color=fillcolor_) 
    plt.xlabel(r'$' + xlabel_ + '$',fontsize=12)
    plt.ylabel(r'$' + ylabel_ + '$', fontsize=12)
    plt.grid(color='0.85', linestyle='--')

    plt.ylim([0,max(yData_)+0.05*max(yData_)])

    plt.savefig(outfilenamepdf_, bbox_inches='tight')
    plt.figure()



def main():

    os.system('cls' if os.name=='nt' else 'clear')

    if os.path.exists(ResultsDir):
        for file in os.listdir(ResultsDir):
            file_path = os.path.join(ResultsDir, file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception, e:
                print e

    else:
        os.makedirs(ResultsDir)


    cpu_cores=multiprocessing.cpu_count()
    start = time.clock()

    pool = Pool(processes=cpu_cores)   

    N = [1,2,4,6,8,10]

    print "\nExtendable Freely Jointed Chain 1D"
    print "----------------------------------\n"
    print "Number of Links\t\t:", N
    print "Link Length\t\t:", link_length
    print "Extendable Link Length\t:", extendable_length
    print  
    
    Z0_data()

    pool.map(ZR_data, N)

    print 
    print 'Runtime (' + str(cpu_cores) + ' cpu cores) : ' + str(time.clock()-start)
    print

if __name__ == '__main__':
    main()

