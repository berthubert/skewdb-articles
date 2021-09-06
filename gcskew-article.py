#!/usr/bin/env python
# coding: utf-8

# In[104]:



from IPython.display import set_matplotlib_formats
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

import matplotlib
import math
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [9.5, 7]
plt.rcParams['font.size'] = '13'

#plt.rcParams['figure.figsize'] = [15, 6]
import datetime
import pandas
from random import shuffle
from scipy.stats import poisson
import statsmodels.api as sm
import statsmodels.formula.api as smf


# In[105]:


# The goal of this notebook is to recreate all the graphs from the skewdb article


# In[106]:



prefix="/home/ahu/git/antonie/genomes/ncbi-genomes-2021-08-19/"
#data=pandas.read_csv(prefix+"/skplot.csv")
#data.describe()


# In[107]:


results = pandas.read_csv(prefix+"/results.csv")
results.sort_values(["rmsGC","rmsTA"], inplace=True)
results["gccontent"]=results.gccount/results.siz
results.describe()
gennames=pandas.read_csv(prefix+"/genomes.csv", sep=';')
codongc=pandas.read_csv(prefix+"/codongc.csv", sep=';')


m=gennames.merge(results, on="name").merge(codongc, on="name")
m
m.to_csv(prefix+"/gcskewdb.csv", float_format='%g')
print(len(m))


# In[108]:


plt.figure()
logbins = np.geomspace(0.005, m.rmsGC.max(), 100)

plt.hist(m.rmsGC, bins=logbins, cumulative=True, density=True)
plt.xscale('log')


# In[109]:





# In[146]:


#chromoname="NZ_CP032359.1"  # this is one of the neatest in the whole database
chromoname="NZ_CP012122.1" # also very neat
#chromoname="NZ_AP023049.1" #test
#chromoname="NZ_CP044177.1" # anomalous
print(chromoname)
fitted = pandas.read_csv(prefix+"/"+chromoname+"_fit.csv")

plt.figure()
us=results[results.name==chromoname]
plt.plot(fitted.pos, fitted.gcskew, label="Cumulative GC skew")
plt.plot(fitted.pos, fitted.predgcskew, label="Fitted GC skew")

leshift=us["shift"].item()
if leshift > 0:
    plt.axvline(leshift, ls=':', color='black')
else:
    plt.axvline(fitted.pos.max() + leshift, ls=':', color='black')

plt.xlabel("Locus")
plt.ylabel("Skew")
plt.legend()
plt.title(chromoname + " alpha1 " + str(round(us["alpha1gc"].item(),3)) + " alpha2 " + str(round(us["alpha2gc"].item(),3)) + 
          " div " + str(round(us["div"].item(),3)) + " rmsGC "+str(round(us["rmsGC"].item(),3)) )
plt.grid()
print(fitted.tail(1).gcskew.item())
plt.hlines([0, fitted.tail(1).gcskew.item()], 
           [0, fitted.tail(1).pos.item()] ,
           [fitted.pos.max() * m[m.name==chromoname]["div"].item(), 
            fitted.pos.max() * m[m.name==chromoname]["div"].item()])

plt.annotate('alpha1',
            xy=(190000, 1500), xytext=(190000, 1500), xycoords='data')

plt.annotate('alpha2',
            xy=(fitted.tail(1).pos.item() - 420000, fitted.tail(1).gcskew.item()+1400),
             xytext=(fitted.tail(1).pos.item() - 420000, fitted.tail(1).gcskew.item()+1400),
             xycoords='data')
plt.savefig("test.svg")


#plt.annotate("",
#            xy=(105000, 5000), xycoords='data',
#            xytext=(170000, 0), textcoords='data',
#            arrowprops=dict(arrowstyle="-",
#                            connectionstyle="arc3,rad=0.3")
#            )

plt.savefig("explainer.svg") # needs to be postprocessed with Inkscape, explainer.svg


# In[111]:





# In[124]:


chromoname="NZ_CP044177.1" # anomalous
print(chromoname)
fitted = pandas.read_csv(prefix+"/"+chromoname+"_fit.csv")
results[results.name==chromoname]
plt.figure()
us=results[results.name==chromoname]
plt.plot(fitted.pos, fitted.gcskew, label="Cumulative GC skew")
plt.plot(fitted.pos, fitted.predgcskew, label="Fitted GC skew")
plt.plot(fitted.pos, fitted.taskew, label="Cumulative TA skew")
plt.plot(fitted.pos, fitted.predtaskew, label="Fitted TA skew")

#plt.plot(fitted.pos, fitted.sbskew, label="Cumulative SB skew")
#plt.plot(fitted.pos, fitted.predsbskew, label="Fitted SB skew")


leshift=us["shift"].item()
if leshift > 0:
    plt.axvline(leshift, ls=':', color='black')
else:
    plt.axvline(fitted.pos.max() + leshift, ls=':', color='black')

plt.xlabel("Locus")
plt.ylabel("Skew")
plt.legend()
plt.title(chromoname + " alpha1 " + str(round(us["alpha1gc"].item(),3)) + " alpha2 " + str(round(us["alpha2gc"].item(),3)) + 
          " div " + str(round(us["div"].item(),3)) + " rmsGC "+str(round(us["rmsGC"].item(),3)) )
plt.grid()
plt.savefig("anomalous.pdf")


# In[65]:





# In[143]:


chromoname="NZ_CP019870.1" # c difficile
fitted = pandas.read_csv(prefix+"/"+chromoname+"_fit.csv")
plt.figure()
us=results[results.name==chromoname]
print(us.alpha1ta)
plt.plot(fitted.pos, fitted.gcskew, label="Cumulative GC skew")
plt.plot(fitted.pos, fitted.predgcskew, label="Fitted GC skew")
plt.plot(fitted.pos, fitted.taskew, label="Cumulative TA skew")
plt.plot(fitted.pos, fitted.predtaskew, label="Fitted TA skew")

#plt.plot(fitted.pos, fitted.sbskew, label="Cumulative SB skew")
#plt.plot(fitted.pos, fitted.predsbskew, label="Fitted SB skew")


leshift=us["shift"].item()
if leshift > 0:
    plt.axvline(leshift, ls=':', color='black')
else:
    plt.axvline(fitted.pos.max() + leshift, ls=':', color='black')

plt.xlabel("Locus")
plt.ylabel("Skew")
plt.legend()
plt.title(chromoname + " alpha1 " + str(round(us["alpha1gc"].item(),3)) + " alpha2 " + str(round(us["alpha2gc"].item(),3)) + 
          " div " + str(round(us["div"].item(),3)) + " rmsGC "+str(round(us["rmsGC"].item(),3)) )
plt.grid()
plt.savefig("cdif.pdf")


# In[144]:


plt.figure()

leX=fitted.sbskew.diff().rolling(10, center=True).mean()/4096
leY=fitted.gcskew.diff().rolling(10, center=True).mean()/4096


frame = { 'predleading': fitted.predleading, 'x': leX, 'y': leY }
result = pandas.DataFrame(frame)
#print(result.tail(20))


sub=result[(result.predleading==1) & ((result.x < 0) | (result.x>0)) & ((result.y < 0) | (result.y > 0)) ]
#z = np.polyfit(sub.x, sub.y, 1)
#p = np.poly1d(z)
#xp = np.linspace(leX.min(), leX.max(), 100)
#plt.plot(xp, p(xp), ':', color='red')

model = smf.quantreg('y ~ x', sub).fit(q=0.5)
print(model.summary())

print(model.params['Intercept'], model.params['x'])

get_y = lambda a, b: a + b * sub.x
y = get_y(model.params['Intercept'], model.params['x'])
plt.plot(sub.x, y, color='black')


sub=result[(result.predleading==0) & ((result.x < 0) | (result.x>0)) & ((result.y < 0) | (result.y > 0))]

#z = np.polyfit(sub.x, sub.y, 1)
#p = np.poly1d(z)
#xp = np.linspace(leX.min(), leX.max(), 100)
#plt.plot(xp, p(xp), ':', color='red')

model = smf.quantreg('y ~ x', sub).fit(q=0.5)
print(model.summary())
print(model.params['Intercept'], model.params['x'])
get_y = lambda a, b: a + b * sub.x
y = get_y(model.params['Intercept'], model.params['x'])
plt.plot(sub.x, y, color='black')


plt.scatter(result[result.predleading==0].x, result[result.predleading==0].y, marker="+", color='#1f77b4', label="GC lagging")
plt.scatter(result[result.predleading==1].x, result[result.predleading==1].y, marker="^", color='#1f77b4', label="GC leading")

leX=fitted.sbskew.diff().rolling(10, center=True).mean()/4096
leY=fitted.taskew.diff().rolling(10, center=True).mean()/4096

z = np.polyfit(leX.dropna(), leY.dropna(), 1)
p = np.poly1d(z)
xp = np.linspace(leX.min(), leX.max(), 100)
plt.plot(xp, p(xp), ':', color='red')

plt.scatter(leX, leY,  s=10, label="TA", color='#ff7f0e')
plt.xlabel("SB skew")
plt.ylabel("GC/TA skew")
plt.title(chromoname)
plt.legend()


plt.grid()
plt.savefig("cdif-histo.pdf")


# In[49]:


plt.figure()
plt.hist(m["div"], bins=50, density=True)
plt.grid()
plt.title("Division between leading/lagging strand")


# In[155]:


# how many chromosomes show one strand being 3 times smaller than the other, for good fit quality?
m[(m.rmsGC <0.1) & ((m["div"] <0.25) | (m["div"] > 0.75))]["div"].describe()


# In[163]:


# how many chromosomes show one strand having 3 times more/less skew than the ohther?
m[(m.rmsGC <0.1) & ((m["alpha1gc"]/m["alpha2gc"] < 0.33) | (m["alpha1gc"]/m["alpha2gc"] > 3))].alpha1gc


# In[50]:


plt.figure()
plt.scatter(results.alpha1gc, results.alpha2gc)
plt.grid()


# In[51]:


plt.figure()
plt.hist((results.alpha1gc-results.alpha2gc)/results.alpha1gc, range=(-0.5, 0.5), bins=20)
plt.grid()


# In[941]:


results[results["div"] < 0.35]


# In[225]:


# you can pick if you want flat genomes:
sel=m[(m.rmsGC < 0.035) & (
        ((m["alpha1gc"] < 0.0014) & (m["alpha1gc"] > 0))
    ^
       ((m["alpha2gc"] < 0.0014) & (m["alpha2gc"] > 0))
    ) ]  #  & (m.realm3=="Firmicutes")

#or unequally distributed ones:
#sel=m[(m["div"] < 0.2) & (m.rmsGC < 0.002) ]
print(len(sel))
d=2
fig,axs=plt.subplots(d,d, sharex=False, sharey=False, constrained_layout=True)
a=0
b=0
names=sel.name.unique()
print(len(names))
shuffle(names)
for k in names:
        print(a,b,k)
        fitted = pandas.read_csv(prefix+"/"+k+"_fit.csv")
        t=results[results.name==k]
        #print(t.alpha1gc, t.minpos, t.maxpos)
        axs[b,a].plot(fitted.pos, fitted.gcskew)
        axs[b,a].plot(fitted.pos, fitted.predgcskew)
        #axs[b,a].plot(fitted.pos, fitted.taskew)
        #axs[b,a].plot(fitted.pos, fitted.predtaskew)

        #axs[b,a].get_yaxis().set_ticks([])
        axs[b,a].set_title(k) #  + " "+str(1000*results[results.name==k].rmsGC.mean()), 
        #axs[b,a].grid()
        a=a+1
        if(a>=d):
            a=0
            b=b+1
        if(b == d):
            break
#fig.suptitle("GC/TA skew in random bacterial chromosomes + fit")
plt.savefig("flat-skew.pdf")


# In[224]:


# you can pick if you want flat genomes:
sel=m[((m["div"] < 0.20) | (m["div"] > 0.8)) & (m.rmsGC < 0.035) ]
print(len(sel))
d=2
fig,axs=plt.subplots(d,d, sharex=False, sharey=False, constrained_layout=True)
a=0
b=0
names=sel.name.unique()
print(len(names))
shuffle(names)
for k in names:
        print(a,b,k)
        fitted = pandas.read_csv(prefix+"/"+k+"_fit.csv")
        t=results[results.name==k]
        #print(t.alpha1gc, t.minpos, t.maxpos)
        axs[b,a].plot(fitted.pos, fitted.gcskew)
        axs[b,a].plot(fitted.pos, fitted.predgcskew)
        #axs[b,a].plot(fitted.pos, fitted.taskew)
        #axs[b,a].plot(fitted.pos, fitted.predtaskew)

        #axs[b,a].get_yaxis().set_ticks([])
        axs[b,a].set_title(k) #  + " "+str(1000*results[results.name==k].rmsGC.mean()), 
        #axs[b,a].grid()
        a=a+1
        if(a>=d):
            a=0
            b=b+1
        if(b == d):
            break
#fig.suptitle("GC/TA skew in random bacterial chromosomes + fit")
plt.savefig("strand-div.pdf")


# In[220]:


# you can pick if you want flat genomes:
d=4
fig,axs=plt.subplots(d,d, sharex=False, sharey=False, constrained_layout=True)
a=0
b=0

res,edges=pandas.qcut(m.rmsGC, 16, retbins=True)
edges
for q in edges:
        sel = m[m.rmsGC>q]
        k=sel.sort_values(["rmsGC"]).head(1).name.item()
        fitted = pandas.read_csv(prefix+"/"+k+"_fit.csv")
        t=results[results.name==k]
        #print(t.alpha1gc, t.minpos, t.maxpos)
        axs[b,a].plot(fitted.pos, fitted.gcskew)
        axs[b,a].plot(fitted.pos, fitted.predgcskew)
        axs[b,a].get_xaxis().set_visible(False)
        axs[b,a].get_yaxis().set_visible(False)

        #axs[b,a].plot(fitted.pos, fitted.taskew)
        #axs[b,a].plot(fitted.pos, fitted.predtaskew)

        #axs[b,a].get_yaxis().set_ticks([])
        axs[b,a].set_title(("rmsGC = %.4f" % q )) #  + " "+str(1000*results[results.name==k].rmsGC.mean()), 
        #axs[b,a].grid()
        a=a+1
        if(a>=d):
            a=0
            b=b+1
        if(b == d):
            break
#fig.suptitle("GC/TA skew fit for 16 equally sized quality limits")
plt.savefig("rms-samples.pdf")


# In[223]:


len(m[m.rmsGC<0.1])/len(m)


# In[79]:


for k in m.groupby(["realm2"]).name.count().reset_index().sort_values(["name"], ascending=False).realm2.head(10):
    print(k)


# In[115]:


plt.figure()
for k in m.groupby(["realm2"]).name.count().reset_index().sort_values(["name"], ascending=False).realm2.head(5):
    sel=m[m.realm2==k]
    plt.scatter(sel.alpha1gc, sel.alpha1ta, alpha=0.2, label=k)
    
plt.xlabel("alpha1 of GC skew")
plt.ylabel("alpha1 of TA skew")
leg=plt.legend()
for lh in leg.legendHandles: 
    lh.set_alpha(1)
plt.grid()
plt.xlim((0,0.09))
plt.savefig("phylo-histo.png", dpi=300)


# In[118]:


firmi=m[(m.realm3=="Firmicutes") & (m.rmsGC < 0.3) & (m.rmsTA < 0.3)]
plt.figure()
#leX = firmi.acounts2/(firmi.acounts2 + firmi.ccounts2 + firmi.gcounts2 + firmi.tcounts2)
#leX = firmi.gccount/firmi.siz
leX =  -(firmi.cfrac - firmi.gfrac)*(1-firmi.ngcount/firmi.siz)*firmi.alpha2sb
#leX =  firmi.alpha1sb

leY=firmi.alpha1gc
z1 = np.polyfit(leX, leY, 1)
print(z1)
plt.scatter(leX, leY, s=1.5, label="GC")
p = np.poly1d(z1)
xp = np.linspace(leX.min(), leX.max(), 100)
plt.plot(xp, p(xp), color="red")

#

leX =  -(firmi.tfrac - firmi.afrac)*(1-firmi.ngcount/firmi.siz)*firmi.alpha2sb
#leX =  (firmi.tfrac-0.5)*firmi.alpha1sb

leY=firmi.alpha1ta

z2 = np.polyfit(leX, leY, 1)
print(z2)
plt.scatter(leX, leY, s=1.5, label="TA")
p = np.poly1d(z2)
xp = np.linspace(leX.min(), leX.max(), 100)
plt.plot(xp, p(xp), color="red")

plt.grid()
#plt.scatter(firmi.ccounts2-firmi.tcounts2, firmi.alpha1gc)
print(z1[0]/z2[0])
plt.legend()


# In[226]:


firmi=m[(m.realm3=="Firmicutes") & (m.rmsGC < 0.25) & (m.rmsTA < 0.25)]
plt.figure()
leX =  -(firmi.cfrac - firmi.gfrac)*(1-firmi.ngcount/firmi.siz)*firmi.alpha2sb

leY=firmi.alpha1gc

frame = {  'x': leX, 'y': leY }
sub = pandas.DataFrame(frame)

model = smf.quantreg('y ~ x', sub).fit(q=0.5)
print(model.summary())

print(model.params['Intercept'], model.params['x'])

get_y = lambda a, b: a + b * sub.x
y = get_y(model.params['Intercept'], model.params['x'])

plt.scatter(leX, leY, s=2, label="GC")
plt.plot(sub.x, y, color='black')

#

leX =  -(firmi.tfrac - firmi.afrac)*(1-firmi.ngcount/firmi.siz)*firmi.alpha2sb
leY=firmi.alpha1ta

frame = {  'x': leX, 'y': leY }
sub = pandas.DataFrame(frame)

model = smf.quantreg('y ~ x', sub).fit(q=0.5)
print(model.summary())

print(model.params['Intercept'], model.params['x'])

get_y = lambda a, b: a + b * sub.x
y = get_y(model.params['Intercept'], model.params['x'])
plt.scatter(leX, leY, s=2, label="TA")

plt.plot(sub.x, y, color='black', label='fit')





plt.grid()
#plt.scatter(firmi.ccounts2-firmi.tcounts2, firmi.alpha1gc)
plt.legend()
plt.xlabel("Product of gene strand bias, codon bias skew, percentage coding")
plt.ylabel("GC/TA skew fraction")
#plt.title("Data for "+str(len(firmi))+ " Firmicute chromosomes")
plt.savefig("firmi.pdf")


# In[120]:


firmi=m[(m.realm3=="Firmicutes") & (m.rmsGC < 0.2) ]
plt.figure()
#leX = firmi.acounts2/(firmi.acounts2 + firmi.ccounts2 + firmi.gcounts2 + firmi.tcounts2)
#leX = firmi.gccount/firmi.siz
leX =  -(firmi.cfrac - firmi.gfrac)*(1-firmi.ngcount/firmi.siz)*firmi.alpha1sb


leY=firmi.alpha1gc


z1 = np.polyfit(leX, leY, 1)
print(z1)
plt.scatter(leX, leY, s=1.5, label="GC")
p = np.poly1d(z1)
xp = np.linspace(leX.min(), leX.max(), 100)
plt.plot(xp, p(xp), color="red")

#

leX =  -(firmi.tfrac - firmi.afrac)*(1-firmi.ngcount/firmi.siz)*firmi.alpha1sb
leY=firmi.alpha1ta

z2 = np.polyfit(leX, leY, 1)
print(z2)
plt.scatter(leX, leY, s=1.5, label="TA")
p = np.poly1d(z2)
xp = np.linspace(leX.min(), leX.max(), 100)
plt.plot(xp, p(xp), color="red")

plt.grid()
#plt.scatter(firmi.ccounts2-firmi.tcounts2, firmi.alpha1gc)
print(z1[0]/z2[0])
plt.legend()


# In[ ]:





