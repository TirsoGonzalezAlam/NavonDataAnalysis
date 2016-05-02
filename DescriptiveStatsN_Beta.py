# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 13:36:26 2016

@author: trga500
"""

#CHANGELOG
# do normal tests
# add sd
# add plots

#Import things we need
import os, csv
import numpy as np
import scipy
from scipy import stats
from matplotlib import pyplot

#Implement absolute path and change directory to the script's location
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

#Load the directory to analyse
parent = os.listdir(dname)

#Prepare csv file for output
f = open(r'C:\Users\trga500\Desktop\Data\Collection\Week 6\Navon\results\WeeklyAnalysis.csv', 'wb+')
filewriter = csv.writer(f, delimiter=',')
#filewriter.writerow(['Part_ID', 'Word_RT_Mean', 'Word_SD', 'Word_Acc', 'Pic_RT_Mean', 'Pic_SD', 'Pic_Acc', 'Easy_RT_Mean', 'Easy_SD', 'Easy_Acc', 'Hard_RT_Mean', 'Hard_SD', 'Hard_Acc'])
f.write('Part_ID, Congruent_RT_Mean, Incongruent_RT_Mean, Congruent_ACC, Incongruent_ACC, Congruent_EFF, Incongruent_EFF, Global_RT_Mean, Local_RT_Mean, Global_ACC, Local_ACC, Global_EFF, Local_EFF, Global_Congruent_RT, Global_Incongruent_RT, Local_Congruent_RT, Local_Incongruent_RT, Global_Congruent_ACC, Global_Incongruent_ACC, Local_Congruent_ACC, Local_Incongruent_ACC, Global_Congruent_EFF, Global_Incongruent_EFF, Local_Congruent_EFF, Local_Incongruent_EFF \n')


#Initialise empty variables and lists we need
g_cong_rt = []
g_incong_rt = []
g_glob_rt = []
g_loc_rt = []
g_cong_acc = []
g_incong_acc = []
g_glob_acc = []
g_loc_acc = []
g_cong_eff = []
g_incong_eff = []
g_glob_eff = []
g_loc_eff = []

g_glob_cong_RT_mean = []
g_glob_cong_ACC_mean = []
g_glob_cong_EFF_mean = []
g_glob_incong_RT_mean = []
g_glob_incong_ACC_mean = []
g_glob_incong_EFF_mean = []
g_loc_cong_RT_mean = []
g_loc_cong_ACC_mean = []
g_loc_cong_EFF_mean = []
g_loc_incong_RT_mean = []
g_loc_incong_ACC_mean = []
g_loc_incong_EFF_mean = []

cong_mean = []
incong_mean = []
glob_mean = []
loc_mean = []
acc_cong_mean = []
acc_incong_mean = []
acc_glob_mean = []
acc_loc_mean = []
glob_cong_mean = []
glob_incong_mean = []
loc_cong_mean = []
loc_incong_mean = []
acc_glob_cong_mean = []
acc_glob_incong_mean = []
acc_loc_cong_mean = []
acc_loc_incong_mean = []

g_cong_counter = 0
g_incong_counter = 0
g_glob_counter = 0
g_loc_counter = 0
g_glob_cong_counter = 0
g_glob_incong_counter = 0
g_loc_cong_counter = 0
g_loc_incong_counter = 0



#Now we fill in the rows
for i in parent:
    cong_rt = []
    incong_rt = []
    glob_rt = []
    loc_rt = []
    cong_acc = []
    incong_acc = []
    glob_acc = []
    loc_acc = []
    glob_cong_rt = []
    glob_incong_rt = []
    loc_cong_rt = []
    loc_incong_rt = []
    glob_cong_acc = []
    glob_incong_acc = []
    loc_cong_acc = []
    loc_incong_acc = []
    cong_counter = 0
    incong_counter = 0
    glob_counter = 0
    loc_counter = 0
    glob_cong_counter = 0
    glob_incong_counter = 0
    loc_cong_counter = 0
    loc_incong_counter = 0
    if i.endswith('.csv'):
        g = open(i)
        for line in g:
            #Prepare the lines
            line2 = line.strip()
            data = line2.split(',')
            if data[0] == 'Block':
                continue
            elif 'Done' in data[0]:
                continue
            else:
                
                #Load the RT and ACC data into the appropriate list  
                #First, the interaction
            
                #Global Congruent and Incongruent
                if data[0] == 'glob1' or data[0] == 'glob2':
                    glob_counter = glob_counter + 1
                    glob_acc.append(data[4])
                    if data[5] != 'nan':
                        glob_rt.append(float(data[5]))

                    #For the group statistics
                    if data[1] == 'congruent':
                        g_glob_cong_counter = g_glob_cong_counter + 1
                        glob_cong_acc.append(data[4])

                    elif data[1] == 'incongruent':
                        g_glob_incong_counter = g_glob_incong_counter + 1
                        glob_incong_acc.append(data[4])
                            
                    #and for the participant statistics
                    if data[1] == 'congruent':
                        glob_cong_counter = glob_cong_counter + 1
                        if data[5] != 'nan':
                            glob_cong_rt.append(float(data[5]))

                    elif data[1] == 'incongruent':
                        glob_incong_counter = glob_incong_counter + 1
                        if data[5] != 'nan':
                            glob_incong_rt.append(float(data[5]))
                                     
                #Local Congruent and Incongruent
                if data[0] == 'loc1' or data[0] == 'loc2':
                    loc_counter = loc_counter + 1
                    loc_acc.append(data[4])
                    if data[5] != 'nan':
                        loc_rt.append(float(data[5]))

                    #For the group statistics
                    if data[1] == 'congruent':
                        g_loc_cong_counter = g_loc_cong_counter + 1

                    elif data[1] == 'incongruent':
                        g_loc_incong_counter = g_loc_incong_counter + 1
                            
                    #and for the participant statistics
                    if data[1] == 'congruent':
                        loc_cong_counter = loc_cong_counter + 1
                        loc_cong_acc.append(data[4])
                        if data[5] != 'nan':
                            loc_cong_rt.append(float(data[5]))

                    elif data[1] == 'incongruent':
                        loc_incong_counter = loc_incong_counter + 1
                        loc_incong_acc.append(data[4])
                        if data[5] != 'nan':
                            loc_incong_rt.append(float(data[5]))
 #ADD JUST CONGRUENT
                #Just Congruent
                if data[1] == 'congruent':
                    cong_counter = cong_counter + 1
                    g_cong_counter = g_cong_counter + 1
                    cong_acc.append(data[4])
                    if data[5] != 'nan':
                        cong_rt.append(float(data[5]))

                elif data[1] == 'incongruent':
                    incong_counter = incong_counter + 1
                    g_incong_counter = g_incong_counter + 1
                    incong_acc.append(data[4])
                    if data[5] != 'nan':
                        incong_rt.append(float(data[5])) 
        g.close()

        #Once we have this data, we can calculate individual results
        #Calculate mean and %correct for individual participant
        #Calculate ID
        part_id = i
        
        """Calculate Accuracy"""
        #Global Congruent
        x = stats.itemfreq(glob_cong_acc)
        try:
            corr_perc_GlobCong = (float(x[1][1])/glob_cong_counter)*100
        except IndexError:
            corr_perc_GlobCong = 100

        #Global Incongruent
        y = stats.itemfreq(glob_incong_acc)
        try:
            corr_perc_GlobIncong = (float(y[1][1])/glob_incong_counter)*100
        except IndexError:
            corr_perc_GlobIncong = 100
        
        #Local Congruent
        xx = stats.itemfreq(loc_cong_acc)
        try:
            corr_perc_LocCong = (float(xx[1][1])/loc_cong_counter)*100
        except IndexError:
            corr_perc_LocCong = 100

        #Local Incongruent
        yy = stats.itemfreq(loc_incong_acc)
        try:
            corr_perc_LocIncong = (float(yy[1][1])/loc_incong_counter)*100
        except IndexError:
            corr_perc_LocIncong = 100       
        
        #Global
        aa = stats.itemfreq(glob_acc)
        try:
            corr_perc_Glob = (float(aa[1][1])/glob_counter)*100
        except IndexError:
            corr_perc_Glob = 100

        #Local
        bb = stats.itemfreq(loc_acc)
        try:
            corr_perc_Loc = (float(bb[1][1])/loc_counter)*100
        except IndexError:
            corr_perc_Loc = 100

        #Congruent        
        cc = stats.itemfreq(cong_acc)
        try:
            corr_perc_Cong = (float(cc[1][1])/cong_counter)*100
        except IndexError:
            corr_perc_Cong = 100

        #Incongruent
        dd = stats.itemfreq(incong_acc)
        try:
            corr_perc_Incong = (float(dd[1][1])/incong_counter)*100
        except IndexError:
            corr_perc_Incong = 100
            
        """Calculate RT's"""
        #Global Congruent
        mean_rt_GlobCong = np.mean(glob_cong_rt)
        #Global Incongruent
        mean_rt_GlobIncong = np.mean(glob_incong_rt)
        #Local Congruent
        mean_rt_LocCong = np.mean(loc_cong_rt)
        #Local Incongruent
        mean_rt_LocIncong = np.mean(loc_incong_rt)
        #Global
        mean_rt_Glob = np.mean(glob_rt)
        #Local
        mean_rt_Loc = np.mean(loc_rt)
        #Congruent
        mean_rt_Cong = np.mean(cong_rt)
        #Incongruent
        mean_rt_Incong = np.mean(incong_rt)
        
        """Calculate Efficiency"""
        GlobCong_eff = mean_rt_GlobCong/(corr_perc_GlobCong/100)
        GlobIncong_eff = mean_rt_GlobIncong/(corr_perc_GlobIncong/100)
        LocCong_eff = mean_rt_LocCong/(corr_perc_LocCong/100)
        LocIncong_eff = mean_rt_LocIncong/(corr_perc_LocIncong/100)
        Glob_eff = mean_rt_Glob/(corr_perc_Glob/100)
        Loc_eff = mean_rt_Loc/(corr_perc_Loc/100)
        Cong_eff = mean_rt_Cong/(corr_perc_Cong/100)
        Incong_eff = mean_rt_Incong/(corr_perc_Incong/100)
        
        """Prepare individual participant's results for group stats"""
        #Accuracy
        g_glob_cong_ACC_mean.append(corr_perc_GlobCong)
        g_glob_incong_ACC_mean.append(corr_perc_GlobIncong)
        g_loc_cong_ACC_mean.append(corr_perc_LocCong)
        g_loc_incong_ACC_mean.append(corr_perc_LocIncong)
        g_glob_acc.append(corr_perc_Glob)
        g_loc_acc.append(corr_perc_Loc)
        g_cong_acc.append(corr_perc_Cong)
        g_incong_acc.append(corr_perc_Incong)
        #RT
        g_glob_cong_RT_mean.append(mean_rt_GlobCong)
        g_glob_incong_RT_mean.append(mean_rt_GlobIncong)
        g_loc_cong_RT_mean.append(mean_rt_LocCong)
        g_loc_incong_RT_mean.append(mean_rt_LocIncong)
        g_glob_rt.append(mean_rt_Glob)
        g_loc_rt.append(mean_rt_Loc)
        g_cong_rt.append(mean_rt_Cong)
        g_incong_rt.append(mean_rt_Incong)
        #Efficiency
        g_glob_cong_EFF_mean.append(GlobCong_eff)
        g_glob_incong_EFF_mean.append(GlobIncong_eff)
        g_loc_cong_EFF_mean.append(LocCong_eff)
        g_loc_incong_EFF_mean.append(LocIncong_eff)
        g_glob_eff.append(Glob_eff)
        g_loc_eff.append(Loc_eff)
        g_cong_eff.append(Cong_eff)
        g_incong_eff.append(Incong_eff)

        f.write('%s, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n'%(part_id[5:7].rstrip('_'), mean_rt_Cong, mean_rt_Incong, corr_perc_Cong, corr_perc_Incong, Cong_eff, Incong_eff, mean_rt_Glob, mean_rt_Loc, corr_perc_Glob, corr_perc_Loc, Glob_eff, Loc_eff, mean_rt_GlobCong, mean_rt_GlobIncong, mean_rt_LocCong, mean_rt_LocIncong, corr_perc_GlobCong, corr_perc_GlobIncong, corr_perc_LocCong, corr_perc_LocIncong, GlobCong_eff, GlobIncong_eff, LocCong_eff, LocIncong_eff))

#Calculate z scores per participant
#RT
z_score_rt_c = scipy.stats.zscore(g_cong_rt)
z_score_rt_i = scipy.stats.zscore(g_incong_rt)
z_score_rt_l = scipy.stats.zscore(g_loc_rt)
z_score_rt_g = scipy.stats.zscore(g_glob_rt)
z_score_rt_gc = scipy.stats.zscore(g_glob_cong_RT_mean)
z_score_rt_gi = scipy.stats.zscore(g_glob_incong_RT_mean)
z_score_rt_lc = scipy.stats.zscore(g_loc_cong_RT_mean)
z_score_rt_li = scipy.stats.zscore(g_loc_incong_RT_mean)
#ACC
z_score_acc_c = scipy.stats.zscore(g_cong_acc)
z_score_acc_i = scipy.stats.zscore(g_incong_acc)
z_score_acc_l = scipy.stats.zscore(g_loc_acc)
z_score_acc_g = scipy.stats.zscore(g_glob_acc)
z_score_acc_gc = scipy.stats.zscore(g_glob_cong_ACC_mean)
z_score_acc_gi = scipy.stats.zscore(g_glob_incong_ACC_mean)
z_score_acc_lc = scipy.stats.zscore(g_loc_cong_ACC_mean)
z_score_acc_li = scipy.stats.zscore(g_loc_incong_ACC_mean)
#EFF
z_score_eff_c = scipy.stats.zscore(g_cong_eff)
z_score_eff_i = scipy.stats.zscore(g_incong_eff)
z_score_eff_l = scipy.stats.zscore(g_glob_eff)
z_score_eff_g = scipy.stats.zscore(g_loc_eff)
z_score_eff_gc = scipy.stats.zscore(g_glob_cong_EFF_mean)
z_score_eff_gi = scipy.stats.zscore(g_glob_incong_EFF_mean)
z_score_eff_lc = scipy.stats.zscore(g_loc_cong_EFF_mean)
z_score_eff_li = scipy.stats.zscore(g_loc_incong_EFF_mean)


#Group Statistics
#Global Congruent
Group_Global_Congruent_RT = np.mean(g_glob_cong_RT_mean)
Group_Global_Congruent_ACC = np.mean(g_glob_cong_ACC_mean)
Group_Global_Congruent_EFF = np.mean(g_glob_cong_EFF_mean)
StdErr_Global_Congruent_RT = stats.sem(g_glob_cong_RT_mean)
StdErr_Global_Congruent_ACC = stats.sem(g_glob_cong_ACC_mean)
StdErr_Global_Congruent_EFF = stats.sem(g_glob_cong_EFF_mean)
#Global Incongruent
Group_Global_Incongruent_RT = np.mean(g_glob_incong_RT_mean)
Group_Global_Incongruent_ACC = np.mean(g_glob_incong_ACC_mean)
Group_Global_Incongruent_EFF = np.mean(g_glob_incong_EFF_mean)
StdErr_Global_Incongruent_RT = stats.sem(g_glob_incong_RT_mean)
StdErr_Global_Incongruent_ACC = stats.sem(g_glob_incong_ACC_mean)
StdErr_Global_Incongruent_EFF = stats.sem(g_glob_incong_EFF_mean)
#Local Congruent
Group_Local_Congruent_RT = np.mean(g_loc_cong_RT_mean)
Group_Local_Congruent_ACC = np.mean(g_loc_cong_ACC_mean)
Group_Local_Congruent_EFF = np.mean(g_loc_cong_EFF_mean)
StdErr_Local_Congruent_RT = stats.sem(g_loc_cong_RT_mean)
StdErr_Local_Congruent_ACC = stats.sem(g_loc_cong_ACC_mean)
StdErr_Local_Congruent_EFF = stats.sem(g_loc_cong_EFF_mean)
#Local Incongruent
Group_Local_Incongruent_RT = np.mean(g_loc_incong_RT_mean)
Group_Local_Incongruent_ACC = np.mean(g_loc_incong_ACC_mean)
Group_Local_Incongruent_EFF = np.mean(g_loc_incong_EFF_mean)
StdErr_Local_Incongruent_RT = stats.sem(g_loc_incong_RT_mean)
StdErr_Local_Incongruent_ACC = stats.sem(g_loc_incong_ACC_mean)
StdErr_Local_Incongruent_EFF = stats.sem(g_loc_incong_EFF_mean)
#Global
Group_Global_RT = np.mean(g_glob_rt)
Group_Global_ACC = np.mean(g_glob_acc)
Group_Global_EFF = np.mean(g_glob_eff)
StdErr_Global_RT = stats.sem(g_glob_rt)
StdErr_Global_ACC = stats.sem(g_glob_acc)
StdErr_Global_EFF = stats.sem(g_glob_eff)
#Local
Group_Local_RT = np.mean(g_loc_rt)
Group_Local_ACC = np.mean(g_loc_acc)
Group_Local_EFF = np.mean(g_loc_eff)
StdErr_Local_RT = stats.sem(g_loc_rt)
StdErr_Local_ACC = stats.sem(g_loc_acc)
StdErr_Local_EFF = stats.sem(g_loc_eff)
#Congruent
Group_Congruent_RT = np.mean(g_cong_rt)
Group_Congruent_ACC = np.mean(g_cong_acc)
Group_Congruent_EFF = np.mean(g_cong_eff)
StdErr_Congruent_RT = stats.sem(g_cong_rt)
StdErr_Congruent_ACC = stats.sem(g_cong_acc)
StdErr_Congruent_EFF = stats.sem(g_cong_eff)
#Incongruent
Group_Incongruent_RT = np.mean(g_incong_rt)
Group_Incongruent_ACC = np.mean(g_incong_acc)
Group_Incongruent_EFF = np.mean(g_incong_eff)
StdErr_Incongruent_RT = stats.sem(g_incong_rt)
StdErr_Incongruent_ACC = stats.sem(g_incong_acc)
StdErr_Incongruent_EFF = stats.sem(g_incong_eff)

#Do Mean Comparisons
#RTs
RT_CongVsIncong = scipy.stats.ttest_rel(g_cong_rt, g_incong_rt)
RT_GlobVsLoc = scipy.stats.ttest_rel(g_glob_rt, g_loc_rt)
#ACC
ACC_CongVsIncong = scipy.stats.ttest_rel(g_cong_acc, g_incong_acc)
ACC_GlobVsLoc = scipy.stats.ttest_rel(g_glob_acc, g_loc_acc)
#EFF
EFF_CongVsIncong = scipy.stats.ttest_rel(g_cong_eff, g_incong_eff)
EFF_GlobVsLoc = scipy.stats.ttest_rel(g_glob_eff, g_loc_eff)

#ANOVA
RT_Interaction = scipy.stats.f_oneway(g_glob_cong_RT_mean, \
g_glob_incong_RT_mean, g_loc_cong_RT_mean, g_loc_incong_RT_mean)
ACC_Interaction = scipy.stats.f_oneway(g_glob_cong_ACC_mean, \
g_glob_incong_ACC_mean, g_loc_cong_ACC_mean, g_loc_incong_ACC_mean)
EFF_Interaction = scipy.stats.f_oneway(g_glob_cong_EFF_mean, \
g_glob_incong_EFF_mean, g_loc_cong_EFF_mean, g_loc_incong_EFF_mean)

#Write results to the file
f.write('\n')
f.write('Group Statistics:\n')
f.write('GROUP ANALYSIS, Congruent_RT_Mean, Incongruent_RT_Mean, Congruent_ACC, Incongruent_ACC, Congruent_EFF, Incongruent_EFF, Global_RT_Mean, Local_RT_Mean, Global_ACC, Local_ACC, Global_EFF, Local_EFF, Global_Congruent_RT, Global_Incongruent_RT, Local_Congruent_RT, Local_Incongruent_RT, Global_Congruent_ACC, Global_Incongruent_ACC, Local_Congruent_ACC, Local_Incongruent_ACC, Global_Congruent_EFF, Global_Incongruent_EFF, Local_Congruent_EFF, Local_Incongruent_EFF \n')
f.write('%s, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n'%('GROUP', Group_Congruent_RT, Group_Incongruent_RT, Group_Congruent_ACC, Group_Incongruent_ACC, Group_Congruent_EFF, Group_Incongruent_EFF, Group_Global_RT, Group_Local_RT, Group_Global_ACC, Group_Local_ACC, Group_Global_EFF, Group_Local_EFF, Group_Global_Congruent_RT, Group_Global_Incongruent_RT, Group_Local_Congruent_RT, Group_Local_Incongruent_RT, Group_Global_Congruent_ACC, Group_Global_Incongruent_ACC, Group_Local_Congruent_ACC, Group_Local_Incongruent_ACC, Group_Global_Congruent_EFF, Group_Global_Incongruent_EFF, Group_Local_Congruent_EFF, Group_Local_Incongruent_EFF))


#Write the z scores results in the results file
f.write('\n')
f.write('z Scores \n')
f.write('ID, z_Congruent_RT_Mean, z_Incongruent_RT_Mean, z_Congruent_ACC, \z_Incongruent_ACC, z_Congruent_EFF, z_Incongruent_EFF, z_Global_RT_Mean, z_Local_RT_Mean, z_Global_ACC, z_Local_ACC, z_Global_EFF, z_Local_EFF, z_Global_Congruent_RT, z_Global_Incongruent_RT, z_Local_Congruent_RT, z_Local_Incongruent_RT, z_Global_Congruent_ACC, z_Global_Incongruent_ACC, z_Local_Congruent_ACC, z_Local_Incongruent_ACC, z_Global_Congruent_EFF, z_Global_Incongruent_EFF, z_Local_Congruent_EFF, z_Local_Incongruent_EFF \n')
for i in range(0, len(z_score_rt_c)):
        f.write('%s, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n'%(int(i)+1, z_score_rt_c[i], z_score_rt_i[i], z_score_acc_c[i], z_score_acc_i[i], z_score_eff_c[i], z_score_eff_i[i], z_score_rt_g[i], z_score_rt_l[i], z_score_acc_g[i], z_score_acc_l[i], z_score_eff_g[i], z_score_eff_l[i], z_score_rt_gc[i], z_score_rt_gi[i], z_score_rt_lc[i], z_score_rt_li[i], z_score_acc_gc[i], z_score_acc_gi[i], z_score_acc_lc[i], z_score_acc_li[i], z_score_eff_gc[i], z_score_eff_gi[i], z_score_eff_lc[i], z_score_eff_li[i]))

# Plot 1. Reaction Times
pyplot.figure(1,figsize=(8,6),dpi=300)
pyplot.subplot(221)
pyplot.title('Distribution of the Reaction Times Means')
pyplot.ylabel('Frequency')
pyplot.xlabel('Reaction time')
pyplot.axis([0.3, 0.9, 0, 8])
pyplot.hist(g_glob_cong_RT_mean, histtype='bar', label='Global Congruent', \
color='blue')
pyplot.subplot(222)
pyplot.ylabel('Frequency')
pyplot.xlabel('Reaction time')
pyplot.axis([0.3, 0.9, 0, 8])
pyplot.hist(g_glob_incong_RT_mean, histtype='bar', label='Global Incongruent', \
color='red')
pyplot.subplot(223)
pyplot.ylabel('Frequency')
pyplot.xlabel('Reaction time')
pyplot.axis([0.3, 0.9, 0, 8])
pyplot.hist(g_loc_cong_RT_mean, histtype='bar', label='Local Congruent', \
color='yellow')
pyplot.subplot(224)
pyplot.ylabel('Frequency')
pyplot.xlabel('Reaction time')
pyplot.axis([0.3, 0.9, 0, 10])
pyplot.hist(g_loc_incong_RT_mean, histtype='bar', label='Local Incongruent', \
color='green')
pyplot.show()

# Plot 2. Accuracy
pyplot.figure(2,figsize=(8,6),dpi=300)
pyplot.subplot(221)
pyplot.title('Distribution of the Accuracy Scores')
pyplot.ylabel('Frequency')
pyplot.xlabel('Accuracy')
pyplot.axis([80, 100, 0, 11])
pyplot.hist(g_glob_cong_ACC_mean, histtype='bar', label='Global Congruent', \
color='blue')
pyplot.subplot(222)
pyplot.ylabel('Frequency')
pyplot.xlabel('Accuracy')
pyplot.axis([50, 100, 0, 11])
pyplot.hist(g_glob_incong_ACC_mean, histtype='bar', label='Global Incongruent', \
color='red')
pyplot.subplot(223)
pyplot.ylabel('Frequency')
pyplot.xlabel('Accuracy')
pyplot.axis([80, 100, 0, 14])
pyplot.hist(g_loc_cong_ACC_mean, histtype='bar', label='Local Congruent', \
color='yellow')
pyplot.subplot(224)
pyplot.ylabel('Frequency')
pyplot.xlabel('Accuracy')
pyplot.axis([80, 100, 0, 11])
pyplot.hist(g_loc_incong_ACC_mean, histtype='bar', label='Local Incongruent', \
color='green')
pyplot.show()

# Plot 3. Efficiency
pyplot.figure(3,figsize=(8,6),dpi=300)
pyplot.subplot(221)
pyplot.title('Distribution of the Efficiency Means')
pyplot.ylabel('Frequency')
pyplot.xlabel('Efficiency')
pyplot.axis([0.3, 0.9, 0, 8])
pyplot.hist(g_glob_cong_EFF_mean, histtype='bar', label='Global Congruent', \
color='blue')
pyplot.subplot(222)
pyplot.ylabel('Frequency')
pyplot.xlabel('Efficiency')
pyplot.axis([0.3, 1, 0, 8])
pyplot.hist(g_glob_incong_EFF_mean, histtype='bar', label='Global Incongruent', \
color='red')
pyplot.subplot(223)
pyplot.ylabel('Frequency')
pyplot.xlabel('Efficiency')
pyplot.axis([0.3, 0.9, 0, 8])
pyplot.hist(g_loc_cong_EFF_mean, histtype='bar', label='Local Congruent', \
color='yellow')
pyplot.subplot(224)
pyplot.ylabel('Frequency')
pyplot.xlabel('Efficiency')
pyplot.axis([0.3, 0.9, 0, 8])
pyplot.hist(g_loc_incong_EFF_mean, histtype='bar', label='Local Incongruent', \
color='green')
pyplot.show()

#==============================================================================
# GROUP MEAN PLOTS
#==============================================================================

###INTERACTION###
#RT
# Get means we need to plot
RT_Means = Group_Global_Congruent_RT, Group_Global_Incongruent_RT, \
Group_Local_Congruent_RT, Group_Local_Incongruent_RT
StdErr_RT = [StdErr_Global_Congruent_RT, StdErr_Global_Incongruent_RT, \
StdErr_Local_Congruent_RT, StdErr_Local_Incongruent_RT]
width = 0.7
ind = [1, 2, 3, 4]
labels = ["GC", "GI", "LC", "LI"]
# Plot them
f4 = pyplot.figure(4,figsize=(8,6),dpi=500)
fig, ax = pyplot.subplots()
f4 = pyplot.ylabel('Reaction Time (seconds)')
f4 = pyplot.xlabel('Condition')
f4 = ax.set_ylim([0.47, 0.62])
f4 = ax.set_xlim([0.5, 4.5])
f4 = pyplot.xticks(ind, labels, rotation='horizontal')
f4 = pyplot.tick_params(axis='both', which='both', bottom='on', top='off',\
                   left='on', right='off')
f4 = pyplot.bar([1, 2, 3, 4], RT_Means, width, color='w', \
align='center', yerr=StdErr_RT)
f4 = pyplot.title('Reaction time group means per condition')
f4 = pyplot.savefig(r'Results\RT.png')
f4 = pyplot.show()

#ACC
# Get means we need to plot
ACC_Means = Group_Global_Congruent_ACC, Group_Global_Incongruent_ACC, \
Group_Local_Congruent_ACC, Group_Local_Incongruent_ACC
StdErr_ACC = [StdErr_Global_Congruent_ACC, StdErr_Global_Incongruent_ACC, \
StdErr_Local_Congruent_ACC, StdErr_Local_Incongruent_ACC]
width = 0.7
ind = [1, 2, 3, 4]
labels = ["GC", "GI", "LC", "LI"]
# Plot them
f5 = pyplot.figure(5,figsize=(8,6),dpi=500)
fig, ax = pyplot.subplots()
f5 = pyplot.ylabel('Accuracy (%)')
f5 = pyplot.xlabel('Condition')
f5 = ax.set_ylim([85, 100])
f5 = ax.set_xlim([0.5, 4.5])
f5 = pyplot.xticks(ind, labels, rotation='horizontal')
f5 = pyplot.tick_params(axis='both', which='both', bottom='on', top='off',\
                   left='on', right='off')
f5 = pyplot.bar([1, 2, 3, 4], ACC_Means, width, color='w', \
align='center', yerr=StdErr_ACC)
f5 = pyplot.title('Accuracy group means per condition')
f5 = pyplot.savefig(r'Results\ACC.png')
f5 = pyplot.show()

#EFF
#RT
# Get means we need to plot
EFF_Means = Group_Global_Congruent_EFF, Group_Global_Incongruent_EFF, \
Group_Local_Congruent_EFF, Group_Local_Incongruent_EFF
StdErr_EFF = [StdErr_Global_Congruent_EFF, StdErr_Global_Incongruent_EFF, \
StdErr_Local_Congruent_EFF, StdErr_Local_Incongruent_EFF]
width = 0.7
ind = [1, 2, 3, 4]
labels = ["GC", "GI", "LC", "LI"]
# Plot them
f6 = pyplot.figure(6,figsize=(8,6),dpi=500)
fig, ax = pyplot.subplots()
f6 = pyplot.ylabel('Efficiency scores')
f6 = pyplot.xlabel('Condition')
f6 = ax.set_ylim([0.5, 0.7])
f6 = ax.set_xlim([0.5, 4.5])
f6 = pyplot.xticks(ind, labels, rotation='horizontal')
f6 = pyplot.tick_params(axis='both', which='both', bottom='on', top='off',\
                   left='on', right='off')
f6 = pyplot.bar([1, 2, 3, 4], EFF_Means, width, color='w', \
align='center', yerr=StdErr_EFF)
f6 = pyplot.title('Efficiency score group means per condition')
f6 = pyplot.savefig(r'Results\EFF.png')
f6 = pyplot.show()



###GLOBAL Vs LOCAL###
#RT
# Get means we need to plot
RT_Means = Group_Global_RT, Group_Local_RT
StdErr_RT = [StdErr_Global_RT, StdErr_Local_RT]
width = 0.7
ind = [1, 2]
labels = ["Global", "Local"]
# Plot them
f7 = pyplot.figure(7,figsize=(8,6),dpi=500)
fig, ax = pyplot.subplots()
f7 = pyplot.ylabel('Reaction Time (seconds)')
f7 = pyplot.xlabel('Condition')
f7 = ax.set_ylim([0.50, 0.6])
f7 = ax.set_xlim([0.5, 2.5])
f7 = pyplot.xticks(ind, labels, rotation='horizontal')
f7 = pyplot.tick_params(axis='both', which='both', bottom='on', top='off',\
                   left='on', right='off')
f7 = pyplot.bar([1, 2], RT_Means, width, color='w', \
align='center', yerr=StdErr_RT)
f7 = pyplot.title('Reaction time group means, Global vs. Local')
f7 = pyplot.savefig(r'Results\GlobVsLocRT.png')
f7 = pyplot.show()

#ACC
# Get means we need to plot
ACC_Means = Group_Global_ACC, Group_Local_ACC
StdErr_ACC = [StdErr_Global_ACC, StdErr_Local_ACC]
width = 0.7
ind = [1, 2]
labels = ["Global", "Local"]
# Plot them
f8 = pyplot.figure(8,figsize=(8,6),dpi=500)
fig, ax = pyplot.subplots()
f8 = pyplot.ylabel('Accuracy (%)')
f8 = pyplot.xlabel('Condition')
f8 = ax.set_ylim([85, 100])
f8 = ax.set_xlim([0.5, 2.5])
f8 = pyplot.xticks(ind, labels, rotation='horizontal')
f8 = pyplot.tick_params(axis='both', which='both', bottom='on', top='off',\
                   left='on', right='off')
f8 = pyplot.bar([1, 2], ACC_Means, width, color='w', \
align='center', yerr=StdErr_ACC)
f8 = pyplot.title('Accuracy group means. Global Vs Local.')
f8 = pyplot.savefig(r'Results\GlobVsLocACC.png')
f8 = pyplot.show()

#EFF
#RT
# Get means we need to plot
EFF_Means = Group_Global_EFF, Group_Local_EFF
StdErr_EFF = [StdErr_Global_EFF, StdErr_Local_EFF]
width = 0.7
ind = [1, 2]
labels = ["Global", "Local"]
# Plot them
f9 = pyplot.figure(9,figsize=(8,6),dpi=500)
fig, ax = pyplot.subplots()
f9= pyplot.ylabel('Efficiency scores')
f9 = pyplot.xlabel('Condition')
f9 = ax.set_ylim([0.5, 0.65])
f9 = ax.set_xlim([0.5, 2.5])
f9 = pyplot.xticks(ind, labels, rotation='horizontal')
f9 = pyplot.tick_params(axis='both', which='both', bottom='on', top='off',\
                   left='on', right='off')
f9 = pyplot.bar([1, 2], EFF_Means, width, color='w', \
align='center', yerr=StdErr_EFF)
f9 = pyplot.title('Efficiency score group means. Global Vs Local.')
f9 = pyplot.savefig(r'Results\GlobVsLocEFF.png')
f9 = pyplot.show()



###CONGRUENT Vs INCONGRUENT ###
#RT
# Get means we need to plot
RT_Means = Group_Congruent_RT, Group_Incongruent_RT
StdErr_RT = [StdErr_Congruent_RT, StdErr_Incongruent_RT]
width = 0.7
ind = [1, 2]
labels = ["Congruent", "Incongruent"]
# Plot them
f10 = pyplot.figure(10,figsize=(8,6),dpi=500)
fig, ax = pyplot.subplots()
f10 = pyplot.ylabel('Reaction Time (seconds)')
f10 = pyplot.xlabel('Condition')
f10 = ax.set_ylim([0.50, 0.6])
f10 = ax.set_xlim([0.5, 2.5])
f10 = pyplot.xticks(ind, labels, rotation='horizontal')
f10 = pyplot.tick_params(axis='both', which='both', bottom='on', top='off',\
                   left='on', right='off')
f10 = pyplot.bar([1, 2], RT_Means, width, color='w', \
align='center', yerr=StdErr_RT)
f10 = pyplot.title('Reaction time group means, Congruent vs. Incongruent')
f10 = pyplot.savefig(r'Results\CongVsIncongRT.png')
f10 = pyplot.show()

#ACC
# Get means we need to plot
ACC_Means = Group_Congruent_ACC, Group_Incongruent_ACC
StdErr_ACC = [StdErr_Congruent_ACC, StdErr_Incongruent_ACC]
width = 0.7
ind = [1, 2]
labels = ["Congruent", "Incongruent"]
# Plot them
f11 = pyplot.figure(11,figsize=(8,6),dpi=500)
fig, ax = pyplot.subplots()
f11 = pyplot.ylabel('Accuracy (%)')
f11 = pyplot.xlabel('Condition')
f11 = ax.set_ylim([85, 100])
f11 = ax.set_xlim([0.5, 2.5])
f11 = pyplot.xticks(ind, labels, rotation='horizontal')
f11 = pyplot.tick_params(axis='both', which='both', bottom='on', top='off',\
                   left='on', right='off')
f11 = pyplot.bar([1, 2], ACC_Means, width, color='w', \
align='center', yerr=StdErr_ACC)
f11 = pyplot.title('Accuracy group means. Congruent vs. Incongruent.')
f11 = pyplot.savefig(r'Results\CongVsIncongACC.png')
f11 = pyplot.show()

#EFF
#RT
# Get means we need to plot
EFF_Means = Group_Congruent_EFF, Group_Incongruent_EFF
StdErr_EFF = [StdErr_Congruent_EFF, StdErr_Incongruent_EFF]
width = 0.7
ind = [1, 2]
labels = ["Congruent", "Incongruent"]
# Plot them
f12 = pyplot.figure(9,figsize=(8,6),dpi=500)
fig, ax = pyplot.subplots()
f12 = pyplot.ylabel('Efficiency scores')
f12 = pyplot.xlabel('Condition')
f12 = ax.set_ylim([0.5, 0.65])
f12 = ax.set_xlim([0.5, 2.5])
f12 = pyplot.xticks(ind, labels, rotation='horizontal')
f12 = pyplot.tick_params(axis='both', which='both', bottom='on', top='off',\
                   left='on', right='off')
f12 = pyplot.bar([1, 2], EFF_Means, width, color='w', \
align='center', yerr=StdErr_EFF)
f12 = pyplot.title('Efficiency score group means. Congruent vs. Incongruent.')
f12 = pyplot.savefig(r'Results\CongVsIncongEFF.png')
f12 = pyplot.show()




print 'T TESTS'
print "congruent vs incongruent RT means p = %.4f"%(RT_CongVsIncong[1])
print "congruent vs incongruent ACC means p = %.4f"%(ACC_CongVsIncong[1])
print "congruent vs incongruent EFF means p = %.4f"%(EFF_CongVsIncong[1])
print ''
print "global vs local RT means p = %.4f"%(RT_GlobVsLoc[1])
print "global vs local ACC means p = %.4f"%(ACC_GlobVsLoc[1])
print "global vs local EFF means p = %.4f"%(EFF_GlobVsLoc[1])
print ''
print 'ANOVA (Interaction)'
print "RT Interaction p = %.4f"%(RT_Interaction[1])
print "ACC Interaction p = %.4f"%(ACC_Interaction[1])
print "EFF Interaction p = %.4f"%(EFF_Interaction[1])



print Group_Global_EFF, Group_Local_EFF
print 'Your analysis is done!'
f.close()