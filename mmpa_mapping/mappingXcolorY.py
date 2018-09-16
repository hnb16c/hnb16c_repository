# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 22:20:15 2018

@author: Hiroshi Nakano
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.figure as figure
import matplotlib.cm as cm
from sklearn.manifold import Isomap
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

#import sklearn as sk
#from sklearn.preprocessing import StandardScaler
from gtm_kaneko import gtm
#from gtm_amaotone import GTM







def loadXY(xyfile):
    df0 = pd.read_csv(xyfile)
    #print(df0)
    df1 = df0.drop('id',axis=1)
    XY = df1.as_matrix()
    #print(XY)
    

    idvec = df0['id'].as_matrix()
    #print(idvec)
    nmol=len(idvec)
    id_to_no={}
    
    for i in range(nmol): id_to_no[idvec[i]]=i
    #print(len(id_to_no))
    # skip descriptors with all zero for all molecules.    ncol = XY.shape[1]
    normcolXY = np.linalg.norm(XY,axis=0)
    #print(len(normcolXY))
    XY = XY[:,normcolXY>0]     
    #print(XY.shape)
    return XY, id_to_no

def loadMMP(mmpfile,id_to_no):
    #mmpvec = []
    df2 = pd.read_csv(mmpfile)
    df2 = df2[df2['flagAC']==1]
    df3=df2.loc[:,['id_a','id_b']]
    #XY = df1.as_matrix()
    #print(XY)
    mmpno=np.zeros(df3.shape,int)
    nmol = df3.shape[0]
    #print(df3)
    #print(df3.iat[0,0])
    for i in range(nmol):
        mmpno[i,0] = id_to_no[df3.iat[i,0]]
        mmpno[i,1] = id_to_no[df3.iat[i,1]]
    
    for i in range(nmol):
        if mmpno[i,0]==3025:
            print(mmpno[i,0])
        if mmpno[i,1]==3025:
            print(mmpno[i,1])

    # skip descriptors with all zero for all molecules.    ncol = XY.shape[1]
    #normcolXY = np.linalg.norm(XY,axis=0)
    #print(len(normcolXY))
    #XY = XY[:,normcolXY>0]     
    #print(XY.shape)
    return mmpno

def scalingXY(X,Y):
    scaler_X = StandardScaler()
    scaler_Y = StandardScaler()
    Xs = scaler_X.fit_transform(X)
    Ys = scaler_Y.fit_transform(np.reshape(Y,(len(Y),1)))
    #Xvs = scaler_X.transform(Xv)
    #Yvs = scaler_Y.transform(np.reshape(Yv,(len(Yv),1)))
    return [Xs,Ys,scaler_X, scaler_Y]

def exec_gtm(X,Y,mmpno):
    # settings
    shape_of_map = [30, 30]
    shape_of_rbf_centers = [12, 12]
    variance_of_rbfs = 0.05
    lambda_in_em_algorithm = 0#1.0e-5
    number_of_iterations = 100
    display_flag = 1
    
    model = gtm(shape_of_map, shape_of_rbf_centers, variance_of_rbfs, lambda_in_em_algorithm, number_of_iterations,
            display_flag)
    model.fit(X)

    if model.success_flag:
        # calculate of responsibilities
        responsibilities = model.responsibility(X)

        Ymax = np.max(Y)
        Ymin = np.min(Y)
         
        Y0to1 = (Y-Ymin)/(Ymax-Ymin)
        # plot the mean of responsibilities
        means = responsibilities.dot(model.map_grids)
        plt.rcParams["font.size"] = 18
        plt.rcParams["font.family"] = "Serif"
        plt.figure(figsize=figure.figaspect(1))
        plt.figure(figsize=(10,10))
        plt.scatter(means[:, 0], means[:, 1], c=cm.RdYlGn(1-Y0to1),s=100)
        #plt.ylim(-1.1, 1.1)
        #plt.xlim(-1.1, 1.1)
        plt.xlabel("z1 (mean)")
        plt.ylabel("z2 (mean)")
        plt.show()

        # plot the mean of responsibilities with AC
        means = responsibilities.dot(model.map_grids)
        plt.rcParams["font.size"] = 18
        plt.rcParams["font.family"] = "Serif"        
        plt.figure(figsize=figure.figaspect(1))
        plt.figure(figsize=(10,10))
        plt.scatter(means[:, 0], means[:, 1], c=cm.RdYlGn(1-Y0to1),s=100)
        print(means.shape)
        for i,no_a in enumerate(mmpno[:,0]):
            no_b = mmpno[i,1]
            if no_b >= 3025:
                print(i)
            plt.plot([means[no_a,0],means[no_b,0]],[means[no_a,1],means[no_b,1]],color='blue')
        #plt.ylim(-1.1, 1.1)
        #plt.xlim(-1.1, 1.1)
        plt.xlabel("z1 (mean)")
        plt.ylabel("z2 (mean)")
        plt.show()
        
        # plot the mode of responsibilities
        modes = model.map_grids[responsibilities.argmax(axis=1), :]
        plt.figure(figsize=figure.figaspect(1))
        plt.figure(figsize=(10,10))
        plt.scatter(modes[:, 0], modes[:, 1], c=cm.RdYlGn(1-Y0to1),s=100)
        #plt.ylim(-1.1, 1.1)
        #plt.xlim(-1.1, 1.1)
        plt.xlabel("z1 (mode)")
        plt.ylabel("z2 (mode)")
        plt.show()

        
        # plot the mode of responsibilities
        modes = model.map_grids[responsibilities.argmax(axis=1), :]
        plt.figure(figsize=figure.figaspect(1))
        plt.figure(figsize=(10,10))
        plt.scatter(modes[:, 0], modes[:, 1], c=cm.RdYlGn(1-Y0to1),s=100)
        for i,no_a in enumerate(mmpno[:,0]):
            no_b = mmpno[i,1]
            if no_b >= 3025:
                print(i)
            plt.plot([modes[no_a,0],modes[no_b,0]],[modes[no_a,1],modes[no_b,1]],color='blue')

        #plt.ylim(-1.1, 1.1)
        #plt.xlim(-1.1, 1.1)
        plt.xlabel("z1 (mode)")
        plt.ylabel("z2 (mode)")
        plt.show()

def exec_isomap(X,Y,mmpno):
    #n_neighbors=20
    isomap = Isomap(n_neighbors=10, n_components=2,eigen_solver='dense')
    X_iso = isomap.fit(X).transform(X)
    Ymax = np.max(Y)
    Ymin = np.min(Y)
    Y0to1 = (Y-Ymin)/(Ymax-Ymin)
    

    plt.figure(figsize=(1,8))
    #plt.scatter(, X_iso[:, 1], c=cm.RdYlGn(1-y),s=30)
    plt.show()
    
    plt.figure(figsize=(8,8))
    plt.rcParams["font.size"] = 18
    plt.rcParams["font.family"] = "Serif"
    plt.scatter(X_iso[:, 0], X_iso[:, 1], c=cm.RdYlGn(1-Y0to1),s=30)
    plt.ylim(-20, 20)
    plt.xlim(-20, 20)
    plt.xlabel("z1")
    plt.ylabel("z2")
    plt.show()
    
    
    #plt.figure(figsize=figure.figaspect(1))
    plt.figure(figsize=(8,8))
    plt.scatter(X_iso[:, 0], X_iso[:, 1], c=cm.RdYlGn(1-Y0to1),s=30)
    for i,no_a in enumerate(mmpno[:,0]):
        no_b = mmpno[i,1]
        if no_b >= 3025:
            print(i)
        plt.plot([X_iso[no_a,0],X_iso[no_b,0]],[X_iso[no_a,1],X_iso[no_b,1]],color='blue')

    
    
    plt.ylim(-20, 20)
    plt.xlim(-20, 20)
    plt.xlabel("z1")
    plt.ylabel("z2")
    plt.show()

def colorbar_show(Y):
    # colorbar
    Ymax = np.max(Y)
    Ymin = np.min(Y)
    
    fig,ax = plt.subplots(figsize=(0.5,8))
    gradient = np.linspace(0, 1, 256)
    x05vec=np.ones((256),float)*0.5
    ax.tick_params(labelbottom=False,bottom=False)
    #gradient = np.vstack((gradient, gradient))
    #print(gradient)
    #print(x05vec)
    #ax.set_ylabel('p$K_i$')
    ax.set_title('p$K_i$')
    ax.set_ylim(Ymin, Ymax)
    ax.set_xlim(0, 1)
    #ax.figure(figsize=(1,8))
    ax.scatter(x05vec, gradient*(Ymax-Ymin)+Ymin, c=cm.RdYlGn(1-gradient),s=1000)    
    #plt.show()


if __name__ == '__main__':
    xyfile="st11sample_dataset.csv"
    XY, id_to_no=loadXY(xyfile)
    mmpno=loadMMP("st12sample_dataset.csv",id_to_no)
    Y=XY[:,-1]
    colorbar_show(Y)
    [Xs,Ys,scaler_X, scaler_Y] = scalingXY(XY[:,0:-1],Y)
    pca_model = PCA(n_components=6)
    Xpca=pca_model.fit_transform(Xs)
    print('各次元の寄与率: {0}'.format(pca_model.explained_variance_ratio_))
    print('累積寄与率: {0}'.format(sum(pca_model.explained_variance_ratio_)))

    #print(Xs)
    #exec_gtm(Xs,Y,mmpno)
    #exec_isomap(Xs,Y,mmpno)
    exec_gtm(Xpca,Y,mmpno)




