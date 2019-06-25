# -*- coding: utf-8 -*-
#%% 
import numpy as np
import matplotlib.pylab as plt
import dclab
import pandas as pd
from matplotlib import cm
from scipy.stats import gaussian_kde

######## colormap chan be changed here
cmap_vir = cm.get_cmap('viridis')


def density_scatter( x , y, bins, ax, sort = True, **kwargs )   :

    np.nan_to_num(y, copy=False)
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    
    ax.scatter( x, y, c=z, cmap = cmap_vir, marker = ".", s = 4, picker = True, **kwargs )
    plt.subplots_adjust(wspace = 0.2)
    plt.minorticks_on()
    #plt.grid(b=True, which='both', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.85', linestyle='--')
    plt.grid(b=True, which='major', color='0.85', linestyle='-')
    plt.xticks()
    plt.yticks()
    plt.rcParams["font.size"] = 12
    plt.gcf().set_tight_layout(False)

    return ax


def ds_filter_results (celltype, filepath, ds_child2, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list) :
    
    area_mean = np.mean(ds_child2["area_um"])
    area_median = np.median(ds_child2["area_um"])
    area_std = np.std(ds_child2["area_um"])
    defo_mean = np.mean(ds_child2["deform"])
    defo_median = np.median(ds_child2["deform"])
    defo_std = np.std(ds_child2["deform"])
    bright_mean = np.mean(ds_child2["bright_avg"])
    bright_std = np.std(ds_child2["bright_avg"])
    emodulus_mean = np.nanmean(ds_child2["emodulus"])
    emodulus_std = np.nanstd(ds_child2["emodulus"])
    temperature = ds["temp"][0]
    viscosity = (-3.14318907716840 * temperature) + 136.42041546595024 # in mPa.s
        #emodulus = ds["emodulus"]
    counts_gated = len(ds_child2)
    counts_pct = (100*len(ds_child2)/counts[1])
    
    dict1 = {"Cell_subset":celltype, "Counts": counts[1], "Counts_gated": counts_gated, "Counts_percentage": counts_pct,
             "Area":area_mean, "Area_SD":area_std, "Deformation":defo_mean, "Deformation_SD":defo_std, 
             "Youngs_Modulus":emodulus_mean, "Youngs_Modulus_SD":emodulus_std, "Brightness":bright_mean, "Brightness_SD":bright_std}
    #dict1.update(blah..) 
    rows_list.append(dict1)  # in loop

    def onpick(event):
        ind = event.ind
        print('onpick scatter event number:', ind)
        print('Shown index', ind[0])
        print('length of index', len(ind))
        print('area of event', ds_child2["area_um"][ind[0]])

        plt.figure(figsize=(10,5))
        ax1 = plt.subplot(211, title="image")
        ax2 = plt.subplot(212, title="mask")
        ax1.imshow(ds_child2["image"][ind[0]], cmap="gray")
        ax2.imshow(ds_child2["mask"][ind[0]])

        print(ds_child2["trace"][ind[0]])
#        plt.figure(figsize=(10,5))
#        ax1 = plt.subplot(211, title="trace 1")
#        ax2 = plt.subplot(212, title="trace 2")
#        

        
    figure = plt.figure(figsize=(20,10))

    ax = plt.subplot(231, xlabel = FL1, xlim = (1,10000), xscale = "log", ylabel = FL2, ylim = (0.1, 10000), yscale = "log")
    density_scatter(ds_child2["fl1_max"], ds_child2["fl2_max"], bins = [1000,100], ax = ax)
    ax = plt.subplot(232, xlabel = FL1, xlim = (1,10000), xscale = "log", ylabel = FL3, ylim = (0.1, 10000), yscale = "log")
    density_scatter(ds_child2["fl1_max"], ds_child2["fl3_max"], bins = [1000,100], ax = ax)
    ax = plt.subplot(233, xlabel = FL2, xlim = (1,10000), xscale = "log", ylabel = FL3, ylim = (0.1, 10000), yscale = "log")
    density_scatter(ds_child2["fl2_max"], ds_child2["fl3_max"], bins = [1000,100], ax = ax)

    ax = plt.subplot(234, xlabel = 'Area [um]', xlim = (10,130), ylabel = 'Deformation [a.u.]', ylim = (0, 0.4))
    density_scatter(ds_child2["area_um"], ds_child2["deform"], bins = [1000,100], ax = ax)
    ax = plt.subplot(235, xlabel = 'Area [um]', xlim = (10,130), ylabel = 'Brightness [a.u.]', ylim = (90, 160))
    density_scatter(ds_child2["area_um"], ds_child2["bright_avg"], bins = [1000,100], ax = ax)
    ax = plt.subplot(236, xlabel = 'Area [um]', xlim = (10,130), ylabel = 'Youngs modulus [kPa]', ylim = (3.5, 12))
    density_scatter(ds_child2["area_um"], ds_child2["emodulus"], bins = [1000,100], ax = ax)

    plt.text(0.02, 0.8, 'Patient: ' + patient + '\n' + exp_date + '\n' + exp_time +
             '\n \n Filters: \n' + celltype, fontsize=14, transform=plt.gcf().transFigure)
    plt.text(0.02, 0.6, 'Temperature: \n %.2f \nViscosity: \n %.2f \n' %(temperature, viscosity),
             fontsize=12, transform=plt.gcf().transFigure)
    plt.text(0.02, 0.1, 'Area: \n Mean = \n  %.3f \n Std = \n %.3f \n \
             \n Deformation: \n Mean = \n  %.3f \n Std = \n %.3f \n \
             \n E modulus: \n Mean = \n  %.3f \n Std = \n %.3f \
             \n \n \n Counts: \n Total = \n %.2f \n Gated = \n %.2f \n Perc. = \n %.2f'
             %(area_mean, area_std, defo_mean, defo_std, emodulus_mean, emodulus_std, counts[1], counts_gated, counts_pct),
             fontsize=12, transform=plt.gcf().transFigure) 
    # str(area_mean)
    
    #show maximised window
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()
    
    figure.canvas.mpl_connect('pick_event', onpick)

    plt.savefig(filepath + "_" + celltype , bbox_inches='tight')
    
    return (rows_list)

 
#%%  
# for each new patient click here and then press 'Ctrl + Enter'
    # this will open 12 figures with graphs, it will also save them as .pngs in the specified folder
    # if you later want to close all figures at once, type this in Spyder command line: 
    # plt.close('all')    

###################################################################################################################
############################################## LOAD DATA ###########################################################
###################################################################################################################

# HERE SPECIFY FILE PATH, PATIENT INFO, PANEL 1 OR 2, THRESHOLDS OF FLUORESCENCE
filepath = r"E:\UK ERLANGEN DATA\20190521_Marketa_Alex_PBMC\A1 good\P2\0p06\M001_data"
ds = dclab.new_dataset(filepath + ".rtdc")
patient = 'A1 \n spondyloarthritis' #here specify patient name, disease, whatever (\n is new line)
panel = 2 #define if panel 1 or panel 2
FL1_threshold = 250
FL2_threshold = 200
FL3_threshold = 100
#print(ds.features)

# plt.close('all')  #### to close all figures

########################################## IMAGES #######################################################
#
## store images and masks in arrays
#image=[]
#mask =[]
#for ii in range(len(ds)):
#    image.append(ds["image"][ii])
#    mask.append(ds["mask"][ii])
#    ## this is equivalent to ds["bright_avg"][ii]
#    #bright_avg = np.mean(image[mask])
#    #print("average brightness of event {}: {:.1f}".format(ii, bright_avg))
#
## show a particular image and mask
##ax1 = plt.subplot(211, title="image")
##ax2 = plt.subplot(212, title="mask")
##ax1.imshow(ds["image"][10], cmap="gray")
##ax2.imshow(ds["mask"][10])

############################### BASIC FILTERS (first plot = fig 1) #############################################


exp_date = ds.config["experiment"]["date"]
exp_time = ds.config["experiment"]["time"]
ds.config["calculation"]["emodulus temperature"] = ds["temp"]
ds.config["calculation"]["emodulus viscosity"] = 54.7 # 0.8 MC
ds.config["filtering"]["area_um min"] = 25
ds.config["filtering"]["area_um max"] = 600
ds.config["filtering"]["area_ratio min"] = 1
ds.config["filtering"]["area_ratio max"] = 1.05
ds.config["filtering"]["bright_avg min"] = 115   # to exclude too dark / too bright events
ds.config["filtering"]["bright_avg max"] = 300
ds.config["filtering"]["y_pos min"] = 13   # the y_pos limits have to be changed if the channel was off-center
ds.config["filtering"]["y_pos max"] = 16
ds.apply_filter()  # this step is important!
#area_um_filtered = ds["area_um"][ds.filter.all] #filtered parameter
ds_child = dclab.new_dataset(ds) #dataset filtered according to the set filter
ds_child.apply_filter() ### will change ds_child when the filters for ds are changed
counts = (len(ds), len(ds_child))
elasticity = dclab.features.emodulus.get_emodulus(ds_child["area_um"], ds_child["deform"], medium = 54.7, channel_width = 20.0, 
                                     flow_rate = ds_child.config["setup"]["flow rate"], px_um =0.34, temperature = ds_child["temp"][0], copy = False )
ds_child_dict = {"emodulus":elasticity, "area_um":ds_child["area_um"], "deform":ds_child["deform"],
                 "area_ratio":ds_child["area_ratio"], "bright_avg":ds_child["bright_avg"], "bright_sd":ds_child["bright_sd"], "image":ds_child["image"], 
                 "mask":ds_child["mask"], "fl1_max":ds_child["fl1_max"], "fl2_max":ds_child["fl2_max"], "fl3_max":ds_child["fl3_max"]}
ds_child = dclab.rtdc_dataset.RTDC_Dict(ds_child_dict)
print(ds_child.features)
#"emodulus" in ds_child #test feature availability
#ds_child.config["emodulus"] = elasticity


emodulus_mean = np.mean(ds_child["emodulus"])
emodulus_std = np.std(ds_child["emodulus"])



#"image" in ds ##### test feature availability
#plt.imshow(ds["image"][10])
#plt.figure()  
#density_scatter(ds_child["pos_x"], ds_child["pos_y"], bins = [1000,100], ax =None)


################################## ADVANCED FILTERS (fig 2 to fig 12) ############################################

rows_list = [] # na zacatku
del rows_list
rows_list = [] 

if panel == 1 :
    FL1 = 'CD4'
    FL2 = 'CD8'
    FL3 = 'CD3'
    rows_list = ds_filter_results ('BASIC', filepath, ds_child, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    #plt.savefig(filepath)

    ds_child.config["filtering"]["area_um min"] = 30
    ds_child.config["filtering"]["area_um max"] = 120
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = -100
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild01 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild01.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('Area_30_120', filepath, ds_grandchild01, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)

    ds_child.config["filtering"]["area_um min"] = 30
    ds_child.config["filtering"]["area_um max"] = 50
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = -100
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild02 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild02.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('Area_30_50', filepath, ds_grandchild02, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["area_um min"] = 50
    ds_child.config["filtering"]["area_um max"] = 120
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = -100
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild03 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild03.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('Area 50_120', filepath, ds_grandchild03, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)

    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = -100
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = FL3_threshold
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild1 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild1.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD3p', filepath, ds_grandchild1, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["area_um min"] = 30
    ds_child.config["filtering"]["area_um max"] = 50
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = -100
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = FL3_threshold
    ds_grandchild2 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild2.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD3n_Area_30_50', filepath, ds_grandchild2, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = FL2_threshold
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = FL3_threshold
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild3 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild3.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD3p_CD8p', filepath, ds_grandchild3, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["fl1_max min"] = FL1_threshold
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = -100
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = FL3_threshold
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild4 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild4.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD3p_CD4p', filepath, ds_grandchild4, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = FL1_threshold
    ds_child.config["filtering"]["fl2_max min"] = -100
    ds_child.config["filtering"]["fl2_max max"] = FL2_threshold
    ds_child.config["filtering"]["fl3_max min"] = FL3_threshold
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild5 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild5.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD3p_CD4n_CD8n', filepath, ds_grandchild5, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["fl1_max min"] = FL1_threshold
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = FL2_threshold
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = FL3_threshold
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild6 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild6.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD3p_CD4p_CD8p', filepath, ds_grandchild6, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["fl1_max min"] = FL1_threshold
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = -100
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild7 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild7.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD4p', filepath, ds_grandchild7, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = FL2_threshold
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild8 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild8.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD8p', filepath, ds_grandchild8, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
else:
    FL1 = 'CD16'
    FL2 = 'CD14'
    FL3 = 'CD19'
    ds_filter_results ('Filters: \n BASIC', filepath, ds_child, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)

    ds_child.config["filtering"]["area_um min"] = 30
    ds_child.config["filtering"]["area_um max"] = 120
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = -100
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild01 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild01.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('Area_30_120', filepath, ds_grandchild01, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["area_um min"] = 30
    ds_child.config["filtering"]["area_um max"] = 50
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = -100
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild02 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild02.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('Area_30_50', filepath, ds_grandchild02, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["area_um min"] = 50
    ds_child.config["filtering"]["area_um max"] = 120
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = -100
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild03 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild03.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('Area_50_120', filepath, ds_grandchild03, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = -100
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = FL3_threshold
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild1 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild1.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD19p', filepath, ds_grandchild1, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["area_um min"] = 30
    ds_child.config["filtering"]["area_um max"] = 50
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = -100
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = FL3_threshold
    ds_grandchild2 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild2.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD19n_Area_30_50', filepath, ds_grandchild2, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = FL2_threshold
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild3 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild3.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD14p', filepath, ds_grandchild3, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["area_um min"] = 30
    ds_child.config["filtering"]["area_um max"] = 50
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = FL2_threshold
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild4 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild4.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD14p_Area_30_50', filepath, ds_grandchild4, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["area_um min"] = 50
    ds_child.config["filtering"]["area_um max"] = 120
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = FL2_threshold
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild5 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild5.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD14p_Area_50_120', filepath, ds_grandchild5, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["fl1_max min"] = -100
    ds_child.config["filtering"]["fl1_max max"] = FL1_threshold
    ds_child.config["filtering"]["fl2_max min"] = FL2_threshold
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild6 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild6.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD14p_CD16n', filepath, ds_grandchild6, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["fl1_max min"] = FL1_threshold
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = FL2_threshold
    ds_child.config["filtering"]["fl2_max max"] = 20000
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild7 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild7.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD14p_CD16p', filepath, ds_grandchild7, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    
    ds_child.config["filtering"]["area_um min"] = 50
    ds_child.config["filtering"]["area_um max"] = 120
    ds_child.config["filtering"]["fl1_max min"] = FL1_threshold
    ds_child.config["filtering"]["fl1_max max"] = 20000
    ds_child.config["filtering"]["fl2_max min"] = -100
    ds_child.config["filtering"]["fl2_max max"] = FL2_threshold
    ds_child.config["filtering"]["fl3_max min"] = -100
    ds_child.config["filtering"]["fl3_max max"] = 20000
    ds_grandchild8 = dclab.new_dataset(ds_child) #dataset filtered according to the set filter
    ds_grandchild8.apply_filter() ### will change ds_child when the filters for ds are changed
    ds_filter_results ('CD14p_CD16p_Area_50_120', filepath, ds_grandchild8, patient, exp_date, exp_time, counts, FL1, FL2, FL3, rows_list)
    

############################################### SAVE TO EXCEL FILE #################################################
df = pd.DataFrame(rows_list, columns=['Cell_subset', 'Counts', 'Counts_gated', 'Counts_percentage', 'Area', 'Area_SD', 'Deformation', 'Deformation_SD', 'Youngs_Modulus', 'Youngs_Modulus_SD', 'Brightness', 'Brightness_SD'])
df.to_csv(filepath + ".csv")
#df.to_csv(r'Path where you want to store the exported CSV file\File Name.csv')