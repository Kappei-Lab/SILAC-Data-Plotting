import pandas as pd
import numpy as np
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import skew
import json    
import os
import math
import ast
import argparse


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
pd.options.display.max_rows = 200
pd.options.display.max_columns = 200

def create_plot_dict(plot_info_path):
    
    """ Create a dictionary from text file for plottting figures
    Dictionary key will be a variable and dictionary value will be value for that variable
    
    Args:
    plot_info_path : location of the plot_info file
    Returns:
    plot_dict: cleaned dictionary
    """
    

    plot_dict={} # Create an empty dictionary
    for line in open(plot_info_path,'r'): # read the plot info line by line
        key,value=line.split(':') # split the line by ':' part before ':' becomes the key and part after becomes the value 
        value=value.strip('\n')     # remove the new line character from the end of value
        plot_dict[key]=value        # create an entry in the dictionary

    # All the value are now strings. We need to convert into correct format    
        
    plot_dict['left_sig']=ast.literal_eval(plot_dict['left_sig']) # Convert string to boolean
    plot_dict['right_sig']=ast.literal_eval(plot_dict['right_sig']) # Convert strint to boolean
    plot_dict['cut_off']=float(plot_dict['cut_off']) # convert string to float
    plot_dict['legend_size']=float(plot_dict['legend_size']) # convert string to float
    plot_dict['title_size']=float(plot_dict['title_size']) # convert string to float
    plot_dict['xylabel_size']=float(plot_dict['xylabel_size']) # convert string to float
    plot_dict['bait_protein']=plot_dict['bait_protein'].split(';') # convert string to list of protein names
    
    # Format lower limit and upper limit values
    
    if (plot_dict['x_lower_limit']!='') and (plot_dict['x_upper_limit']!='') and (plot_dict['y_lower_limit']!='') and (plot_dict['y_upper_limit']!=''):
        
        plot_dict['x_lower_limit']=float(plot_dict['x_lower_limit'])
        plot_dict['x_upper_limit']=float(plot_dict['x_upper_limit'])
        plot_dict['y_lower_limit']=float(plot_dict['y_lower_limit'])
        plot_dict['y_upper_limit']=float(plot_dict['y_upper_limit'])
        
    else:
        plot_dict['x_lower_limit']=None
        plot_dict['x_upper_limit']=None
        plot_dict['y_lower_limit']=None
        plot_dict['y_upper_limit']=None
    
    
    
    return plot_dict


def read_protein_groups(proteingroups_path,directory):
    
    """ Read the protein group file and create a directory for results
    Args:
    proteingroups_path: location of the protein groups file
    Returns:
    pgroups: dataframe containining the proteins groups data
    """

    # Read the file
    pgroups=pd.read_table(proteingroups_path,sep="\t")
    
    # create a folder name for the results folder
    folder_no=1   
    new_directory=directory+"Analysis_"+str(folder_no)
    
    # Check whether the folder already exists. If yes,increase suffix and check again
    
    while os.path.isdir(new_directory):
        folder_no+=1
        new_directory=directory+"Analysis_"+str(folder_no)
    
    
    
    # create a new folder with the finalised name
    os.mkdir(new_directory)

    print("Dimensions of the table are ",pgroups.shape)

    # Get the name of the directory where protein group files are located. This will be used in excel and image 
    # file names
    ls=directory.split("/")
    folder_name=ls[len(ls)-2]
    print("Folder name is",folder_name)
    print("Columns names are :",pgroups.columns.values)
    return pgroups,folder_name,new_directory

def remove_unreliable_proteins(pgroups):
    """
    This function removes the rows containing '+' in any of the columns: 'Only identified by site',
    'Reverse','Potential contaminant'
    Args:
    pgroups: dataframe to be processed
    Returns:
    pgroup2: cleaned dataframe
    
    """
    pgroups2=pgroups.copy()
    print("Initially dimensions of the table are ",pgroups2.shape)
    
    if 'Only identified by site' in pgroups2.columns.values: # Check
        if len(pgroups2['Only identified by site'])!=(pgroups2['Only identified by site'].isnull().sum()):
            pgroups2=pgroups2[~(pgroups2['Only identified by site']=="+")]
    if 'Reverse' in pgroups2.columns.values:
        if len(pgroups2['Reverse'])!=(pgroups2['Reverse'].isnull().sum()):
            pgroups2=pgroups2[~(pgroups2['Reverse']=="+")]
    if 'Potential contaminant' in pgroups2.columns.values:
        if len(pgroups2['Potential contaminant'])!=(pgroups2['Potential contaminant'].isnull().sum()):
            pgroups2=pgroups2[~(pgroups2['Potential contaminant']=="+")]

    print("\nAfter cleaning dimensions of the table are ",pgroups2.shape)
    
    return pgroups2

def filter_for_rev(pgroups2):
    
    print("Fitering rows ....\n")
    
    print("Dimensions of the table before this step are :",pgroups2.shape,"\n")
    
    
    if 'Ratio H/L normalized for'  in pgroups2.columns.values:
        
        forward=pgroups2['Ratio H/L normalized for']
        print('Number of zeros in the column "Ratio H/L normalized for"',sum(forward==0))
        print('Number of null values in the column "Ratio H/L normalized for"',sum(forward.isnull()))
    else:
        print('Column not found')
    if  'Ratio H/L normalized rev'  in pgroups2.columns.values:
        
        reverse=pgroups2['Ratio H/L normalized rev']
        print('Number of zeros in the column "Ratio H/L normalized rev"',sum(reverse==0))
        print('Number of null values in the column "Ratio H/L normalized rev"',sum(reverse.isnull()))
    else:    
        print("Column not found")
        
    rows_to_remove=(pgroups2['Ratio H/L normalized for']==0)|(pgroups2['Ratio H/L normalized for'].isnull())|(pgroups2['Ratio H/L normalized rev']==0) |(pgroups2['Ratio H/L normalized rev'].isnull())
    pgroups3=pgroups2[~rows_to_remove]
    
    print("\nDimensions of the table after this step are :",pgroups3.shape)
    
    return pgroups3

def select_cols(pgroups3,cols):
    
    """ Select important columns from the dataframe
    Args:
    pgroups3: input dataframe
    cols: Columns to be selected
    Returns:
    pgroups4: Dataframe with selected columns
    
    """
    cols_present=[]
    for item in cols:
        if item in pgroups3.columns.values:
            cols_present.append(item)
        else:
            print("{} not found".format(item))


    pgroups4=pgroups3[cols_present]
    
    return pgroups4

def create_plot_df(pgroups4):
    
    """ Create dataframe for plotting with log2 transformed values
    Args:
    pgroups4: final cleaned dataframe
    
    Returns:
    plot_df: a dataframe with three columns: Gene names,log2(forward Ratio),log2(reverse ratio)
    
    """
    
    log2_for=np.log2(pgroups4['Ratio H/L normalized for'])
    log2_rev=np.log2(pgroups4['Ratio H/L normalized rev'])
    plot_df=pd.DataFrame({'Gene_names':pgroups4['Gene names'],
                         'log2_for':log2_for,'log2_rev':log2_rev},columns=['Gene_names','log2_for','log2_rev'])
    
    return plot_df

def create_plot_df_nNorm(pgroups4):
    
    """ Create dataframe for plotting with log2 transformed values
    Args:
    pgroups4: final cleaned dataframe
    
    Returns:
    plot_df: a dataframe with three columns: Gene names,log2(forward Ratio),log2(reverse ratio)
    
    """
    
    log2_for=np.log2(pgroups4['Ratio H/L for'])
    log2_rev=np.log2(pgroups4['Ratio H/L rev'])
    plot_df=pd.DataFrame({'Gene_names':pgroups4['Gene names'],
                         'log2_for':log2_for,'log2_rev':log2_rev},columns=['Gene_names','log2_for','log2_rev'])
    
    return plot_df

def create_legend_df(plot_df,plot_dict):
    
    """ Create the information needed for writing the legend from plot_df file
    
    Args:
    plot_df: dataframe with three columns: Gene_names, log2_for and log2_rev
    Returns:
    legend_df: dataframe containing information for writing the legends
    
    """
    
    def clean_protein_names(col):
        newnames=[]
        for item in col:
            if type(item) is str:
                x=item.split(';',1)[0]
                newnames.append(x)
            else:
                newnames.append('None')
        return newnames

    legend_df=pd.DataFrame(columns=['protein','legend','color'])
    
    if plot_dict['left_sig']:
        color=plot_dict['color_1']
        global up_df
        up_df=plot_df[(plot_df['log2_for']>plot_dict['cut_off']) & (plot_df['log2_rev']< -plot_dict['cut_off'])].sort_values(by='log2_for',                                                                                               ascending=False)
        up_df['rank']=np.arange(1,up_df.shape[0]+1,1)
        up_df['Gene_names']=clean_protein_names(up_df['Gene_names'])
        
        for i in range(up_df.shape[0]):
            legend=str(up_df['rank'].iloc[i])+"."+up_df['Gene_names'].iloc[i]
            legend_df=legend_df.append({'protein':up_df['Gene_names'].iloc[i],'legend':legend,'color':color},ignore_index=True)
    
    if plot_dict['right_sig']:
        if plot_dict['left_sig']:
            color=plot_dict['color_2']
        else:
            color=plot_dict['color_1']
        global down_df    
        down_df=plot_df[(plot_df['log2_for']<-plot_dict['cut_off']) & (plot_df['log2_rev']> plot_dict['cut_off'])].sort_values(by='log2_rev',ascending=False)
        down_df['rank']=np.arange(1,down_df.shape[0]+1,1)    
        down_df['Gene_names']=clean_protein_names(down_df['Gene_names'])
   
        for i in range(down_df.shape[0]):
            legend=str(down_df['rank'].iloc[i])+"."+down_df['Gene_names'].iloc[i]
            legend_df=legend_df.append({'protein':down_df['Gene_names'].iloc[i],'legend':legend,'color':color},ignore_index=True)


    for gene in plot_dict['bait_protein']:
        if sum(legend_df['protein']==gene)>0:
            legend_df['color'][legend_df['protein']==gene]=plot_dict['color_bait']
        else:
            legend_df=legend_df.append({'protein':gene,'legend':gene,'color':plot_dict['color_bait']},ignore_index=True)       

    row=[]
    col=[]
    for i in range(legend_df.shape[0]):
        row.append(i%33)
        col.append(i//33)
    legend_df['row']=row
    legend_df['col']=col

    return legend_df

def save_fig(plot_df,plot_dict,legend_df,new_directory,folder_name):
    
    """ Create and save the figures
    Args:
    plot_df: dataframe with three columns: Gene_names,log2_for,log2_rev
    plot_dict: Dictionary containing values for different parameters
    legend_df: dataframe containing information for plotting the legends
    new_directory: directory to save the results
    
    Returns:
    None
    """

    fig=plt.figure(figsize=(10,10))
    ax=fig.add_subplot(1,1,1)

    up_bool=(plot_df['log2_for']>plot_dict['cut_off']) & (plot_df['log2_rev']< -plot_dict['cut_off'])
    down_bool=(plot_df['log2_for']<-plot_dict['cut_off']) & (plot_df['log2_rev']> plot_dict['cut_off'])
    gene_names = plot_df['Gene_names']
    all_gene_names= plot_df['Gene_names'].str.split(';', expand=True)
    all_gene_names_bool = all_gene_names.isin(plot_dict['bait_protein'])
    bait_bool=all_gene_names_bool.any(axis=1)


    if plot_dict['left_sig'] and plot_dict['right_sig']:
        sel_1=up_bool
        sel_2=down_bool
        sel_3=bait_bool
        sel_4=~(sel_1|sel_2|bait_bool)
        color_1,color_2,color_3=plot_dict['color_1'],plot_dict['color_2'],plot_dict['color_bait']

    elif plot_dict['left_sig']:
        sel_1=up_bool
        sel_2=[False]*plot_df.shape[0]
        sel_3=bait_bool
        sel_4=~(sel_1|sel_2|bait_bool) 
        color_1,color_2,color_3=plot_dict['color_1'],plot_dict['color_1'],plot_dict['color_bait']

    elif plot_dict['right_sig']:
        sel_1=[False]*plot_df.shape[0]
        sel_2=down_bool
        sel_3=bait_bool
        sel_4=~(sel_1|sel_2|bait_bool)
        color_1,color_2,color_3=plot_dict['color_1'],plot_dict['color_1'],plot_dict['color_bait']
    else:
        sel_1=[False]*plot_df.shape[0]
        sel_2=[False]*plot_df.shape[0]
        sel_3=bait_bool
        sel_4=[True]*plot_df.shape[0]
        color_1,color_2,color_3=plot_dict['color_1'],plot_dict['color_1'],plot_dict['color_bait']


    ax.scatter(plot_df['log2_rev'][sel_1],plot_df['log2_for'][sel_1],edgecolors="black",
               color=color_1,s=250,alpha=1,zorder=2)
    ax.scatter(plot_df['log2_rev'][sel_2],plot_df['log2_for'][sel_2],edgecolors="black",
           color=color_2,s=250,alpha=1,zorder=2)

    ax.scatter(plot_df['log2_rev'][sel_3],plot_df['log2_for'][sel_3],edgecolors="black",
               color=color_3,s=250,alpha=1,zorder=2)
    ax.scatter(plot_df['log2_rev'][sel_4],plot_df['log2_for'][sel_4],edgecolors="black",
           color=color_1,s=50,alpha=.3,zorder=1)



    ### number the points #####

    if plot_dict['left_sig']:
        for i in range(0,len(up_df)):
            ax.annotate(up_df['rank'].iloc[i],xy=(up_df['log2_rev'].iloc[i],up_df['log2_for'].iloc[i]),horizontalalignment='center',
                    verticalalignment='center',fontsize=8)
    if plot_dict['right_sig']:
        for i in range(0,len(down_df)):
            ax.annotate(down_df['rank'].iloc[i],xy=(down_df['log2_rev'].iloc[i],down_df['log2_for'].iloc[i]),horizontalalignment='center',
                    verticalalignment='center',fontsize=8)
    ### 

    ### Plot the legends ###

    for i in range(legend_df.shape[0]):
        ax.text(1.01+(.20*legend_df['col'].iloc[i]),1.0-(.03*legend_df['row'].iloc[i]),legend_df['legend'].iloc[i],
                fontsize=plot_dict['legend_size'],color=legend_df['color'].iloc[i],transform=ax.transAxes,verticalalignment='top',
                family=plot_dict['legend_font'])

    ####



    ax.axhline(linewidth=1,color="black",zorder=0)
    ax.axvline(linewidth=1,color="black",zorder=0)


    ax.axhline(linewidth=.5,alpha=0.5,color="black",zorder=0,linestyle='--',y=plot_dict['cut_off'],dashes=(10,10))
    ax.axhline(linewidth=.5,alpha=0.5,color="black",zorder=0,linestyle='--',dashes=(10,10),y=-plot_dict['cut_off'])
    ax.axvline(linewidth=.5,alpha=0.5,color="black",zorder=0,linestyle='--',dashes=(10,10),x=plot_dict['cut_off'])
    ax.axvline(linewidth=.5,alpha=0.5,color="black",zorder=0,linestyle='--',dashes=(10,10),x=-plot_dict['cut_off'])

    f_max=plot_df['log2_for'].max()
    f_min=plot_df['log2_for'].min()
    r_max=plot_df['log2_rev'].max()
    r_min=plot_df['log2_rev'].min()

    if (plot_dict['x_lower_limit']!=None):
        
        ax.set_xlim(plot_dict['x_lower_limit'],plot_dict['x_upper_limit'])
        ax.set_ylim(plot_dict['y_lower_limit'],plot_dict['y_upper_limit'])
        
    else:
        abs_max=abs(max([f_max,f_min,r_max,r_min],key=abs))
        abs_max=math.ceil(abs_max)
        ax.set_xlim(-abs_max-0.5,abs_max+0.5)
        ax.set_ylim(-abs_max-0.5,abs_max+0.5)


    ax.tick_params(axis='both',labelsize=15)

    ax.set_xlabel(plot_dict['xlabel'],color='black',fontsize=plot_dict['xylabel_size'],family='arial')
    ax.set_ylabel(plot_dict['ylabel'],color='black',fontsize=plot_dict['xylabel_size'],family='arial')
    ax.set_title(plot_dict['title'],color='black',fontsize=plot_dict['title_size'],family='arial')

    no=1
    figure_name_1=new_directory+'/'+folder_name+'_silac_plot_'+str(no)+'_normalized_ratios.pdf'
    figure_name_2=new_directory+'/'+folder_name+'_silac_plot_'+str(no)+'_normalized_ratios.png'

    while os.path.isfile(figure_name_1) or os.path.isfile(figure_name_1) :
        no+=1
        figure_name_1=new_directory+'/'+folder_name+'_silac_plot_'+str(no)+'_non-normalized_ratios.pdf'
        figure_name_2=new_directory+'/'+folder_name+'_silac_plot_'+str(no)+'_non-normalized_ratios.png'


    plt.savefig(figure_name_1,transparent=True,bbox_inches="tight",dpi=300)
    plt.savefig(figure_name_2,transparent=True,bbox_inches="tight",dpi=300)



def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("protein_groups_file", help="Location of the protein group file",
                    type=str)
    parser.add_argument("plot_info", help="Location of the plot_info",
                    type=str)
    args = parser.parse_args() 

    directory=os.getcwd()+'/'
    print('Current directory is \n',directory)

    plot_dict=create_plot_dict(args.plot_info) 
    pgroups,folder_name,new_directory=read_protein_groups(args.protein_groups_file,directory)
    writer=pd.ExcelWriter(new_directory+'/'+folder_name+"_Processed_ProteinGroups.xlsx")
    pgroups.to_excel(writer,'Original')
    pgroups2=remove_unreliable_proteins(pgroups)
    pgroups2=pgroups2.sort_values(by='Ratio H/L normalized for',ascending=False,inplace=False)
    pgroups2.to_excel(writer,'Filtered')
    pgroups3=filter_for_rev(pgroups2) 
    
    cols=['Protein names', 'Gene names', 'Number of proteins', 'Peptides',  'Unique peptides', 'Peptides for', 
      'Peptides rev', 'Unique peptides for', 'Unique peptides rev', 'Sequence coverage [%]','Mol. weight [kDa]', 
      'Sequence length','Ratio H/L for', 'Ratio H/L normalized for','Ratio H/L rev', 'Ratio H/L normalized rev',
     'Intensity','Intensity L for', 'Intensity H for', 'Intensity L rev', 'Intensity H rev']
    
    pgroups4=select_cols(pgroups3,cols)
    pgroups4=pgroups4.sort_values(by=['Ratio H/L normalized for'],ascending=False)
    pgroups4.to_excel(writer,'Final')
    writer.save()
    plot_df=create_plot_df(pgroups4)
    legend_df=create_legend_df(plot_df,plot_dict)
    save_fig(plot_df,plot_dict,legend_df,new_directory,folder_name)

    plot_df_nNorm = create_plot_df_nNorm(pgroups4)
    legend_df_nNorm=create_legend_df(plot_df_nNorm,plot_dict)
    save_fig(plot_df_nNorm,plot_dict,legend_df_nNorm,new_directory,folder_name)

if __name__ == "__main__":
    main()