#!/usr/bin/env python
# coding: utf-8

# In[1]:


# 1 - LOAD THE PACKAGES

import pandas as pd
import os
import numpy as np
import sys


# In[2]:


# 2 - FILE LOCATION

files = [f for f in os.listdir('.') if os.path.isfile(f)]
files # Take the name of the file


# In[3]:


# 3 - LOADING ALL THE FILES INTO A DICTIONARY

dict_test = pd.read_excel('Technical Test - Data Wrangling.xlsx',sheet_name=None)


# In[4]:


# 4 - VISUALIZING THE SHEET NAMES

display(dict_test.keys())


# In[5]:


######################################## 5 DATA WRANGLING/SERUM PROTEIN DATA #########################################
# 5.1 - Melting Serum Protein data

SerumProteinDF = pd.melt(dict_test['Serum Protein data'], id_vars=['Patient','Sample'], value_vars=['Serum IL-6 (g/L)', 'Serum IL-6 Receptor (mg/L)'])


# In[6]:


# 5.2 - Curating the annotation of IL-6 and IL-6 Receptor

SerumProteinDF['variable'] = SerumProteinDF['variable'].str.replace(" Receptor", "R")
SerumProteinDF['variable'] = SerumProteinDF['variable'].str.replace("-", "")


# In[7]:


# 5.3 - Curating the annotation of measurement units

SerumAnno = SerumProteinDF['variable'].str.split(" ",n=2,expand=True)
SerumAnno.columns = ['Material','Gene_Symbol','Result_Units']

SerumAnno['Result_Units'] = SerumAnno['Result_Units'].str.replace("(", "")
SerumAnno['Result_Units'] = SerumAnno['Result_Units'].str.replace(")", "")


# In[8]:


# 5.4 - Updating Serum protein data with curated Gene_Symbol and Result_Units 

SerumProteinDF_2 = pd.concat([SerumProteinDF,SerumAnno],axis=1)
mapping = {'Serum':'SERUM'}
SerumProteinDF_2['Material'] = SerumProteinDF_2['Material'].map(mapping)
SerumProteinDF_2.head()


# In[9]:


# 5.5 - Subsetting the columns of interest

cols_interest = ['Sample','Gene_Symbol', 'value','Material','Result_Units']
SerumProteinDF_2[cols_interest].head()


# In[10]:


######################################## 6 DATA WRANGLING/ RNA-seq (RPKM) ############################################
# 6.1 - Exploring data
dict_test['RNA-seq (RPKM)'].head()


# In[11]:


# 6.1 - Transposing RNA-seq (RPKM) dataframe
RNAseq = dict_test['RNA-seq (RPKM)'].T
new_header = RNAseq.iloc[0]
RNAseq = RNAseq[1:]
RNAseq.columns = new_header
RNAseq.reset_index(inplace=True)
RNAseq.head()


# In[12]:


# 6.2 Melt RNAseq dataframe
RNAseq2 = pd.melt(RNAseq, id_vars='index', value_vars=['ICAM1', 'IL6','IL6R','VCAM1','SELE'])
RNAseq2.head()


# In[13]:


# 6.3 - Create Material and Result units columns
RNAseq2['Material'] = 'RNA'
RNAseq2['Result_Units'] = 'RPKM'
RNAseq2.head()


# In[14]:


# 6.4 - Renaming columns
RNAseq2 = RNAseq2.rename(columns={'index': 'Sample','GeneID':'Gene_Symbol'})
RNAseq2.head()


# In[15]:


######################################## 7 Merging RNAseq2 with SerumProteinDF_2 #####################################

RNASerum_df = pd.concat([RNAseq2,SerumProteinDF_2[cols_interest]])
RNASerum_df.head()


# In[16]:


############################################ 8 Create PatSample_ID #####################################################
# 8.1 - Renaming Tissue Sample Metadata columns
dict_test['Tissue Sample Metadata'] = dict_test['Tissue Sample Metadata'].rename(columns={'Patient  Number': 'Patient'})


# In[17]:


# 8.2 - Select Patient and Sample from Tissue Sample Metadata and SerumProteinDF_2 to create PatSample_id DF

pat_int = ['Patient','Sample']
PatSample_id = pd.concat([dict_test['Tissue Sample Metadata'][pat_int],SerumProteinDF_2[pat_int]])
PatSample_id


# In[18]:


############################################ 9 Create RNASerum_Patient DF ############################################
# 9.1 - Merge PatSample_id with RNASerum_df

RNASerum_Patient = pd.merge(PatSample_id,RNASerum_df,on='Sample')
RNASerum_Patient.head()


# In[19]:


############################################ 10 RNASerum_Patient_Study DF ##########################################################
# 10.1 - Merge Patient_clinical_data with RNASerum_Patient to create RNASerum_Patient_Study
RNASerum_Patient_Study = pd.merge(dict_test['Patient_clinical_data'],RNASerum_Patient,left_on='Patient  Number',right_on = 'Patient')
RNASerum_Patient_Study.head()


# In[20]:


# 10.2 - Create the Unique_Patient_ID
RNASerum_Patient_Study['Unique_Patient_ID'] = RNASerum_Patient_Study['Study_ID'] + '_' + RNASerum_Patient_Study['Patient  Number'].astype('str')
RNASerum_Patient_Study.head()


# In[21]:


####################### 11 - Creation of a dataframe in the specified format #########################################
# 11.1 - Selecting columns of interest in Tissue Sample Metadata

cols_interest_RNA = ['Patient','Sample','Sample type','Material']
dict_test['Tissue Sample Metadata'][cols_interest_RNA].head()


# In[22]:


# 11.2 - Merge RNASerum_Patient_Study with Tissue Sample Metadata to create the final_df

final_df = pd.merge(RNASerum_Patient_Study, dict_test['Tissue Sample Metadata'][cols_interest_RNA],  how='left', 
                    on=['Patient','Sample','Material'])
final_df.head()


# In[23]:


####################################### 12 - Customizing final_Df ####################################################
# 12.1 - Selecting the columns of interest
cols_interest = ['Study_ID','Patient','Unique_Patient_ID','Sex','Age','Sample','Sample type',
                 'Material','Gene_Symbol','value','Result_Units']

final_df_2 = final_df[cols_interest]
final_df_2.head()


# In[24]:


# 12.2 - Changing colnames
final_df_2.rename(columns={'Patient': 'Patient_ID','Sample':'Sample_ID','Sample type':'Sample_General_Pathology',
                          'Material':'Material_type','value':'Result'},inplace=True)
final_df_2.head()


# In[25]:


######################################### 13 - Modifying data types ##################################################
# 13.1 Transforming Result into a numeric variable
final_df_2['Result'] = pd.to_numeric(final_df_2['Result'], errors='coerce')


# In[26]:


# 13.2 Transforming Age into a numeric variable
final_df_2['Age'] = final_df_2['Age'].round().astype('int')


# In[27]:


######################################### 14 - Modifying feature entries #############################################
# 14.1 - Changing the sex entries
mapping = {'M':'MALE', 'F':'FEMALE'}
final_df_2['Sex'] = final_df_2['Sex'].map(mapping)


# In[28]:


# 14.2 - Changing Sample_General_Pathology entries
mapping2 = {'Normal':'NORMAL', 'Liver Tumor':'PRIMARY','Metastic Lung':'METASTATIC'}
final_df_2['Sample_General_Pathology'] = final_df_2['Sample_General_Pathology'].map(mapping2)
final_df_2['Sample_General_Pathology'].fillna('NA',inplace = True)


# In[29]:


# 14.3 - Curate Result_Units. Transform entries with mg/L as g/L

new_results = []
for res,units in zip(final_df_2['Result'],final_df_2['Result_Units']):
    if units == 'mg/L':
        new_results.append(res/1000)
    elif units == 'g/L':
        new_results.append(res)
    elif units == 'RPKM':
        new_results.append(res)

final_df_2['Result'] = new_results


# In[30]:


# 14.4 - Change the mg/L to g/L in Result_Units
mapping3 = {'mg/L':'g/L','RPKM':'RPKM','g/L':'g/L'}
final_df_2['Result_Units'] = final_df_2['Result_Units'].map(mapping3)


# In[31]:


# 14.5 - Creating the columns status

status = []
for result in final_df_2['Result']:
    if np.isnan(result):
        status.append('NOT DONE')
    else:
        status.append('NA')
final_df_2['Status'] = status


# In[32]:


# 14.6 - Capitalize Sample_ID
final_df_2['Sample_ID'] = final_df_2['Sample_ID'].str.upper()


# In[33]:


# 14.7 - Print the result
final_df_2


# In[34]:


################################################ 15 - FINAL REMARKS ##################################################
# 15.1 - Version of packages used
print('\n'.join(f'{m.__name__}=={m.__version__}' for m in globals().values() if getattr(m, '__version__', None)))


# In[35]:


# 15.2 - Python version
print(sys.version)


# In[36]:


################################################ 16 - OUTPUT TO CSV FILE ##############################################
final_df_2.to_csv('TechnicalTest_result.csv')


# In[ ]:




