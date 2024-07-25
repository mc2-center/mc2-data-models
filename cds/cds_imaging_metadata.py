import yaml
import json
import os
import re
import requests
import pandas as pd
import numpy as np

import synapseclient
from google.cloud import bigquery
from openpyxl import load_workbook


def synapse_path(syn,synID,fileview): 
    # walk up parent chain to create Synapse_file_path
    synapse_path = fileview.loc[fileview['id'] == synID,'name'].iloc[0]
    parent_id = fileview.loc[fileview['id'] == synID,'parentId'].iloc[0]
    
    while parent_id not in list(htan_centers.values()):
        try:
            parent = syn.get(parent_id,downloadFile=False)
            synapse_path = parent.name+'/' + synapse_path
            parent_id = parent.parentId
        
        except TypeError as e:
            logging.error("entityId {} not found.".format(entity))
            return None
    
    return synapse_path


def get_channel_url(syn,package_id,cloud_path,syn_path):
    f = re.sub(r'[^a-zA-Z0-9\-\.\_\/]', '_', syn_path)
    bucket = cloud_path.split('/')[2:3][0]
    
    if 'proddata' in cloud_path:
        dest_folder = 'htan-dcc-srrs'
    else:
        dest_folder = bucket
    
    dest_url = 's3://cds-290-htan-x6gj-xr09/{0}/{1}/{2}'.format(
            package_id,dest_folder,f)
    
    return dest_url


def get_file_url(package_id,entityId,syn_path,fileview):
    f = re.sub(r'[^a-zA-Z0-9\-\.\_\/]', '_', syn_path)
    fv = fileview[fileview['id'] == entityId]
    bucket = fv['dataFileBucket'].values[0]
    
    if 'proddata' in fv['dataFileBucket'].values[0]:
        dest_folder = 'htan-dcc-srrs'
    else:
        dest_folder = bucket
    
    dest_url = 's3://cds-290-htan-x6gj-xr09/{0}/{1}/{2}'.format(
            package_id,dest_folder,f)
    
    return dest_url


def get_cds_value_set(cde_pvs, cds_attribute):
    '''
    Function to map HTAN attributes to similar CDS attributes containing commas
    '''

    pv_map = pd.DataFrame()
    pv_map['CDS'] = cde_pvs
    pv_map['HTAN'] = [x.replace(',','').replace(';','') for x in cde_pvs]
    
    return pv_map


def make_table(template_type, table, package_id, mappings, pv_df):

    num_files = len(table)
    cds = pd.DataFrame()

    for a in mappings[template_type]['attributes']:
        target_attribute = a['target_attribute']
        source_attribute = a['source_attribute']
        print(target_attribute)
        
        if 'fixed_value' in a:
            cds = cds.assign(**{
                target_attribute: [a['fixed_value']] * num_files
            })
        
        elif 'map' in a:

            r = eval(a['map'], {
                'source_col': table[
                    source_attribute.replace(' ','_')],
                }
            )

            if target_attribute in ['primary_diagnosis',
                'tissue_or_organ_of_origin','site_of_resection_or_biopsy']:
                i_s = int(pv_df.index[pv_df['Value Set Name'] == target_attribute].tolist()[0]) - 1 
                
                gen = (i for i, v in enumerate(pv_df['Term']) if v == None and i > i_s)
                i_e = next(gen)
                template_pv_list = pv_df['Term'][i_s:i_e]

                pv_map = get_cds_value_set(template_pv_list, target_attribute)
                r = [pv_map[pv_map['HTAN']==x]['CDS'].values[0] if x in list(
                    pv_map['HTAN']) else x for x in r]

            cds = cds.assign(**{target_attribute: r})
        
        else: 
            cds = cds.assign(**{target_attribute: None})

    cds = cds.loc[cds.astype(str).drop_duplicates().index]

    cds.to_csv(
        f'../tables/{package_id}/{template_type.replace(" ", "_")}_{package_id}.csv', 
        index=False
    )

def main():

    package_id = 'v24.3.1.img'

    print( ' ' )
    print( 'Preparing metadata tables for %s' % package_id )
    print( ' ' )

    # instantiate bigquery client
    client = bigquery.Client()

    syn=synapseclient.Synapse()
    syn.login()

    # read in mappings yaml file
    with open('../mapping/mappings.yaml', 'r') as f:
        mappings = yaml.safe_load(f)

    with open('./config.yaml', 'r') as file:
        config = yaml.safe_load(file)
   
    htan_centers = config['centers']

    # Open workbook
    wb = load_workbook(filename='../templates_and_dictionaries/CDS Imaging Templates-CDE-2023-08-31.xlsx')
    value_set = wb['Terms and Valuesets']

    pv_df = pd.DataFrame(value_set.values)
    pv_df = pv_df.rename(columns=pv_df.iloc[0]).loc[1:]

    url = 'https://raw.githubusercontent.com/ncihtan/htan-portal/'\
        'c6280bfc962b1c672c3b3244ffb6d03419b2f78e/data/human-organ-mappings.json' 
    resp = requests.get(url)
    primary_site_mappings = json.loads(resp.text)

    htan_center_mappings = json.load(
        open('../mapping/htan_center_mappings.json'))

    # create a directory for the current submission
    directory = os.getcwd().replace('scripts','tables')
    folder_path = os.path.join(directory, package_id)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    # ------------------------------------------------------------------------

    # pull in files and metadata from BQ
    img_metadata = client.query("""
        SELECT * EXCEPT(HTAN_Participant_ID)
        FROM `htan-dcc.combined_assays.ImagingLevel2`
        WHERE CDS_Release IS NULL
        AND Data_Release IS NOT NULL
    """).result().to_dataframe()

    srrs_img = client.query("""
        SELECT * EXCEPT(HTAN_Participant_ID)
        FROM `htan-dcc.combined_assays.SRRSImagingLevel2`
        WHERE CDS_Release IS NULL
        AND Data_Release IS NOT NULL
    """).result().to_dataframe()

    img_all = pd.concat([img_metadata, srrs_img])

    img_all = img_all.assign(
        **{'HTAN_Parent_Biospecimen_ID':img_all[
        'HTAN_Parent_Biospecimen_ID'].str.split(',')})
    img_all = img_all.explode('HTAN_Parent_Biospecimen_ID')

    fileview = syn.tableQuery(
        "SELECT id,name,parentId,dataFileConcreteType,dataFileBucket,dataFileKey \
        FROM syn20446927 WHERE type = 'file' \
        AND projectId NOT IN \
        ('syn21989705','syn20977135','syn20687304','syn32596076','syn52929270')"
    ).asDataFrame()

    released = client.query("""
        SELECT * FROM `htan-dcc.released.entities`
    """).result().to_dataframe()

    # get biospecimen and clinical metadata
    clinical = client.query("""
        SELECT DISTINCT * EXCEPT(HTAN_Center)
        FROM `htan-dcc.metadata.cds_clinical`
    """).result().to_dataframe()

    # merge assay and clinical data
    full_metadata = img_all.merge(clinical, 
        how = 'left', left_on ='HTAN_Parent_Biospecimen_ID', 
        right_on = 'HTAN_Biospecimen_ID'
    ).drop_duplicates('entityId', keep='last')

    for n in ['pi_first', 'pi_last', 'pi_email']:
        full_metadata[n] = [
        htan_center_mappings[y][n] for y in [
        x.split('_')[0] for x in list(full_metadata['HTAN_Participant_ID'])]
    ]

    # ------------------------------------------------------------------------

    full_metadata['transfer_id'] = package_id

    map_syn_channel = {}
    syn_channel_ids = full_metadata[full_metadata[
        'Channel_Metadata_Filename'].str.startswith('syn')]['Channel_Metadata_Filename']
    
    syn_channel_ids = list(set(syn_channel_ids))
    for c in syn_channel_ids:
        syn_path = synapse_path(syn, c, fileview)
        map_syn_channel[c] = syn_path
    
    full_metadata['Channel_Metadata_Filename'] = full_metadata[
        'Channel_Metadata_Filename'].replace(map_syn_channel)

    file_url_in_cds = []
    channel_metadata_url = []

    for i,r in full_metadata.iterrows():
        entityId = r['entityId']
        file_path = r['Filename']
        channel_path = r['Channel_Metadata_Filename']
        cloud_path = r['Cloud_Storage_Path']

        file_url = get_file_url(package_id, entityId, file_path, fileview)
        file_url_in_cds.append(file_url)

        channel_url  = get_channel_url(syn, package_id, cloud_path, channel_path)
        channel_metadata_url.append(channel_url)

    full_metadata['file_url_in_cds'] = file_url_in_cds
    full_metadata['channel_metadata_url'] = channel_metadata_url

    full_metadata['cancer_type'] = [(
        next((k for k,v in primary_site_mappings.items()
        if x in v['byTissueOrOrganOfOrigin']),
        'Not Reported')) for x in list(full_metadata[
        'Tissue_or_Organ_of_Origin'].apply(
        lambda x: str(x).replace(" NOS", ", NOS")))
    ]

    # convert age at diagnosis to years
    full_metadata['age_at_diagnosis_years'] = [int(float(x)/365)
        if (x not in [None, np.nan, 'unknown', 'Not Applicable', 'Not Reported'])
        else None for x in list(full_metadata['Age_at_Diagnosis'])]

    # mask ages <18 and >89
    full_metadata['age_at_diagnosis_years'] = [90 if y>89 else 18
        if y<18 else y for y in list(full_metadata['age_at_diagnosis_years'])]

    # For WorkingDistance, PhysicalSizeX, PhysicalSizeY, and PhysicalSizeZ
    # attributes, include only if <attribute>Unit is in microns
    for s in ['WorkingDistance','PhysicalSizeX','PhysicalSizeY','PhysicalSizeZ']:
        full_metadata[s] = [
        x if y in ['um','µm','<U+00B5>m'] else None for x,y in zip(
        list(full_metadata[s]),
        list(full_metadata[s+'Unit']))
    ]

    full_metadata.to_csv('../tables/%s/full_metadata.csv' % package_id , index=False)

    # subset by imaging type
    pathology = full_metadata[
        (full_metadata['Imaging_Assay_Type']=='H&E') |
        (full_metadata['Channel_Metadata_Filename']=='Not Applicable') |
        (full_metadata['HTAN_Center']=='HTAN SRRS')]

    multiplex = full_metadata[~full_metadata['entityId'].isin(pathology['entityId'].values)]

    print( ' ' )
    print( ' This package contains: ' )
    print( '    %s pathology files' % len(pathology) )
    print( '    %s multiplex microscopy files' % len(multiplex) )
    print( ' ' )

    templates = [
        'CDS Imaging Collection Information',
        'CDS Imaging Participant Information',
        'CDS Imaging File Specific General Information'
    ]
    for t in templates:
        print(f'Preparing {t}')
        make_table(t, full_metadata, package_id, mappings, pv_df)

    if len(pathology) > 0:
        t_path = 'CDS Imaging Non-DICOM Pathology'
        print(f'Preparing {t_path}')
        make_table(t_path, pathology, package_id, mappings, pv_df)

    if len(multiplex) > 0:
        t_mm = 'CDS Imaging Multiplex Microscopy'
        print(f'Preparing {t_mm}')
        make_table(t_mm, multiplex, package_id, mappings, pv_df)

    print( '' )
    print( 'Done!' )
    print( '' )


if __name__ == "__main__":

    main()