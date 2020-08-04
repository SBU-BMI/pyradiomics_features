import json
import os
import numpy as np
import csv
import cv2
import sys
import shutil
import subprocess
import skimage.io
import time
from datetime import datetime
import logging
import concurrent.futures 
from multiprocessing import Process
from multiprocessing import cpu_count
import collections
from pymongo import MongoClient
import shapely.geometry
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
#import openslide 
from cal_pyradiomics import *
import glob
import pandas as pd



##################################################################
def read_poly_coor(poly_file, poly_field_name, label_field_name):
    poly_arr = [];
    label_arr = [];
    with open(poly_file,"rb") as csv_file:
        csv_reader = csv.reader(csv_file)
        line_count = 0
        for line_count, row in enumerate(csv_reader):
            if line_count == 0:
                # This is the title line
                polygonid = row.index(poly_field_name)
                labelid   = row.index(label_field_name)
                continue;

            # Read Polygons
            poly_str = row[polygonid];
            poly_str = poly_str[1:-1];  # This is to remove brackets at the two ends of the string
            poly_str = poly_str.split(':');

            # Map to long and subtract the offset
            poly_int = map(float, poly_str);
            poly_int = np.array(poly_int);
            poly_x = poly_int[range(0, len(poly_int), 2)];
            poly_y = poly_int[range(1, len(poly_int), 2)];
            #poly_x = poly_x.astype(np.int64);
            #poly_y = poly_y.astype(np.int64);

            # Combine x and y into array with format [[x1, y1], [x2, y2],...]
            poly = zip(poly_x, poly_y);
            poly = [list(t) for t in poly];
            if (len(poly_arr) == 0):
                poly_arr = [poly];
            else:
                poly_arr.append(poly);

            # Read labels
            lbl_str = int(row[labelid]);
            label_arr.append(lbl_str)

    return poly_arr, label_arr;
###############################################

    
####################################
def cal_patch(mask):
    # Nuclei ratio
    bin_mask = (mask > 0).astype(np.float);
    nuclei_area = np.sum(bin_mask);
    nuclei_ratio = nuclei_area / mask.size;

    # Nuclei perimeter
    polyidx_max = np.amax(mask);
    peri_arr = [];
    area_arr = []
    for poly_idx in range(polyidx_max):
        mask_nucleus = (mask == poly_idx+1).astype(np.uint8)*255;

        if (np.amax(mask_nucleus) == 0):
            # Meaning that this nucleus is not in this patch
            continue;

        # Get edge to compute perimeter
        edges = cv2.Canny(mask_nucleus, 100, 200);
        edges = (edges > 127).astype(np.uint16);
        peri_arr.append(np.sum(edges));

        # Compute area of this nucleus
        mask_nucleus_bin = (mask_nucleus > 127).astype(np.uint16);
        area_arr.append(np.sum(mask_nucleus_bin));

    if (len(area_arr) > 0):
        area_arr = np.array(area_arr);
        peri_arr = np.array(peri_arr);
        return nuclei_ratio, np.mean(area_arr), np.mean(peri_arr);
    else:
        return nuclei_ratio, 0, 0;
##################################################################


############################################################
def cal_patch_tissue(mask, tissue):
    # Mask the nuclei mask with tissue mask
    mask = np.multiply(mask, tissue);

    # Nuclei ratio
    bin_mask = (mask > 0).astype(np.float);
    nuclei_area = np.sum(bin_mask);
    if np.sum(tissue) > 0:
        nuclei_ratio = nuclei_area / np.sum(tissue);
    else:
        nuclei_ratio = 0.0;

    # Nuclei perimeter
    polyidx_max = np.amax(mask);
    peri_arr = [];
    area_arr = []
    for poly_idx in range(polyidx_max):
        mask_nucleus = (mask == poly_idx+1).astype(np.uint8)*255;

        if (np.amax(mask_nucleus) == 0):
            # Meaning that this nucleus is not in this patch
            continue;

        # Get edge to compute perimeter
        edges = cv2.Canny(mask_nucleus, 100, 200);
        edges = (edges > 127).astype(np.uint16);
        peri_arr.append(np.sum(edges));

        # Compute area of this nucleus
        mask_nucleus_bin = (mask_nucleus > 127).astype(np.uint16);
        area_arr.append(np.sum(mask_nucleus_bin));

    if (len(area_arr) > 0):
        area_arr = np.array(area_arr);
        peri_arr = np.array(peri_arr);
        return nuclei_ratio, np.mean(area_arr), np.mean(peri_arr);
    else:
        return nuclei_ratio, 0, 0;
########################################################################


##########################################################################################
def process_tile(csvfile, jsonfile, svsfile, outfiledir, outtiledir, outtissuedir):    
    # For pyradiomics
    feature_setting_file = 'pyrad_features.txt'
    feature_setting_file_path = os.path.join(code_base, feature_setting_file);    
    
    pyradiomics_ins = PyRadiomics_Features(feature_setting_file_path);                 
    
    # Param
    patch_width_org = int(PATCH_SIZE);
    patch_height_org = int(PATCH_SIZE);    
    tile_name = os.path.basename(csvfile).split('-features')[0];
    outfilepath = os.path.join(outfiledir, tile_name+'-stat.csv');

    # Read json file   
    with open(jsonfile) as f:
        datajson = json.load(f);        
        img_width = int(datajson['image_width']);
        img_height = int(datajson['image_height']);
        tile_x = int(datajson['tile_minx']);
        tile_y = int(datajson['tile_miny']);
        tile_width = int(datajson['tile_width']);
        tile_height = int(datajson['tile_height']);
        mpp = datajson['mpp'];
        caseid = datajson['subject_id'];      
      
    # Read csv file    
    poly_arr = [];    
    with open(csvfile, newline='') as csv_file:
        csv_reader = csv.reader(csv_file)
        line_count = 0
        for line_count, row in enumerate(csv_reader):
            if line_count == 0:
                # This is the title line
                polygonid =  row.index('Polygon')
                continue;
            poly_str = row[polygonid];
            poly_str = poly_str[1:-1];  # This is to remove brackets at the two ends of the string
            poly_str = poly_str.split(':');
            # Map to long and subtract the offset            
            poly_int = map(float, poly_str);            
            poly_int = list(poly_int)           
            poly_int = np.array(poly_int);           
            poly_x = poly_int[range(0, len(poly_int), 2)];
            poly_y = poly_int[range(1, len(poly_int), 2)];
            poly_x = (poly_x - tile_x).astype(np.int64);
            poly_y = (poly_y - tile_y).astype(np.int64);

            # Combine x and y into array with format [[x1, y1], [x2, y2],...]
            poly = zip(poly_x, poly_y);
            poly = [list(t) for t in poly];
            if len(poly_arr) == 0:
                poly_arr = [poly];
            else:
                poly_arr.append(poly);      
      
    # Convert polygons to mask    
    mask = np.zeros((tile_height, tile_width), dtype=np.uint16)
    for poly_idx, single_poly in enumerate(poly_arr):
      cv2.fillConvexPoly(mask, np.array(single_poly), (poly_idx+1));        
    
    image,tissue_mask = get_tissue_mask(jsonfile, svsfile, outtiledir,tile_x,tile_y, tile_width, tile_height, outtissuedir);    
         
    # Break the tile into patches
    with open(outfilepath, 'w') as outfile:
        writer = csv.writer(outfile);
        wrtline = ['case_id', 'image_width', 'image_height', 'mpp_x', 'mpp_y', 'patch_x', 'patch_y', 'patch_width', 'patch_height',
                   'patch_area_micro', 'nuclei_area_micro', 'nuclei_ratio', 'nuclei_average_area', 'nuclei_average_perimeter', 'pseudo-feat'];
        
        # For pyradiomics
        wrtline = wrtline + ['fg_' + s for s in pyradiomics_ins.get_feature_name_list()] + ['bg_' + s for s in pyradiomics_ins.get_feature_name_list()]       
          
        writer.writerow(wrtline);
        
        feature_name_array=pyradiomics_ins.get_feature_name_list();
        
        for patch_x in range(0,tile_width,patch_width_org):
            for patch_y in range(0, tile_height, patch_height_org):
                patch_width  = min(tile_width - patch_x, patch_width_org);
                patch_height = min(tile_height - patch_y, patch_height_org);
                if ((patch_width == 0) or (patch_height == 0)):
                    continue;
                patch = mask[patch_y : patch_y + patch_height, patch_x : patch_x + patch_width];
                tissue_patch = tissue_mask[patch_y: patch_y + patch_height, patch_x: patch_x + patch_width];
                image_patch = image[patch_y: patch_y + patch_height, patch_x: patch_x + patch_width];
                nu_ratio, nu_area, nu_peri = cal_patch_tissue(patch, tissue_patch);
                patch_area_micro = patch_width * patch_height * mpp * mpp;
                nu_area_micro = patch_area_micro * nu_ratio;
                nucleus_material_percentage= nu_ratio*100.0;         
                              
                wrtline = [caseid, img_width, img_height, mpp, mpp, patch_x+tile_x, patch_y+tile_y, patch_width, patch_height,
                           patch_area_micro, nu_area_micro, nu_ratio, nu_area, nu_peri];
                
                # Pseudo Feat
                pseu_feat = (patch_y + tile_y) * img_width + (patch_x + tile_x);
                wrtline = wrtline + [str(pseu_feat)];
                pseud_feat = str(pseu_feat); 
                
                # For pyradiomics 
                fg_patch = (patch > 0).astype(np.uint8)
                fg_patch = np.multiply(fg_patch, tissue_patch)
                pyrad_feat_fg = pyradiomics_ins.cal_pyradiomics(image_patch, fg_patch);
                wrtline = wrtline + pyrad_feat_fg;                
                bg_patch =  1-(patch > 0).astype(np.uint8)
                bg_patch = np.multiply(bg_patch, tissue_patch)
                pyrad_feat_bg = pyradiomics_ins.cal_pyradiomics(image_patch, bg_patch);
                wrtline = wrtline + pyrad_feat_bg;                 
                writer.writerow(wrtline); 
      
                radiomics_json={};
                for index,feature_name in enumerate(feature_name_array):
                  tmp_name_fg='fg_' + str(feature_name); 
                  tmp_name_bg='bg_' + str(feature_name);
                  tmp_fg_value=pyrad_feat_fg[index];
                  tmp_bg_value=pyrad_feat_bg[index];
                  radiomics_json[tmp_name_fg]=tmp_fg_value;
                  radiomics_json[tmp_name_bg]=tmp_bg_value;                         
##############################################################################################


#######################################################################################################
def get_tissue_mask(jsonfile, svsfile, outtiledir,tile_x,tile_y, tile_width, tile_height, outtissuedir):
    # Param
    verification_patch_size = 128;

    # Check and extract tile
    tile_name = os.path.basename(jsonfile).split('-algmeta.json')[0];
    rgb_tile_path = os.path.join(outtiledir, tile_name + '-tile.png');    
    #xcoor = int(tile_name.split('_')[-2][1:]);
    #ycoor = int(tile_name.split('_')[-1].split('-')[0][1:]);
    xcoor=tile_x;
    ycoor=tile_y;
    #print "xcoor,ycoor";
    #print xcoor,ycoor;
    if not os.path.isfile(rgb_tile_path):
        cmd = 'openslide-write-png {} {} {} 0 {} {} {}'.format(svsfile, xcoor, ycoor, tile_width, tile_height, rgb_tile_path);
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE);
        process.wait();
    
    #print rgb_tile_path; 
    # Read the rgb file
    org_img = skimage.io.imread(rgb_tile_path).astype(np.float);
    img = org_img / 255;  # Divide by 255 to convert to [0,1]
    tissue_mask = np.zeros((img.shape[0], img.shape[1]), dtype=np.uint8);

    # Identify glass for each patch
    for patch_x in range(0, tile_width, verification_patch_size):
        for patch_y in range(0, tile_height, verification_patch_size):
            patch_width = min(tile_width - patch_x, verification_patch_size);
            patch_height = min(tile_height - patch_y, verification_patch_size);
            if (patch_width == 0) or (patch_height == 0):
                continue;
            patch = img[patch_y: patch_y + patch_height, patch_x: patch_x + patch_width, :];
            wh = patch[..., 0].std() + patch[..., 1].std() + patch[..., 2].std();
            if wh >= 0.20:
                tissue_mask[patch_y: patch_y + patch_height, patch_x: patch_x + patch_width] = 1;

    # Save tissue mask
    if outtissuedir is not None:
        # Save this tissue mask
        tissuefile = os.path.join(outtissuedir, tile_name + '-tissue.png');
        skimage.io.imsave(tissuefile, tissue_mask*255);

    return org_img,tissue_mask;
##########################################################################################


#######################################
def handMultiPolygon(polygon_list):   
  tmp_polygon_list=[];  
  for poly in polygon_list:    
    if poly.geom_type  == 'Polygon':
      tmp_polygon_list.append(poly);
    elif poly.geom_type  == 'MultiPolygon':      
      new_shape = shapely.ops.cascaded_union(poly);      
      if new_shape.geom_type == 'MultiPolygon':   
        #print  "-- find MultiPolygon --";    
        for polygon in new_shape:       
          ext_polygon_points =list(zip(*polygon.exterior.coords.xy)); 
          newPoly = Polygon(ext_polygon_points);
          tmp_polygon_list.append(newPoly);          
      elif new_shape.geom_type == 'Polygon':
        tmp_polygon_list.append(new_shape);         
  return tmp_polygon_list;  
######################################

################################################
def getCsvRecordCount(csv_file_path):
  line_counts=0;
  #print "csv_file_path"+ str(csv_file_path);
  with open(csv_file_path,newline='') as csv_file:
    csv_reader = csv.reader(csv_file) ;    
    try:
       for line_count, row in enumerate(csv_reader):
         line_counts=line_count+1; 
         #print  line_count;   
    except Exception as e: 
      print(e);  
      print ("error occurs while dealing this csv file.");
      print (csv_file_path);     
  #print  line_counts 
  return  line_counts;     
#################################################

##################################################
def isTileRelated2TumorRegion(json_file_path,humanMarkupList_tumor):
  with open(json_file_path) as f:
    datajson = json.load(f);
    img_width = datajson['image_width'];
    img_height = datajson['image_height'];
    tile_x = datajson['tile_minx'];
    tile_y = datajson['tile_miny'];
    tile_width = datajson['tile_width'];
    tile_height = datajson['tile_height'];
    
    x1=float(tile_x)/float(image_width);
    y1=float(tile_y)/float(image_height); 
    x2=float(tile_x+tile_width)/float(image_width);
    y2=float(tile_y+tile_height)/float(image_height);
    if x1>1.0:
      x1=1.0;
    if x1<0.0:
      x1=0.0; 
    if x2>1.0:
      x2=1.0;            
    if x2<0.0:
      x2=0.0;            
    if y1>1.0:
      y1=1.0;
    if y1<0.0:
      y1=0.0; 
    if y2>1.0:
      y2=1.0;
    if y2<0.0:
      y2=0.0;  
    tile_polygon_0=[[x1,y1],[x2,y1],[x2,y2],[x1,y2],[x1,y1]];  
    tmp_poly=[tuple(i) for i in tile_polygon_0];
    tmp_polygon = Polygon(tmp_poly);
    tile_polygon = tmp_polygon.buffer(0);
    
    tileHumanMarkupRelation_tumor="disjoin";   
    
    for humanMarkup in humanMarkupList_tumor:                         
      if (tile_polygon.within(humanMarkup)):              
        tileHumanMarkupRelation_tumor="within";        
        break;
      elif (tile_polygon.intersects(humanMarkup)):                
        tileHumanMarkupRelation_tumor="intersect";          
        break;
      else:               
        tileHumanMarkupRelation_tumor="disjoin";            
                      
    #only calculate features within/intersect tumor region           
    if(tileHumanMarkupRelation_tumor=="disjoin"):                     
      return False;            
    else: 
      return True;      
#################################################

################################################
def getImageMetaData(segmentResultPath):
  json_filename_list = [f for f in os.listdir(segmentResultPath) if f.endswith('.json')] ; 
  for json_filename in json_filename_list[0:1]:                              
      json_file_path=os.path.join(segmentResultPath, json_filename); 
      # Read json file
      with open(json_file_path) as f:
        datajson = json.load(f);
        img_width = datajson['image_width'];
        img_height = datajson['image_height'];
        
        return img_width, img_height;
###############################################


################################################
def findTumor_NonTumorRegions(case_id):
  tumor_region_json_path="/tumor_region/";            
  humanMarkupList_tumor=[];
  tmp_tumor_markup_list0=[];  
  		
  json_filename_list = [f for f in os.listdir(tumor_region_json_path) if f.endswith('.json')] ; 
  if len(json_filename_list) > 0:    	
    print ("tumor region json files are available.");  
    for json_file in json_filename_list:
      json_file_detail = os.path.join(tumor_region_json_path, json_file);       
      with open(json_file_detail) as f:
        file_data = json.load(f); 
        tmp_case_id=file_data['provenance']['image']['case_id'];
        if tmp_case_id == case_id:
          humarkup_polygon_tmp=file_data['geometry']['coordinates'][0];
          tmp_polygon=[tuple(i) for i in humarkup_polygon_tmp];
          tmp_polygon=Polygon(tmp_polygon); 
          original_area=tmp_polygon.area;         
          tmp_polygon2=tmp_polygon.buffer(0); 
          area=tmp_polygon2.area;          
          factor=float(area)/float(original_area);
          if factor<0.5 or factor>2.0:      
            tmp_polygon2=tmp_polygon.convex_hull;                                                     
          tmp_tumor_markup_list0.append(tmp_polygon2);    
        else:
          continue;
                    
  #handle MultiPolygon
  tmp_tumor_markup_list=handMultiPolygon(tmp_tumor_markup_list0);   
   
  #merge polygons if applicale            
  index_intersected=[];                                
  for index1 in range(0, len(tmp_tumor_markup_list)):  
    if index1 in index_intersected :#skip polygon,which is been merged to another one
      continue;
    humarkup_polygon1=tmp_tumor_markup_list[index1];         
    is_within=False;
    is_intersect=False;
    for index2 in range(0, len(tmp_tumor_markup_list)):  
      humarkup_polygon2=tmp_tumor_markup_list[index2];
      if (index1 != index2):
        if (humarkup_polygon1.within(humarkup_polygon2)):    
          is_within=True;            
          break;              
        if (humarkup_polygon1.intersects(humarkup_polygon2)):
          humarkup_polygon1=humarkup_polygon1.union(humarkup_polygon2);           
          is_intersect=True;
          index_intersected.append(index2);                
    if(not is_within and not is_intersect):
      humanMarkupList_tumor.append(humarkup_polygon1);          
    if(is_within):
      continue;         
    if(is_intersect):          
      humanMarkupList_tumor.append(humarkup_polygon1); 
      
  return  humanMarkupList_tumor;        
###############################################


##################################################################### 
if __name__ == "__main__":    
  if len(sys.argv)<1:
    print ("usage:python run_radiomics_all_patches.py case_id");
    exit();    
  
  case_id= sys.argv[1];
  print ("case_id");
  print (case_id);  
  
  current_time=datetime.now();
  LOG_FILENAME = 'error_jonas.log'
  logging.basicConfig(filename=LOG_FILENAME,level=logging.DEBUG) 
  
  code_base=os.getcwd();   
  
  csv.field_size_limit(sys.maxsize);
  max_workers = cpu_count() - 1  # Number of processors to use, keep 1 processores free for other work
  if max_workers < 1:  # in case only one processor is available, ensure that it is used
    max_workers = 1    
    
  #get env variables
  PATCH_SIZE = os.environ['PATCH_SIZE']; 
  image_file_format=os.environ['IMAGE_FORMAT'];  
  HAS_TUMOR_REGION= os.environ['HAS_TUMOR_REGION'];     

  #mapping the host input/output folders to local variables
  imagefiles_path="/image_files";
  segment_results_path="/segment_results";  
  output_path="/output/";  
  
  image_file_format=image_file_format.replace("'", "");
  file_end='.'+str(image_file_format);  
  imageFileName=imagefiles_path + '/'+str(case_id)+file_end       
  #segmentResultPath=segment_results_path + '/'+str(case_id)+'/';
  segmentResultPath=segment_results_path + '/'+str(case_id)+ file_end +'/';
  #print(image_file_format,file_end,case_id,imageFileName);
  
  csv_folder = os.path.join(output_path, 'patch_level_csv'); 
  if not os.path.exists(csv_folder):
    print ('%s folder do not exist, then create it.' % csv_folder);
    os.makedirs(csv_folder); 

  output_folder = os.path.join(csv_folder, case_id); 
  if not os.path.exists(output_folder):
    print ('%s folder do not exist, then create it.' % output_folder);
    os.makedirs(output_folder);
        
  output_tile_folder = os.path.join(output_folder, 'tiles');
  if not os.path.exists(output_tile_folder):
    print ('%s folder do not exist, then create it.' % output_tile_folder);
    os.makedirs(output_tile_folder);
    
  output_tissue_folder=output_tile_folder;    

  if os.path.isdir(segmentResultPath) and len(os.listdir(segmentResultPath)) > 0: 
    print ( " all csv and json files of this image are here.");
  else:      
    print ("No segment results is found."); 
    exit();  
         
  if not os.path.isfile(imageFileName):
    print ("image svs file is not available.");            
    exit()   
  
  image_width, image_height = getImageMetaData(segmentResultPath);  
  
  if HAS_TUMOR_REGION=="yes":
    #get human markup data    
    humanMarkupList_tumor=findTumor_NonTumorRegions(case_id);   
    if(len(humanMarkupList_tumor) ==0):
      print ("No tumor regions has been marked in this image.");
      exit();      
    
  with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:    
    json_filename_list = [f for f in os.listdir(segmentResultPath) if f.endswith('.json')] ; 
    for json_filename in json_filename_list:                              
      json_file_path=os.path.join(segmentResultPath, json_filename);        
      csv_file=json_filename.replace('-algmeta.json','-features.csv');
      csv_file_path=os.path.join(segmentResultPath, csv_file);
     
      if not os.path.isfile(csv_file_path):
        print ("can NOT find matching csv file from json file", case_id,json_filename);
        continue;
        
      csv_record_count=getCsvRecordCount(csv_file_path);           
      if csv_record_count==1: #no nucleus material in this tile, only header available, so skip it
        print ("no nucleus material in this tile, only header available, so skip it");
        continue;   
                 
      if csv_record_count==0: #error occurs while dealing with csv data file, so skip it
        print ( "error occurs while dealing with csv data file, so skip it");
        continue;    
         
      #only calculate tile related to tumor or non tumor region      
      if HAS_TUMOR_REGION=="yes": 
        if not isTileRelated2TumorRegion(json_file_path,humanMarkupList_tumor): 
          print ("This tile is not within/intersect the tumor region.So skip it.");           
          continue;  
        
      print (" --- is processing csv file " + str(csv_file));
      try:    
        executor.submit(process_tile,csv_file_path, json_file_path, imageFileName, output_folder, output_tile_folder, output_tissue_folder);         
      except Exception as e: 
        print(e);          
        continue;     
  
  #merge cvs files from each tile to a single csv file
  os.chdir(output_folder);  
  extension = 'csv'
  all_filenames = [i for i in glob.glob('*.{}'.format(extension))]
  #combine all files in the list
  combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames ])
  #export to csv
  combined_csv.to_csv( "patch_level_radiomics_features.csv", index=False, encoding='utf-8-sig')
                                                 
  exit(); 
