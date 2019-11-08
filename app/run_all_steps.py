from pymongo import MongoClient
import subprocess
import os
import sys
import csv
import json
import collections
import datetime
import numpy as np
import matplotlib.pyplot as plt; plt.rcdefaults()
from textwrap import wrap
import fnmatch
import operator


########################################
def generate_Patch_Level_Radiomics_Feature_Plots(feature_array):      
  for case_id in os.listdir(patch_level_csv_folder):
    case_id_folder=os.path.join(patch_level_csv_folder, case_id);
    csv_file=os.path.join(case_id_folder, 'patch_level_radiomics_features.csv');  
    if not os.path.isfile(csv_file):
      continue;
    picture_folder = os.path.join(patch_level_plot_folder, case_id); 
    if not os.path.exists(picture_folder):
      print ('%s folder do not exist, then create it.' % picture_folder);
      os.makedirs(picture_folder);
    elif len(os.listdir(picture_folder)) ==0:
      print ("No plot is available, so continue");  
    else:
      continue;#done it before,skip this image
      
    for feature in feature_array: 
      with open(csv_file,newline='') as csvFile:
        csv_reader = csv.reader(csvFile) ;  
        feature_value_array=[]; 
        feature_id=0;
        for line_count, row in enumerate(csv_reader):
          if line_count == 0:# This is the title row
            feature_id=row.index(feature);                        
          else:# rest of data rows
            value= row[feature_id];
            if value!='None' and value!='0.0':              
              try:
                float_value=float(value);
                feature_value_array.append(float_value); 
              except Exception as e: 
                print (e);
                print (value);
                continue;         
        total_patch_count=len(feature_value_array); 
        print (case_id,feature,total_patch_count);
        if len(feature_value_array) >0: 
          fig, ax = plt.subplots();    
          n, bins, patches = plt.hist(feature_value_array, bins='auto',facecolor='blue');      
          plt.xlabel(feature)
          plt.ylabel('Patch Count')
          plt.title("\n".join(wrap("patch level "+ feature+ ' Histogram of image '+ str(case_id))))
          #Tweak spacing to prevent clipping of ylabel
          plt.subplots_adjust(left=0.15)
          plt.grid(True);           
          # place a text box in upper left in axes coords
          props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
          total_patch_count="{:,}".format(total_patch_count)
          textstr="Total patch count: " + str(total_patch_count);
          ax.text(0.6, 0.95, textstr, transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props);      
          #plt.show();
          file_name="patch_level_histogram_"+case_id+"_"+feature+".png";  
          graphic_file_path = os.path.join(picture_folder, file_name);
          plt.savefig(graphic_file_path); 
          plt.gcf().clear(); 
          plt.close('all');          
#####################################



#####################################
def generate_image_level_radiomics_feature_plots(feature_array):
  for case_id in os.listdir(patch_level_csv_folder):
    case_id_folder=os.path.join(patch_level_csv_folder, case_id);
    csv_file=os.path.join(case_id_folder, 'patch_level_radiomics_features.csv');
    if not os.path.isfile(csv_file):
      continue;
    picture_folder = os.path.join(image_level_plot_folder, case_id); 
    if not os.path.exists(picture_folder):
      print ('%s folder do not exist, then create it.' % picture_folder);
      os.makedirs(picture_folder);
    elif len(os.listdir(picture_folder)) ==0:
      print ("No plot is available, so continue");  
    else:
      continue;#done it before,skip this image
      
    for feature in feature_array: 
      with open(csv_file,newline='') as csvFile:
        csv_reader = csv.reader(csvFile) ;  
        feature_value_array=[]; 
        feature_id=0;
        for line_count, row in enumerate(csv_reader):
          if line_count == 0:# This is the title row
            feature_id=row.index(feature);                        
          else:# rest of data rows
            value= row[feature_id];
            if value!='None' and value!='0.0':              
              try:
                float_value=float(value);
                feature_value_array.append(float_value); 
              except Exception as e: 
                print (e);
                print (value);
                continue;
                            
      percentile_10th = np.percentile(feature_value_array,10);
      percentile_25th = np.percentile(feature_value_array,25);
      percentile_50th = np.percentile(feature_value_array,50); 
      percentile_75th = np.percentile(feature_value_array,75);
      percentile_90th = np.percentile(feature_value_array,90);   
      objects = ('10th', '25th', '50th', '75th', '90th')
      y_pos = np.arange(len(objects))
      percentile = [percentile_10th,percentile_25th,percentile_50th,percentile_75th,percentile_90th]; 
      plt.bar(y_pos, percentile, align='center', alpha=0.5)
      plt.xticks(y_pos, objects)
      plt.ylabel(feature)       
      plt.title("\n".join(wrap("patient level "+ feature+ ' percentile of image '+ str(case_id))))
      plt.subplots_adjust(left=0.15)
      plt.grid(True);
      #plt.show()
      file_name="image_level_percentile_"+case_id+"_"+feature+".png";  
      graphic_file_path = os.path.join(picture_folder, file_name);
      plt.savefig(graphic_file_path); 
      plt.gcf().clear();  
      plt.close('all');
#####################################


#######################################
def generate_image_level_radiomics_feature_csv_file(feature_array): 
  if len(os.listdir(image_level_csv_folder)) ==0:
    print ("No image level csv file is available, so continue");  
  else:
    return;#done it before,skip this image      
    
  feature_array2=[];     
  feature_array2.append("case_id") ;
          
  for feature in  feature_array: 
    tmp_feature1=str(feature)+ "(10th)" 
    feature_array2.append(tmp_feature1) ;
    tmp_feature2=str(feature)+ "(25th)" 
    feature_array2.append(tmp_feature2) ;
    tmp_feature3=str(feature)+ "(50th)" 
    feature_array2.append(tmp_feature3) ;
    tmp_feature4=str(feature)+ "(75th)" 
    feature_array2.append(tmp_feature4) ;
    tmp_feature5=str(feature)+ "(90th)" 
    feature_array2.append(tmp_feature5) ;      
    
  for case_id in os.listdir(patch_level_csv_folder):
    case_id_folder=os.path.join(patch_level_csv_folder, case_id);
    csv_file=os.path.join(case_id_folder, 'patch_level_radiomics_features.csv');
    if not os.path.isfile(csv_file):
      continue;
    file_name="image_level_radiomics_feature_"+case_id+".csv";
    print ("create file "+str(file_name));
    is_new_file=False;    
    if os.path.isdir(image_level_csv_folder) and len(os.listdir(image_level_csv_folder)) > 0:                            
      csv_filename_list = [f for f in os.listdir(image_level_csv_folder) if f.endswith('.csv')] ;  
      if file_name not in csv_filename_list: 
        is_new_file=True;
    elif os.path.isdir(image_level_csv_folder) and len(os.listdir(image_level_csv_folder)) == 0:    
      is_new_file=True;
    
    if is_new_file:          
      csv_dest_file = os.path.join(image_level_csv_folder, file_name);
      with open(csv_dest_file, 'w') as csv_write: 
        csv_writer = csv.writer(csv_write);
        csv_writer.writerow(feature_array2); 
        content_row=[];
        content_row.append(case_id);
        for feature in feature_array: 
          with open(csv_file,newline='') as csvFile:
            csv_reader = csv.reader(csvFile) ;  
            feature_value_array=[]; 
            feature_id=0;
            for line_count, row in enumerate(csv_reader):
              if line_count == 0:# This is the title row
                feature_id=row.index(feature);                        
              else:# rest of data rows
                value= row[feature_id];
                if value!='None' and value!='0.0':              
                  try:
                    float_value=float(value);
                    feature_value_array.append(float_value); 
                  except Exception as e: 
                    print (e);
                    print (value);
                    continue;
                            
          percentile_10th = np.percentile(feature_value_array,10);
          percentile_25th = np.percentile(feature_value_array,25);
          percentile_50th = np.percentile(feature_value_array,50); 
          percentile_75th = np.percentile(feature_value_array,75);
          percentile_90th = np.percentile(feature_value_array,90); 
          #print (feature,percentile_10th,percentile_25th,percentile_50th,percentile_75th,percentile_90th);
          content_row.append(percentile_10th);  
          content_row.append(percentile_25th); 
          content_row.append(percentile_50th); 
          content_row.append(percentile_75th); 
          content_row.append(percentile_90th); 
          #print (content_row)
        csv_writer.writerow(content_row) 
#######################################



##########################################
def generate_image_level_nuclei_ratio_csv_file(): 
  file_name="image_level_nuclei_ratio.csv";     
  feature_array0=["case_id","total_nuclei_area_micro","total_tumor_area_micro","nuclei_tumor_ratio"]; 
  csv_dest_file = os.path.join(image_level_csv_folder, file_name);
  if os.path.exists(csv_dest_file):
    print ("image_level_nuclei_ratio.csv existes, exit.");
    return;
    
  with open(csv_dest_file, 'w') as csv_write: 
    csv_writer = csv.writer(csv_write);
    csv_writer.writerow(feature_array0);             
    
    for case_id in os.listdir(patch_level_csv_folder):
      case_id_folder=os.path.join(patch_level_csv_folder, case_id);
      csv_file=os.path.join(case_id_folder, 'patch_level_radiomics_features.csv');  
      if not os.path.isfile(csv_file):
        continue;
      total_nuclei_area_micro=0.0;
      total_tumor_area_micro=0.0;
      total_patch_area_micro=0.0;  
      
      with open(csv_file,newline='') as csvFile:
        csv_reader = csv.reader(csvFile) ;  
        feature_value_array=[]; 
        nuclei_ratio_id=0;
        nuclei_area_micro_id=0;        
        for line_count, row in enumerate(csv_reader):
          if line_count == 0:# This is the title row
            nuclei_ratio_id=row.index('nuclei_ratio'); 
            nuclei_area_micro_id=row.index('nuclei_area_micro'); 
            patch_area_micro_id=row.index('patch_area_micro');                                              
          else:# rest of data rows
            nuclei_ratio_value= row[nuclei_ratio_id];
            nuclei_area_micro_value= row[nuclei_area_micro_id]; 
            patch_area_micro_value= row[patch_area_micro_id];           
            if nuclei_ratio_value!='None' and nuclei_ratio_value!='0.0':              
              try:
                nuclei_ratio=float(nuclei_ratio_value);
                nuclei_area_micro=float(nuclei_area_micro_value);
                patch_area_micro=float(patch_area_micro_value);                
                tumor_area_micro=nuclei_area_micro/nuclei_ratio;                 
                total_nuclei_area_micro=total_nuclei_area_micro+nuclei_area_micro;  
                total_tumor_area_micro=total_tumor_area_micro+ tumor_area_micro;
                total_patch_area_micro=total_patch_area_micro+patch_area_micro; 
              except Exception as e: 
                print (e);
                print (value);
                continue;  
      
        nuclei_tumor_ratio = total_nuclei_area_micro/ total_tumor_area_micro;
        nuclei_patch_ratio = total_nuclei_area_micro / total_patch_area_micro;
        row=[];
        row.append(case_id); 
        row.append(total_nuclei_area_micro);
        row.append(total_tumor_area_micro);        
        row.append(nuclei_tumor_ratio);        
        csv_writer.writerow(row);
########################################


if __name__ == '__main__':
  if len(sys.argv)<0:
    print ("usage:python run_all_steps.py");
    exit(); 
  
  #get env variables
  PATCH_SIZE = os.environ['PATCH_SIZE'];
  IMAGE_FORMAT= os.environ['IMAGE_FORMAT']; 
  HAS_TUMOR_REGION= os.environ['HAS_TUMOR_REGION'];     
   
  #get code base  
  code_base = os.getcwd();   
    
  #mapping the host input/output folders to local variables    
  imagefiles_path="/image_files/";
  segment_results_path="/segment_results/";  
  output_path="/output/"; 
  
  #create subfolder for output if not existing
  patch_level_csv_folder = os.path.join(output_path,'patch_level_csv');    
  patch_level_plot_folder= os.path.join(output_path,'patch_level_plot');
  image_level_csv_folder= os.path.join(output_path,'image_level_csv');
  image_level_plot_folder= os.path.join(output_path,'image_level_plot');
    
  if not os.path.exists(patch_level_csv_folder):
    print ('%s folder do not exist, then create it.' % patch_level_csv_folder);
    os.makedirs(patch_level_csv_folder); 
    
  if not os.path.exists(patch_level_plot_folder):
    print ('%s folder do not exist, then create it.' % patch_level_plot_folder);
    os.makedirs(patch_level_plot_folder);
    
  if not os.path.exists(image_level_csv_folder):
    print ('%s folder do not exist, then create it.' % image_level_csv_folder);
    os.makedirs(image_level_csv_folder);
    
  if not os.path.exists(image_level_plot_folder):
    print ('%s folder do not exist, then create it.' % image_level_plot_folder);
    os.makedirs(image_level_plot_folder);         
          
  #get image list from image path
  print (" --- get image list from image file path ---" );
  image_list=[];
  file_end="";
  IMAGE_FORMAT=IMAGE_FORMAT.replace("'", "");  
  file_end='.'+str(IMAGE_FORMAT);
  
  imagefile_list = [f for f in os.listdir(imagefiles_path) if f.endswith(file_end)] ;   
  for filename in imagefile_list:    
    case_id=filename.rstrip(file_end);    
    image_list.append(case_id);
  print ("total rows from image_list file is %d " % len(image_list) );   
  
  #read radiomics_feature_selected.txt file
  feature_array=[]; 
  radiomics_feature_selected = "radiomics_features_selected.txt"
  feature_selected_file_path = os.path.join(code_base, radiomics_feature_selected); 
  with open(feature_selected_file_path,) as f:
    reader = csv.reader(f, delimiter=',')
    my_list = list(reader);    
    for each_row in my_list:                      
      feature=str(each_row[0]);
      yes_no=str(each_row[1]);      
      if yes_no == "yes":   
        feature_array.append(feature); 
          
  #generate patch level radiomics feature csv files
  for case_id in image_list: 
    cmd = "python run_radiomics_all_patches.py  " + case_id ;     
    proc = subprocess.Popen(cmd, shell=True); 
    status = proc.wait() ;
 
  generate_Patch_Level_Radiomics_Feature_Plots(feature_array);              
  generate_image_level_radiomics_feature_plots(feature_array);  
  generate_image_level_radiomics_feature_csv_file(feature_array);               
  generate_image_level_nuclei_ratio_csv_file();     
                                             
  exit(); 
