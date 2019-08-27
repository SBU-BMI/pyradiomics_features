# pyradiomics_features
   A docker version of pyradiomics feature pipeline, which generates patch level pyradiomics feature dataset and image level
aggregational result. It also create histogram plot and percentile plot for each pyradiomics feature; 

  This develop branch code base uses one docker container.

Step 1: Clone source code from github:
```
git clone -b develop https://github.com/SBU-BMI/pyradiomics_features.git
```

Step 2: CD to ./pyradiomics_features/ folder

Step 3: Edit pyradiomics_features.yml of parameters IMAGE_FORMAT, PATCH_SIZE and HAS_TUMOR_REGION; The value of IMAGE_FORMAT could be svs and tif; The value of HAS_TUMOR_REGION could be yes and no; 

environment:

      - IMAGE_FORMAT=svs
      - PATCH_SIZE=500 
      - HAS_TUMOR_REGION=yes
      
Step 4: Build docker image in ./pyradiomics_features/ folder 
```
 docker-compose -f pyradiomics_features.yml build
```
Step 5: Run docker container in ./pyradiomics_features/ folder
  ```
  docker-compose -f pyradiomics_features.yml up -d
  ```
Step 6: Copy image svs or tif file to ./pyradiomics_features/image_files/ folder;

Step 7: Copy  segmentation results to ./pyradiomics_features/segment_results/ folder; First create subfolder under  /pyradiomics_features/segment_results folder using image file name without extension; for example, 
  if image name is my_image_name.svs, then subfolder name should be my_image_name.
  Create subfolder for each image, then copy segment results (json and csv pair files) to each subfolder.

Step 8: if applicable, copy tumor region json file to ./pyradiomics_features/tumor_region/ folder;

Step 9: Execute docker container
```
docker exec radiomics-core  /app/computing_patch_level_radiomics_features.sh
```
Step 10; verify results. All results are in ./pyradiomics_features/output/ folder. There should be four subfolders available after complete running, named "patch_level_csv", "patch_level_plot","image_level_csv","image_level_plot" .etc.

Step 11: Stop docker container  in ./pyradiomics_features/ folder
```
docker-compose -f pyradiomics_features.yml down
```

Step 12: Remove docker images if necessary
```
docker rmi pyradiomics_features_radiomics-core
```

Step 13: Remove all unused volumes if necessary
```
docker volume prune
```
