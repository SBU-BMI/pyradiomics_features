import SimpleITK as sitk
from radiomics import featureextractor
import numpy as np
import skimage
import logging

class PyRadiomics_Features(object):
    def __init__(self, feature_setting_file):
        with open(feature_setting_file) as f:
            content = f.readlines()
        content = [x.strip() for x in content]

        settings = {}
        settings['binWidth'] = 25
        settings['resampledPixelSpacing'] = None  # [3,3,3] is an example for defining resampling (voxels with size 3x3x3mm)
        settings['interpolator'] = sitk.sitkBSpline

        self.extractor = featureextractor.RadiomicsFeatureExtractor(**settings)
        self.extractor.disableAllFeatures()
        self.extractor.enableFeaturesByName(firstorder=['Mean', 'Skewness'])

        feature_dict = {}
        for line in content:
            feature_type, feature_name = line.split(':')
            if feature_type in feature_dict:
                feature_dict[feature_type].append(feature_name)
            else:
                feature_dict[feature_type] = [feature_name];

        for feature_type, feature_name_list in feature_dict.items():
            if feature_type == 'firstorder':
                self.extractor.enableFeaturesByName(firstorder=feature_name_list);
            if feature_type == 'glcm':
                self.extractor.enableFeaturesByName(glcm=feature_name_list);
            if feature_type == 'gldm':
                self.extractor.enableFeaturesByName(gldm=feature_name_list);
            if feature_type == 'shape':
                self.extractor.enableFeaturesByName(shape=feature_name_list);
            if feature_type == 'glrlm':
                self.extractor.enableFeaturesByName(glrlm=feature_name_list);
            if feature_type == 'glszm':
                self.extractor.enableFeaturesByName(glszm=feature_name_list);
            if feature_type == 'ngtdm':
                self.extractor.enableFeaturesByName(ngtdm=feature_name_list);

        self.generate_feature_name_list();

    def get_feature_name_list(self):
        return self.feat_name_list;

    def generate_feature_name_list(self):
        # Create pseudo data
        data = np.zeros(shape=(4, 4), dtype=np.uint8);
        data[1:3, 1:3] = 1;
        data[1, 1] = 0;
        mask = np.zeros(shape=(4, 4), dtype=np.uint8);
        mask[1:3, 1:3] = 1;
        
        data = np.expand_dims(data, axis=2)
        mask = np.expand_dims(mask, axis=2)
        
        im = sitk.GetImageFromArray(data)  # ensure len(im_arr.shape) == 3
        ma = sitk.GetImageFromArray(mask)  # ensure len(ma_arr.shape) == 3
        im.SetSpacing((1, 1, 1))  # (x, y, z)
        ma.SetSpacing((1, 1, 1))  # (x, y, z)

        featureVector = self.extractor.execute(im, ma);

        self.feat_name_list = [];
        self.feat_name_index = {};
        feat_idx = 0;
        for rawfeatureName in featureVector.keys():
            if rawfeatureName.startswith('original_'):
                featureName = rawfeatureName.split('original_')[1]
                self.feat_name_list.append(featureName);
                self.feat_name_index[featureName] = feat_idx;
                feat_idx = feat_idx + 1;


    def cal_pyradiomics(self, image, mask):
        '''
        # Create simple data
        data = np.zeros(shape=(4, 4), dtype=np.uint8);
        data[1:3, 1:3] = 1;
        data[1, 1] = 0;
        # data[1,2] =0;
        mask = np.zeros(shape=(4, 4), dtype=np.uint8);
        mask[1:3, 1:3] = 1;
        '''

        #gray_image = skimage.color.rgb2gray(image)
        #bin_mask = (mask > 0).astype(np.uint8)
        gray_image = np.expand_dims(skimage.color.rgb2gray(image), axis=2)
        bin_mask = np.expand_dims((mask > 0).astype(np.uint8), axis=2)
        
        # Confirm mask is not all black
        if np.sum(bin_mask) == 0:
          res_features = ['None'] * len(self.feat_name_list)
          return res_features;
            
        im = sitk.GetImageFromArray(gray_image)  # ensure len(im_arr.shape) == 3
        ma = sitk.GetImageFromArray(bin_mask)  # ensure len(ma_arr.shape) == 3
        im.SetSpacing((1, 1, 1))  # (x, y, z)
        ma.SetSpacing((1, 1, 1))  # (x, y, z)

        # Calculating features
        res_features = [''] * len(self.feat_name_list)
        
        try:
            featureVector = self.extractor.execute(im, ma);

            for rawfeatureName in featureVector.keys():
                if not rawfeatureName.startswith('original_'):
                    continue;
                featureName = rawfeatureName.split('original_')[1]
                res_features[self.feat_name_index[featureName]] = featureVector[rawfeatureName];
                #print("Computed %s: %s" % (featureName, featureVector[featureName]))
                #res_features.append((featureName, featureVector[featureName]));
        except Exception as e:
            print ('error')
            print (e);
            logging.debug(e);
            for idx, featName in enumerate(self.feat_name_list):
                res_features[idx] = 'None';

        return res_features;
