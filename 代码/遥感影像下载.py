'''下载遥感影像到本地，并将土地分类原始数据转为RGB数据'''
from osgeo import gdal
from PIL import Image
from config import decode_segmap
from config import rs_classes

import os
import ee
import geemap
import argparse


class RS_Download(object):
    def __init__(self,year,download_path,file_name):
        assert file_name in rs_classes, f'目前支持下载{rs_classes}'
        self.year = year
        self.download_path = download_path
        self.file_name = file_name
        directory = os.path.join(self.download_path,self.file_name)
        if not os.path.exists(directory):
            os.makedirs(directory)
        self.file_path = os.path.join(directory,self.file_name+f'_{self.year}.tif')

    @staticmethod
    def return_properties(file_name,property):
        '''各类遥感影像对应的GEE源和相应信息，后续更改或添加遥感影像请维护这个函数'''
        source_mapping = {
            'land_labels':"projects/lulc-datase/assets/LULC_HuangXin",
            'ET':"MODIS/061/MOD16A2GF",
            'NPP':"MODIS/061/MOD17A3HGF",
            'rain':"UCSB-CHG/CHIRPS/DAILY"
        }
        scale_mapping = {
            'land_labels':30,
            'ET':500,
            'NPP':500,
            'rain':500
        }
        dtype_mapping = {
            'land_labels':"uint8",
            'ET':"float32",
            'NPP':"float32",
            'rain':"float32"
        }
        if property == 'source':
            return source_mapping[file_name]
        elif property == 'scale':
            return scale_mapping[file_name]
        elif property == 'dtype':
            return dtype_mapping[file_name]

    def download(self):
        # 下载所需属性
        rs_sorce = self.return_properties(self.file_name,property='source') # 下载源
        scale = self.return_properties(self.file_name,property='scale') # 下载源
        img_dtype = self.return_properties(self.file_name,property='dtype') # 下载源
        # 设置网络代理
        os.environ['HTTP_PROXY'] = 'http://127.0.0.1:1081'
        os.environ['HTTPS_PROXY'] = 'http://127.0.0.1:1081'
        # GEE初始化
        ee.Authenticate()
        ee.Initialize(project='ee-kk87232433')
        # 指定研究区域
        roi = ee.FeatureCollection('projects/ee-kk87232433/assets/city').filter(ee.Filter.eq('市','鄂尔多斯市')).geometry()
        # 不同遥感影像进行不同预处理
        if self.file_name == 'land_labels':
            img = ee.Image(rs_sorce+f'/CLCD_v01_{self.year}').clip(roi)
        elif self.file_name == 'ET':
            imgs =  ee.ImageCollection(rs_sorce)\
                    .filterDate(ee.Date.fromYMD(self.year, 1, 1), ee.Date.fromYMD(self.year+1, 1, 1))\
                    .filterBounds(roi)
            img = imgs.select('ET').sum().multiply(0.1).clip(roi)
        elif self.file_name == 'NPP':
            imgs =  ee.ImageCollection(rs_sorce)\
                    .filterDate(ee.Date.fromYMD(self.year, 1, 1), ee.Date.fromYMD(self.year+1, 1, 1))\
                    .filterBounds(roi)
            img = imgs.select('Npp').sum().multiply(0.0001).multiply(10).clip(roi)
        elif self.file_name == 'rain':
            imgs =  ee.ImageCollection(rs_sorce)\
                    .filterDate(ee.Date.fromYMD(self.year, 1, 1), ee.Date.fromYMD(self.year+1, 1, 1))\
                    .filterBounds(roi)
            img = imgs.select('precipitation').sum().clip(roi)
        # 下载
        geemap.download_ee_image(
            image=img,
            filename=self.file_path,
            region=roi.bounds(),
            scale=scale,
            dtype=img_dtype,
            crs="EPSG:4326" # 使用WGS84坐标系，一般存数据所用；EPSG:3857是各大互联网地图所用
        )
        return

    def convert_rgb(self):
        rgb_path = self.file_path.replace('_labels_','_rgb_')
        classMap = gdal.Open(self.file_path).ReadAsArray().astype('uint8') # shape:(H*W)
        classMap = decode_segmap(classMap)
        classMap = Image.fromarray(classMap, 'RGB') # 太大了，55419 x 36202，需要等比例缩小
        w,h = classMap.size
        classMap = classMap.resize((w//10,h//10))
        classMap.save(rgb_path)
        print(f'Image:Land_Classification_{self.year}, done.')
        return

    def __call__(self):
        self.download()
        if self.file_name == 'land_labels':
            self.convert_rgb()

def main(args):
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    for year in args.years:
        rs = RS_Download(year,args.rs_folder,args.rs_class)
        rs()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='download RS images')
    parser.add_argument('--years',
                        default=[2023],
                        help='指定年份的鄂尔多斯土地分类数据')
    parser.add_argument('--rs-folder',
                        default='遥感数据',
                        help='遥感影像存放目录路径')
    parser.add_argument('--rs-class',
                        default='rain',
                        help='要下载的遥感影像是什么')
    args = parser.parse_args()

    main(args)