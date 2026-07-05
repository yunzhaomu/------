'''剪切MOD13A1中对应城市的NDVI影像，然后计算均值并储存'''
from osgeo import gdal
from PIL import Image
from lxml import etree

import os
import ee
import requests
import numpy as np
import pandas as pd


class RS_Download(object):
    def __init__(self,years,city):
        '''返回指定年指定城市的NDVI均值'''
        self.years = years
        self.city = city

    def __call__(self):
        # 设置网络代理
        os.environ['HTTP_PROXY'] = 'http://127.0.0.1:1081'
        os.environ['HTTPS_PROXY'] = 'http://127.0.0.1:1081'
        # GEE初始化
        ee.Authenticate()
        ee.Initialize(project='ee-kk87232433')
        # 指定研究区域
        fc = ee.FeatureCollection('projects/ee-kk87232433/assets/city')
        cities = {city.replace('市',''):city for city in fc.aggregate_array('市').getInfo()}
        if self.city in cities.keys():
            roi = fc.filter(ee.Filter.eq('市',cities[self.city])).geometry()
        else:
            return [np.nan]*len(self.years)
        # 目前只看6月1日到9月30日的NDVI
        ndvis = []
        for year in self.years:
            imgs =  ee.ImageCollection('MODIS/061/MOD13A1')\
                    .filterDate(ee.Date.fromYMD(year, 6, 1), ee.Date.fromYMD(year, 9, 30))\
                    .filterBounds(roi)\
                    .select('NDVI')
            img = imgs.mean().multiply(0.0001).clip(roi)
            img_mean = img.reduceRegion(**{
                'reducer': ee.Reducer.mean(),
                'geometry': roi,
                'scale': 500,
                'maxPixels': 1e9
            })
            if img_mean.size().getInfo() > 0:
                ndvis.append(img_mean.values().getInfo()[0])
            else:
                ndvis.append(np.nan)
        return ndvis

def get_cities():
    # headers = {
    #     'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/90.0.4430.85 Safari/537.36 Edg/90.0.818.49'
    #     }
    # url = 'https://cyrdebr.sass.org.cn/_s20/2020/1120/c5517a99217/page.psp'
    # cities = []
    # page_text = requests.get(url=url,headers=headers).text
    # tree = etree.HTML(page_text)
    # tr_list = tree.xpath('//*[@id="d-container"]/div/div/div/div/div/div/table[1]/tbody/tr')
    # for tr in tr_list[1:]:
    #     for i in [1,4,7]:
    #         cities.append(tr.xpath(f'./td[{i}]/p/span')[0].text.strip())
    # cities = cities[:-1]
    cities = pd.read_excel('data/城市信息.xlsx')['城市'].tolist()
    return cities

def main(years,cities,output_path):
    if os.path.exists(output_path):
        ndvi = pd.read_excel(output_path,index_col=0)
    else:
        ndvi = pd.DataFrame()
    years = list(set(years)-set(ndvi.columns))
    for city in cities:
        if city in ndvi.index and len(years) == 0:
            continue
        rs = RS_Download(years,city)
        ndvi.loc[city,years] = rs()
        print(f'{city}已录完')
        ndvi.to_excel(output_path)


if __name__ == '__main__':
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    years = range(2000,2024)
    cities = get_cities()
    output_path = 'data/NDVI.xlsx'
    main(years,cities,output_path)
