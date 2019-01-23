from urllib.parse import urlparse
import urllib.request
import re
import requests as req
from PIL import Image
from io import BytesIO
import pandas as pd 


kegg_out = pd.read_table('/home/zluna/Work/kegg/2.xls')



headers = {'User-Agent' : 'Mozilla/5.0 (Linux;Android 6.0;Nexus 5 Build/MRA58N) AppleWebKit/537.36 (KHTML,like Gecko) Chrome/58.0.3029.96 Mobile Safari/537.36'}
url = "https://www.kegg.jp/pathway/map00361+K00462%09red+K01061%09red+K01563%09red+K03119%09red+K06912%09red+K10621%09red+K14583%09red+K14751%09red+K16049%09red"


#list(kegg_out.Pathway)[0].split(':',1)[1].strip(']').lstrip('ko')
#str(kegg_out.KOID[0]).replace(', ', '%09purple+') + '%09purple'
for i in range(len(kegg_out)):
    path_name = list(kegg_out.Pathway)[i].split(':',1)[1].strip(']').lstrip('ko')
    path_gene = str(kegg_out.KOID[i]).replace(', ', '%09purple+') + '%09purple'
    url = "https://www.kegg.jp/pathway/map" + path_name + '+' + path_gene
    response = req.get(url, headers=headers)
    content = response.text
    img_url = re.search(r'<img src=.*pathwayimage.*>', str(content)).group()
    img_url.split('"', 3)[1]
    img_url = "https://www.kegg.jp" + img_url.split('"', 3)[1]
    img = req.get(img_url)
    #print(img)
    img = Image.open(BytesIO(img.content))
    img.save(fp = "Path:" + path_name + '.png')
