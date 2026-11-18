"""Extract Fig. 2A (observations only) and the two 140-mm text endpoints."""
import csv
import hashlib
from pathlib import Path
import urllib.request
import numpy as np
from PIL import Image

HERE=Path(__file__).resolve().parent
CACHE=HERE.parents[1]/'tmp'/'lyon'
URL='https://pmc-oa-opendata.s3.amazonaws.com/PMC10876241.1/painreports-9-e1134-g002.jpg'
SHA='4399895bbe8b2f4c37c2750c56899c65d14b1f91f26fce793e17cb92ce9b46aa'

def main():
    CACHE.mkdir(parents=True,exist_ok=True)
    source=CACHE/'fig2.jpg'
    if not source.is_file():
        urllib.request.urlretrieve(URL,source)
    assert hashlib.sha256(source.read_bytes()).hexdigest()==SHA
    img=np.asarray(Image.open(source).convert('RGB'),dtype=float)
    assert img.shape==(577,578,3)
    result=[]
    for coil,channel,contact,far in [('Cool-B65',2,.88,.0057),('Cool-D-B80',0,.77,.00993)]:
        strength=img[:,:,channel]-img[:,:,2-channel]
        for distance in range(0,61,5):
            x=int(round(94.1+(549.0-94.1)*distance/60))
            strip=strength[20:256,x-2:x+3].copy()
            if distance>=30: strip[:65,:]=0 # published legend
            iy,ix=np.where(strip>42)
            weights=strip[iy,ix]
            assert len(iy)>=7
            # In the marker's core, weigh saturated samples over the dotted line.
            y=20+np.average(iy,weights=weights)
            raw=(259.2-y)*.9/(259.2-23.5)
            value=contact if distance==0 else raw
            result.append([coil,distance,round(value,7),'text_and_fig2a' if distance==0 else 'fig2a',round(raw,7)])
        result.append([coil,140,far,'paper_text',''])
    with (HERE/'air_observations.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['coil','distance_mm','B_T','source','figure_reading_T']);w.writerows(result)

if __name__=='__main__':
    main()
