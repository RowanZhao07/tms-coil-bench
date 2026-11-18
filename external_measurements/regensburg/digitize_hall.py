"""Read Fig. 2 marker centres from the official Schecklmann et al. (2020) PNG.

Input is a published raster graph, not the laboratory's original samples.
Pixel calibration and marker radius were inspected on the 995 x 797 graphic.
"""
from pathlib import Path
import csv
import hashlib
import urllib.request
import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates

HERE = Path(__file__).resolve().parent
CACHE = HERE.parents[1] / "tmp" / "regensburg"
URL = ("https://media.springernature.com/full/springer-static/image/"
       "art%3A10.1038%2Fs41598-020-58034-2/MediaObjects/41598_2020_58034_Fig2_HTML.png")
SHA = "576d70880defdbffc9a3bdd213c4db10cae6e4fa9128b05bc7c66403a43e0635"
X0, X100 = 88.15, 971.55
Y0, Y07 = 708.35, 10.90

def centre_score(signal, x, y_candidates):
    # Marker rings are visible even with capped error bars and connecting lines.
    # Ignore sectors closest to the horizontal/vertical background strokes.
    alpha = np.arange(5, 360, 5) * np.pi / 180
    alpha = alpha[(abs(np.sin(alpha))>.24) & (abs(np.cos(alpha))>.24)]
    best = np.full(len(y_candidates),-np.inf)
    for radius in (6, 7, 8):
        def sample(r):
            yy = y_candidates[:,None] + r*np.sin(alpha)[None,:]
            xx = np.broadcast_to(x + r*np.cos(alpha),yy.shape)
            return map_coordinates(signal,[yy,xx],order=1,mode="nearest").mean(axis=1)
        best = np.maximum(best,sample(radius) - .45*sample(radius-3) - .40*sample(radius+3))
    return best

def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    src = CACHE / "fig2.png"
    if not src.is_file():
        urllib.request.urlretrieve(URL,src)
    assert hashlib.sha256(src.read_bytes()).hexdigest()==SHA
    rgb = np.asarray(Image.open(src).convert("RGB"),dtype=float)
    assert rgb.shape==(797,995,3)
    output = []
    for coil,c,last in [("C-B60",2,81),("D-B80",0,100)]:
        signal = np.clip(np.minimum(rgb[:,:,c]-rgb[:,:,2-c],rgb[:,:,c]-rgb[:,:,1])/85,0,1)
        for z in range(1,last+1):
            x = X0 + z*(X100-X0)/100
            strip = signal[:,round(x)-2:round(x)+3]
            yy,_ = np.where(strip[100:720,:]>.36)
            if not len(yy):
                raise RuntimeError(f"No coloured trace at {coil} {z} mm")
            initial = np.median(yy)+100
            # The colour median of the centre strip gives a neighbourhood; a
            # circular-marker search keeps this independent of the plotted error bars.
            cand = np.arange(max(103,initial-11),min(711,initial+11)+.01,.25)
            score = centre_score(signal,x,cand)
            centre = cand[np.argmax(score)]
            tesla = (Y0-centre)*.7/(Y0-Y07)
            output.append(dict(coil=coil,distance_mm=z,B_peak_T=round(tesla,6),
                               marker_y_px=round(centre,2),strip_median_y_px=round(initial,2)))
    with (HERE/"hall_fig2.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=output[0].keys());writer.writeheader();writer.writerows(output)
    print("Digitized",len(output),"markers")

if __name__=="__main__":
    main()
