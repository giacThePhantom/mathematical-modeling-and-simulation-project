import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tqdm import tqdm

# argv[1] e argv[2] path ai dati .csv
# argv[3] cartella per salvataggio immagini
# opzionali argv[4] e argv[5] per settaggio bins e numero canali totali

rtc = pd.read_csv(os.path.join("..", "data", sys.argv[1]), header=None)
gill = pd.read_csv(os.path.join("..", "data", sys.argv[2]), header=None)

rtc.rename(columns={0: "V", 1: "K", 2: "Ca"}, inplace=True)
gill.rename(columns={0: "V", 1: "K", 2: "Ca"}, inplace=True)

try: 
    bins = int(sys.argv[4])
except:
    bins = 1000

try: 
    n_tot = int(sys.argv[5]) + 1
except:
    n_tot = 101


name_to_df = {
    "RTC": rtc,
    "GILLESPIE": gill
}

interval = {
    "RTC": (max(rtc["V"]) - min(rtc["V"]))/bins,
    "GILLESPIE": (max(gill["V"]) - min(gill["V"]))/bins
}

for key in name_to_df:
    for channel in ["K", "Ca"]:
        
        print(f"Generating plot for {channel} channels, {key} simulation...")
        
        start = min(rtc["V"])
        end = min(rtc["V"]) + interval[key]
        heatmap = np.zeros((n_tot, bins))

        for i in tqdm(range(bins)):
            filtered_v = name_to_df[key].loc[(name_to_df[key]["V"] >= start) & (name_to_df[key]["V"] < end)]
            for y in range(n_tot):
                heatmap[y, i] += len(filtered_v.loc[(filtered_v[channel] == y)])
            start = end
            end = end + interval[key]

        x_axis = [x for x in range(0, bins, int(bins/20))]
        x_labels = [round(min(rtc["V"]) + x*interval[key], 2) for x in x_axis] 

        plt.figure(figsize = (30, 15))
        plt.title(f"{key} - {channel} CHANNELS")
        plt.imshow(heatmap, cmap = "inferno")
        plt.xticks(x_axis, x_labels)
        plt.gca().invert_yaxis()
        img_folder = os.path.join("..", "images", sys.argv[3])
        if not os.path.exists(img_folder):
            os.makedirs(img_folder)
        plt.savefig(os.path.join(img_folder, f"{channel}_{key}_.png"), dpi = 300)