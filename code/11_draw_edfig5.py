#!/usr/bin/env python3
"""ED Fig 5：160国年均温趋势地图（2003-2019线性趋势，°C/十年）"""
import os; os.makedirs("out", exist_ok=True)
import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib import font_manager

for f in ["/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
          "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"]:
    try: font_manager.fontManager.addfont(f)
    except Exception: pass
plt.rcParams['font.family'] = 'Liberation Sans'

tr = pd.read_csv('out/table_climate_trends.csv').rename(columns={'iso3': 'ISO_A3'})
tr['trend_decade'] = tr['slope'] * 10

world = gpd.read_file('data/ne_110m_admin_0_countries.shp')
world.loc[world['ISO_A3'] == 'TWN', 'ISO_A3'] = 'CHN'
world = world.dissolve(by='ISO_A3').reset_index()
world = world[world['ISO_A3'] != 'ATA'].to_crs('+proj=eqearth')
m = world.merge(tr, on='ISO_A3', how='left')

cmap = LinearSegmentedColormap.from_list('div', ['#2166AC', '#92C5DE', '#F7F7F7', '#F4A582', '#B35806'])
vmax = float(np.nanmax(np.abs(m['trend_decade'])))
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

fig, ax = plt.subplots(figsize=(7.2, 3.6))
m.plot(column='trend_decade', ax=ax, cmap=cmap, norm=norm, edgecolor='white', linewidths=0.15,
       missing_kwds={'color': '#D9D9D9'})
ax.set_axis_off()
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cb = fig.colorbar(sm, ax=ax, orientation='horizontal', fraction=0.045, pad=0.02, aspect=42)
cb.set_label('Linear trend in annual mean temperature, 2003-2019 (°C per decade)', fontsize=9, color='#2D3748')
cb.ax.tick_params(labelsize=8, colors='#4A5568'); cb.outline.set_visible(False)
ax.text(0.995, 0.03, '152 of 160 countries positive (93 at P < 0.05); none significantly negative',
        transform=ax.transAxes, fontsize=9.5, ha='right', color='#1A202C', fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#CBD5E0', alpha=0.9))
plt.tight_layout()
plt.savefig('ED_Fig5_climate_trends.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('ED_Fig5_climate_trends.pdf', bbox_inches='tight', facecolor='white')
print('ED Fig 5 saved, vmax =', round(vmax, 2))
