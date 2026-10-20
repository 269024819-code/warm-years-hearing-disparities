#!/usr/bin/env python3
"""Fig 3 v9.2：三列地图。修复：统计文字位置、橙蓝发散色、sqrt标度、标签避让、去南极"""
import os; os.makedirs("out", exist_ok=True)
import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, Normalize
from matplotlib import font_manager

for f in ["/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
          "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"]:
    try: font_manager.fontManager.addfont(f)
    except Exception: pass
plt.rcParams['font.family'] = 'Liberation Sans'

cty = pd.read_csv('out/fig3_country_values.csv')
cty['burden_rate_change'] = (cty['B_ll']/cty['fpop55_l'] - cty['B_ee']/cty['fpop55_e']) * 1e5
cty = cty.rename(columns={'iso3': 'ISO_A3'})

world = gpd.read_file('data/ne_110m_admin_0_countries.shp')
world.loc[world['ISO_A3'] == 'TWN', 'ISO_A3'] = 'CHN'
world = world.dissolve(by='ISO_A3').reset_index()
world = world[world['ISO_A3'] != 'ATA']  # 去南极
world = world.to_crs('+proj=eqearth')
m = world.merge(cty, on='ISO_A3', how='left')

cmap_seq  = LinearSegmentedColormap.from_list('seq', ['#FFF7EC', '#FDBE85', '#E08214', '#8C2D04'])
cmap_div  = LinearSegmentedColormap.from_list('div', ['#2166AC', '#92C5DE', '#F7F7F7', '#F4A582', '#B35806'])
cmap_soft = LinearSegmentedColormap.from_list('soft', ['#67A9CF', '#D1E5F0', '#F7F7F7', '#FDDBC7', '#EF8A62'])
cmap_seq.set_under('#B0C4DE')  # 暴露下降国（4个）

def ssqrt(x): return np.sign(x) * np.sqrt(np.abs(x))

fig, axes = plt.subplots(3, 1, figsize=(7.2, 9.2))

# ── 面板 a：暴露变化（sequential orange，负值浅蓝）─────────────────────────
ax = axes[0]
vals = m['exposure_change'] * 100
norm = Normalize(vmin=0, vmax=80)
m.plot(column=vals, ax=ax, cmap=cmap_seq, norm=norm, edgecolor='white', linewidths=0.15,
       missing_kwds={'color': '#D9D9D9'})
ax.set_axis_off()
ax.set_title('a   Change in hot-year share (percentage points, late − early)',
             fontsize=10.5, loc='left', pad=3, color='#1A202C', fontweight='bold')
sm = plt.cm.ScalarMappable(cmap=cmap_seq, norm=norm)
cb = fig.colorbar(sm, ax=ax, orientation='horizontal', fraction=0.032, pad=0.005, aspect=40)
cb.ax.tick_params(labelsize=8, colors='#4A5568'); cb.outline.set_visible(False)
cb.ax.set_title('blue = decrease (4 countries)', fontsize=7.5, color='#4A5568', pad=2)
ax.text(0.995, 0.10, '148 of 160 countries increased', transform=ax.transAxes,
        fontsize=9.5, ha='right', color='#1A202C', fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#CBD5E0', alpha=0.9))

# ── 面板 b：绝对负担变化（橙蓝发散，signed sqrt）────────────────────────────
ax = axes[1]
vals = m['abs_change']
vmax = np.nanmax(np.abs(vals))
norm = TwoSlopeNorm(vmin=-ssqrt(vmax), vcenter=0, vmax=ssqrt(vmax))
m_plot = m.assign(v=ssqrt(vals))
m_plot.plot(column='v', ax=ax, cmap=cmap_div, norm=norm, edgecolor='white', linewidths=0.15,
            missing_kwds={'color': '#D9D9D9'})
ax.set_axis_off()
ax.set_title('b   Change in absolute burden (thousand YLDs, signed square-root scale)',
             fontsize=10.5, loc='left', pad=3, color='#1A202C', fontweight='bold')
sm = plt.cm.ScalarMappable(cmap=cmap_div, norm=norm)
cb = fig.colorbar(sm, ax=ax, orientation='horizontal', fraction=0.032, pad=0.005, aspect=40)
ticks = [-100, -25, 0, 25, 100, 400, 800]
cb.set_ticks([ssqrt(t) for t in ticks]); cb.set_ticklabels([str(t) for t in ticks])
cb.ax.tick_params(labelsize=8, colors='#4A5568'); cb.outline.set_visible(False)
ax.text(0.995, 0.10, '156 of 160 countries increased', transform=ax.transAxes,
        fontsize=9.5, ha='right', color='#1A202C', fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#CBD5E0', alpha=0.9))
# 代表国家标注（手动偏移避免重叠）
offsets = {'CHN': (0, 14), 'USA': (0, 16), 'JPN': (26, 6), 'IND': (-4, -18), 'BRA': (0, -16)}
names = {'CHN': 'China', 'USA': 'United States', 'JPN': 'Japan', 'IND': 'India', 'BRA': 'Brazil'}
for iso, (dx, dy) in offsets.items():
    geo = m[m['ISO_A3'] == iso]
    if len(geo) == 0: continue
    c = geo.geometry.iloc[0].centroid
    ax.annotate(names[iso], xy=(c.x, c.y), xytext=(dx, dy), textcoords='offset points',
                fontsize=8.5, color='#1A202C', fontweight='bold', ha='center',
                bbox=dict(boxstyle='round,pad=0.18', fc='white', ec='#CBD5E0', alpha=0.85),
                arrowprops=dict(arrowstyle='-', color='#4A5568', lw=0.6))

# ── 面板 c：负担率变化（柔和发散，signed sqrt）──────────────────────────────
ax = axes[2]
vals = m['burden_rate_change']
vmax = np.nanmax(np.abs(vals))
norm = TwoSlopeNorm(vmin=-ssqrt(vmax), vcenter=0, vmax=ssqrt(vmax))
m_plot = m.assign(v=ssqrt(vals))
m_plot.plot(column='v', ax=ax, cmap=cmap_soft, norm=norm, edgecolor='white', linewidths=0.15,
            missing_kwds={'color': '#D9D9D9'})
ax.set_axis_off()
ax.set_title('c   Change in burden rate (YLDs per 100,000 women aged ≥55, signed square-root scale)',
             fontsize=10.5, loc='left', pad=3, color='#1A202C', fontweight='bold')
sm = plt.cm.ScalarMappable(cmap=cmap_soft, norm=norm)
cb = fig.colorbar(sm, ax=ax, orientation='horizontal', fraction=0.032, pad=0.005, aspect=40)
ticks = [-0.5, -0.25, 0, 0.25, 0.5, 1.0, 1.6]
cb.set_ticks([ssqrt(t) for t in ticks]); cb.set_ticklabels([str(t) for t in ticks])
cb.ax.tick_params(labelsize=8, colors='#4A5568'); cb.outline.set_visible(False)
ax.text(0.995, 0.10, '152 of 160 countries increased', transform=ax.transAxes,
        fontsize=9.5, ha='right', color='#1A202C', fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#CBD5E0', alpha=0.9))

plt.tight_layout(h_pad=0.6)
plt.savefig('Figure3_main.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('Figure3_main.pdf', bbox_inches='tight', facecolor='white')
print('Fig 3 v9.2 saved')
