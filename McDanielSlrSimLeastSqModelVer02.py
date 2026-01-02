# McDanielSlrSimLeastSqModelVer01.py
# Last Update 2025-11-18T13:22:00Z
# Author Email: mark.mcdaniel@markmcdanielmd.com
# Program Simulates SLR Tracking
# Developed For Education and Testing Purposes Only
# Selected NASA-GPT Assitance Including: Recommendations, Optimization and Testing
# Discription: Program Simulates Satellite Laser Ranging (SLR)

#-----------------------------------------------------------------------------
# Load Dependencies
#-----------------------------------------------------------------------------
import sys
import numpy as np
from skyfield.api import load, wgs84, EarthSatellite
from scipy.stats import poisson
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches
from PyQt5.QtWidgets import QApplication, QMainWindow, QScrollArea, QVBoxLayout, QWidget
from PyQt5.QtCore import Qt

#-----------------------------------------------------------------------------
# Sets Up Theme  = Dark Mode (LLM Assisted)
#-----------------------------------------------------------------------------
plt.style.use('dark_background')
plt.rcParams['axes.facecolor'] = '#0d1117'
plt.rcParams['figure.facecolor'] = '#0d1117'
plt.rcParams['text.color'] = '#c9d1d9'
plt.rcParams['axes.labelcolor'] = '#c9d1d9'
plt.rcParams['xtick.color'] = '#c9d1d9'
plt.rcParams['ytick.color'] = '#c9d1d9'
plt.rcParams['grid.color'] = '#303030'
plt.rcParams['grid.alpha'] = 0.4

#-----------------------------------------------------------------------------
# Load Satellite Two Line Elements In Standard TLE Format
#-----------------------------------------------------------------------------

tle_ln1 = "1 08820U 76039A 24217.55176019 .00000033 00000-0 00000-0 0 9990"
tle_ln2 = "2 08820 109.8511 71.2355 0044645 9.7861 9.1707 6.38664971 82982"
# EarthSatellite is a class from Skyfield that tranforms orbital data into an object
# load.timescale() is an object from Skyfield that creates a time-handling object
# that knows time formats leap seconds, UTC, etc and stores in a local cached file
satellite = EarthSatellite(tle_ln1, tle_ln2, "LAGEOS 1", load.timescale())
ts = load.timescale()
station = wgs84.latlon(39.0218, -76.8270, elevation_m=58.0)

#-----------------------------------------------------------------------------
# Setup  Generate timing conditons and station/satellite orientation
#-----------------------------------------------------------------------------
duration = 150
t_start = ts.utc(2024, 8, 5, 20, 0, 0)
# t-start timescale object created by load.timescale()
times = t_start + np.linspace(0, duration, 10000)
# Calc where satellite is from station point of view calculation by Skyfield
diff = satellite - station
# Point laswer at time t and calculate how far the satellite is and what angle.
topo = diff.at(times)
one_way_range_km = topo.distance().km
elevation = topo.altaz()[0].degrees
visible = elevation > 20
if not np.any(visible):
    raise RuntimeError("No visible pass")

t_vis = times[visible]
one_way_vis = one_way_range_km[visible]
sec_vis = np.array([(t - t_start) for t in t_vis])
interp_range = interp1d(sec_vis, one_way_vis, kind='cubic', bounds_error=False, fill_value="extrapolate")

#-----------------------------------------------------------------------------
# Laser & Noise - Atmospheric and Light Conditions
#-----------------------------------------------------------------------------

# Speed of light constant
c = 299_792_458
# Laser characteristics
pulse_rate = 10
n_pulses = int(pulse_rate * duration)
np.random.seed(42)
# Fire times elapsed  (in seconds)
fire_sec = np.sort(np.random.uniform(0, duration, n_pulses))  
one_way_fire = interp_range(fire_sec)
true_tof = 2 * one_way_fire * 1000 / c
measured_tof = true_tof + (0.020 * 2 / c) * np.random.randn(n_pulses) + 50e-12 * np.random.randn(n_pulses)
two_way_range_km = (measured_tof * c) / 1000

mean_photons_per_pulse = 1.15
# Adjust  lighting conditions, meterologic & astmospheric conditions
background_rate = 800_000
gate_width = 100e-9
expected_bg = background_rate * gate_width

#-----------------------------------------------------------------------------
# Animation state
#-----------------------------------------------------------------------------
current_pulse = 0
detected_photons = []
pulse_photon_counts = np.zeros(n_pulses, dtype=int)
def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

#-----------------------------------------------------------------------------
# Main Window (note style and window properties with LLM assist)
#-----------------------------------------------------------------------------

class ScrollWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LAGEOS-1 SLR - NASA GGAO - REAL SLR STYLE")
        self.setGeometry(100, 100, 950, 1500)
        fig = plt.figure(figsize=(9, 30))
        gs = fig.add_gridspec(4, 1, hspace=0.9)

        ax_range = fig.add_subplot(gs[0])
        global scat_range, line_fit, pulse_indicator
        scat_range = ax_range.scatter([], [], s=30, c='gray', edgecolor='#58a6ff', linewidth=0.8)
        line_fit, = ax_range.plot([], [], '#ff4444', lw=3, label='Linear Fit')
        ax_range.set_xlim(0, duration)
        ax_range.set_ylim(5000, 20000)
        ax_range.set_xlabel('Time (s)')
        ax_range.set_ylabel('Two-Way Range (km)')
        ax_range.set_title('LAGEOS-1 Daylight Pass — Clear Sky', color='white', fontsize=14)
        ax_range.grid(True, alpha=0.4)
        ax_range.legend(fontsize=10, loc='upper right')
        pulse_indicator = ax_range.add_patch(patches.Rectangle((0.02, 1.08), 0.06, 0.04,
        facecolor='red', edgecolor='white', transform=ax_range.transAxes, clip_on=False))

        ax_live = fig.add_subplot(gs[1])
        global live_bins, live_bars
        live_bins = np.linspace(-60, 60, 40)
        live_bars = ax_live.bar(live_bins[:-1], np.zeros(39), width=np.diff(live_bins)[0],
        color='#58a6ff', edgecolor='#1f6feb', alpha=0.9)
        ax_live.set_xlim(-60, 60); ax_live.set_ylim(0, 15)
        ax_live.set_title('Live Return Pulse (last 1 s)', color='white')

        ax_cum = fig.add_subplot(gs[2])
        global cum_bars, gauss_line
        cum_bars = ax_cum.bar(live_bins[:-1], np.zeros(39), width=np.diff(live_bins)[0],
        color='#bb9cf8', edgecolor='#8b5cf6')
        gauss_line, = ax_cum.plot([], [], '#ff6b6b', lw=2.5)
        ax_cum.set_xlim(-60, 60); ax_cum.set_ylim(0, 100)
        ax_cum.set_title('Cumulative Return Pulse Distribution', color='white')

        ax_status = fig.add_subplot(gs[3])
        ax_status.axis('off')
        global status_text
        status_text = ax_status.text(0.01, 0.98, '', va='top', ha='left',
        fontsize=11, fontfamily='monospace', color='#58a6ff')

        canvas = FigureCanvas(fig)
        scroll = QScrollArea()
        scroll.setWidget(canvas)
        scroll.setWidgetResizable(True)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.addWidget(scroll)
        self.setCentralWidget(container)

        global ani
        ani = FuncAnimation(fig, update, frames=int(duration*30), interval=33, blit=False, repeat=False)
        fig.suptitle('LAGEOS-1 — NASA GGAO — REAL SLR CONSOLE', fontsize=16, fontweight='bold', color='#58a6ff')
        canvas.draw_idle()

#-----------------------------------------------------------------------------
# Update function — Linear fit + all  globals declared
#-----------------------------------------------------------------------------

def update(frame):
    global current_pulse, scat_range, line_fit, pulse_indicator
    global live_bars, cum_bars, gauss_line, status_text

    t_now = frame / 30.0
    t_current = t_start + t_now
    elev = (satellite - station).at(t_current).altaz()[0].degrees
    tracking = elev > 20

    if tracking:
        new = [i for i in range(current_pulse, n_pulses) if fire_sec[i] <= t_now]
        for i in new:
            # poisson distro for photon arrival
            n_phot = poisson(mean_photons_per_pulse + expected_bg).rvs()
            pulse_photon_counts[i] = n_phot
            for _ in range(n_phot):
                off = np.random.uniform(-50, 50)
                detected_photons.append((t_now, off))
            current_pulse = i + 1

    if tracking:
        past = np.arange(current_pulse)[fire_sec[:current_pulse] <= t_now]
        if len(past) > 0:
            t_plot = fire_sec[past]
            r_plot = two_way_range_km[past]
            colors = ['#303030' if pulse_photon_counts[i]==0 else '#58a6ff' if pulse_photon_counts[i]==1 else '#8ac926' for i in past]
            scat_range.set_offsets(np.c_[t_plot, r_plot])
            scat_range.set_color(colors)

            valid = past[pulse_photon_counts[past] > 0]
            if len(valid) > 10:
                coeffs = np.polyfit(fire_sec[valid], two_way_range_km[valid], 1)
                t_fit = np.linspace(t_plot.min(), t_plot.max(), 500)
                r_fit = np.polyval(coeffs, t_fit)
                line_fit.set_data(t_fit, r_fit)
            else:
                line_fit.set_data([], [])

        recent = [off for t, off in detected_photons if t_now - 1 <= t <= t_now]
        if recent:
            counts, _ = np.histogram(recent, bins=live_bins)
            for bar, h in zip(live_bars, counts):
                bar.set_height(h)

        all_off = [off for t, off in detected_photons if t <= t_now]
        if all_off:
            counts, _ = np.histogram(all_off, bins=live_bins)
            for bar, h in zip(cum_bars, counts):
                bar.set_height(h)

    pulse_indicator.set_facecolor('#ff4444' if tracking and any(abs(fire_sec[i]-t_now)<0.05 for i in range(current_pulse)) else '#303030')

    hits = sum(1 for i in range(current_pulse) if pulse_photon_counts[i]>0 and fire_sec[i] >= t_now-10)
    total = sum(1 for i in range(current_pulse) if fire_sec[i] >= t_now-10)
    rate = (hits / max(total, 1)) * 100

    status_text.set_text(f"DAYLIGHT PASS — Clear Sky\n"
                         f"RETURN RATE (10s): {rate:5.1f}%\n"
                         f"Time: {t_now:5.1f}s | NASA GGAO\n"
                         f"Pulses: {current_pulse:4d} | Photons: {len(detected_photons):5d}")

    plt.gcf().canvas.draw_idle()
    return (scat_range, line_fit, pulse_indicator, *live_bars, *cum_bars, gauss_line, status_text)

#-----------------------------------------------------------------------------
# Execution Block
#-----------------------------------------------------------------------------

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ScrollWindow()
    window.show()
    sys.exit(app.exec_())
