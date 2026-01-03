# SLRSim
Simulates Laser Satellite Tracking
Date: 18 November 2025
Purpose: Code/Program Function Description
Purpose
Satellite Laser Ranging (SLR) is the process of using lasers to illuminate satellites in orbit and measure the reflected laser light.  The data gathered by this process informs scientists of the physical characteristics and behaviors of the Earth and the satellites themselves.  This program is a simulation of a Satellite Laser Ranging (SLR) station’s operation.  The code accomplishes  the simulation of the process of laser ranging and of tracking a real-world satellite from an actual SLR station.  The satellite tracked is LAGEOS-1 and the station is NASA’s Goddard Geophysical and Astronomical Observatory (GGAO) located in Prince Georges County, Maryland.  The program was written to accomplish two goals:  First, to learn orbital mechanics, math, and science needed to understand and overcome the challenges to building and maintaining a working SLR system.  Second, to gain understanding of how conditions such as lighting, terrestrial disruptions, atmospheric, meteorological, and orbital perturbations affect the performance of SLR stations.  

The core features of this program are as follows:

•	Uses real, orbital elements used to predict satellite orbits.

•	Employs the python Skyfield library with millimeter-level geographic and orbital accuracy.

•	Laser characteristics used in real world SLR tracking.

•	Models natural photon behavior

•	Displays four live panel performance indicators used in current SLR station operations.
  1.	Top (Laser Ranging Photon Returns): Two-way range vs time (km) with least squares prediction fit (red line).
  2.	Middle (Live Return Pulse): Live 1-second histogram of TOF residuals (nanoseconds) 
  3.	Bottom (Cumulative Return Pulse Distribution): Gaussian fit 
  4.	Text: Real-time status (return rate, pulses fired, photons detected)

Visual Elements Top Graph

•	Color graded photon return plots: Gray=Missed, Blue=Single, Green=Multi

•	Prediction model: Red straight line with linear prediction of satellite motion

•	Laser fire indicator above top graph


Scientific and Educational Value of the Program

•	Method SLR uses to achieve centimeter to millimeter orbital predictions and tracking.

•	Effects of environmental constraining conditions (light, weather, orbital perturbations).

•	Understanding the nature of collecting and measuring photons using quantum noise detectors.

•	Method of fitting models to short orbital passes of satellites.

Bottom Line This is not a toy — it is a physically correct, visually authentic recreation of what happens when a billion-dollar laser system points at a golf-ball-sized satellite 6,000 km away, firing 10 times per second, and counts individual photons at the quantum/photon level returning from space.
