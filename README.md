# Bioreactor Simulation for Cultivated Meat Production

This repository contains a **two-phase flow simulation of a rocking bioreactor**, implemented in the [Basilisk](http://basilisk.fr/) open-source CFD platform. This collaborative work is conducted by the Radu Cimpeanu Scientific Computing Group at the University of Warwick and the Harris Lab at Brown University. The preprint of this work is available at https://arxiv.org/.

## 📌 Features:
- ✅ liquid-gas dynamics with embedded solid boundaries inside a rocking bioreactor
- ✅ Tracer advection/diffusion for evaluating mixing efficiency
- ✅ Oxygen transport, including advection, diffusion and gas-to-liquid interfacial transfer
- ✅ Body acceleration in a non-inertial frame of reference, accounting for gravity, Coriolis, and centrifugal effects

---

## 🛠️ Installations

### 1. Requirements
- [Basilisk](http://basilisk.fr/) (compiled with `qcc`)
- C compiler
- Gnuplot / FFmpeg (for visualizing results)
```bash
sudo apt install gnuplot imagemagick ffmpeg graphviz valgrind gifsicle pstoedit
```

### 2. Clone the repo
```bash
git clone https://github.com/rcsc-group/BioReactor
cd DriverCodes
```
### 3. Install the Basilisk
The code relies on Basilisk to solve the Navier–Stokes equations. Please refer to the installation page for setup instructions: http://basilisk.fr/src/INSTALL

### 4. Copy and paste the header files and compile the main code
- Copy the four header files (draw3.h, henry_oxy2.h, utils2.h, and view3.h) into the src folder of the Basilisk source directory.
- Run the shell script using: sh BioReactor.sh
- The shell script runs the executable file: ./Bioreactor L_bio ANGLE RPM
  - L_bio: Reference length scale in meters (e.g., 0.25)
  - ANGLE: Rocking angle in degrees (e.g., 7)
  - RPM – Rocking frequency in RPM (e.g., 32.5)


## ⚙️ Key Simulations Configuration Options

Modify flags at the top of main.c to enable features:

- EMBED: Enable embedded boundary for solid geometry
- OXYGEN: Enable oxygen concentration simulation
  - OXYGEN_CIRCLE: Initial distribution (circle) of oxygen (if OXYGEN == 1)
#define OXYGEN_AIR       1   // Initial distribution (air side) of oxygen (if OXYGEN == 1)

// Mixing strategies for tracer release
#define TRACER           1   // Enable passive tracer simulation
#define HORIZONTAL_MIXL  0   // Initial distribution (left side) of tracer: Horizontal mixing (if TRACER == 1)
#define HORIZONTAL_MIXR  0   // Initial distribution (right side) of tracer: Horizontal mixing (if TRACER == 1)
#define VERTICAL_MIXUP   1   // Initial distribution (top side) of tracer: Vertical mixing (if TRACER == 1)
#define VERTICAL_MIXDOWN 0   // Initial distribution (bottom side) of tracer: Vertical mixing (if TRACER == 1)

// Other simulation options
#define ACCELERATION     1   // Enable acceleration (rocking motion)
#define AMR              0   // Enable adaptive mesh refinement
#define REMOVE_DROP      0   // Enable automatic droplet removal
#define CFL_COND         0   // Use custom CFL number
#define DUMP             0   // Save dump output
#define NORMCAL          1   // Calculate statistics (norms)
#define FIGURES          1   // Save figures
#define VIDEOS           1   // Save videos

// Output options
#define OUT_FILES         1   // Output full fields
#define OUT_SPECIFIC_TIME 0   // Output data at specific time ranges
#define OUT_INTERFACE     1   // Save interface geometry

---

## 📁 Folder Structure

```bash
.                        
├── main.c               # Main Basilisk simulation code
├── henry_oxy2.h         # Header for oxygen transport config
├── view3.h, utils2.h    # Visualization and utility functions (customized)
├── Data_all/            # Simulation output (velocity, tracer, oxygen, etc.)
├── Fig_vol/, Fig_tr/, ... # Saved PNG images for different fields
├── logstats.dat         # Performance and runtime log
├── normf.dat            # Velocity/vorticity norms over time
├── vol_frac_interf.dat  # Interface positions and volume fraction summary
├── tr_oxy.dat           # Integrated tracer and oxygen values
```

---

## 📊 Outputs

Generates:
- `Data_all/*.txt`: all field variables and interface geometries for the chosen simulation times
- `.dat` files of statistics (e.g., vorticity, velocity, volume fraction, etc.) and performance logs
- `*.mp4` videos (vorticity, tracer, oxygen, volume fraction)
- `*.png` figures (vorticity, tracer, oxygen, volume fraction)

The sample videos can be found 

---

## 📌 References

If you use this code for research or teaching, please cite Basilisk and include a reference to this repository.

- GitHub Repo: https://github.com/yourusername/bioreactor-basilisk
- Author: Minki Kim, Dan M. Harris, Radu Cimpeneau
- License: MIT License

---

## 🧑 Contributing

Feel free to:
- Fork this repo
- Open issues
- Submit pull requests