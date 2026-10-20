# Stage 2: whole hand (not started)

Planned as a MuJoCo reduced-order hand. Each finger is a chain of capsules joined by stiff, damped hinges, and chamber pressure becomes joint torque. The hinge stiffness, damping, and pressure-to-torque map get fitted to the Stage 1 FEM sweeps and to single-finger bench data.

This stage will read designs from [../models/](../models/) and materials from [../materials/](../materials/), plus shared code from [../core/](../core/). It may read Stage 1 *results* (calibration curves in `simulation/outputs/`), but must not import Stage 1 code.

Work starts after Stage 1 has settled the finger design and material.
