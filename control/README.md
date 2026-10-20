# Control

Turns perception output into pump commands for the pneumatic hand. Not implemented yet.

- **Input:** the per-finger curl (0 to 1) and joint angle from [perception/main.py](../perception/main.py).
- **Output:** one pump command per finger. Tentatively, this is a PWM duty from 0 to 1 for each finger's DS-370-12 pump, because v1 has no valves and deflates by passive leakage. Valve states join the command once valves are added.

The hardware is planned as an ESP32 driving one pump per Ecoflex finger, through a MOSFET with a flyback diode.

The exact formats are in [CONTRACTS.md](../CONTRACTS.md).
