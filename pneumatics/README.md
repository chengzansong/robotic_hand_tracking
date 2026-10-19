# Pneumatics

Maps perception output to commands for the pneumatic hand. Not implemented.

Input, once defined, is the per-finger curl (0 to 1) and joint angle from [perception/main.py](../perception/main.py). Output is an actuator command per finger. The command shape is unspecified. The upstream tracker notes an ESP32 controller for individually actuated Ecoflex fingers; that hardware note is not a dependency yet.

See [contracts/README.md](../contracts/README.md).
