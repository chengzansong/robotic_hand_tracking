# Simulation

Plays pneumatic commands against a physical model and reports the predicted hand state. Not implemented.

Commands arrive from [pneumatics/](../pneumatics/). The model will live in [models/](models/). The simulator choice is unspecified. If a later run needs the GPU, use CUDA 12.8 or newer on the RTX 5080 Mobile, and keep a single heavy GPU job at a time.

See [contracts/README.md](../contracts/README.md).
