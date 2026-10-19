# Soft robotic hand

Control pipeline for an Ecoflex pneumatic soft hand. Three modules share this repo:

| Module | Path | State |
| --- | --- | --- |
| Perception | [perception/](perception/) | Imported webcam tracker (OpenCV + MediaPipe) |
| Pneumatics | [pneumatics/](pneumatics/) | Command mapping, not implemented |
| Simulation | [simulation/](simulation/) | Physical playback, not implemented |

The tracker used to live at the repository root. Its program and original README moved with history into `perception/`. Root [requirements.txt](requirements.txt) is still `opencv-python`, `mediapipe`, and `numpy`.

## Run the tracker

From the repository root:

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
curl -L -o perception/hand_landmarker.task https://storage.googleapis.com/mediapipe-models/hand_landmarker/hand_landmarker/float16/1/hand_landmarker.task
python perception/main.py
```

Press `Q` to close the webcam window. The model file is gitignored and must sit beside `perception/main.py`. The original write-up is [perception/README.md](perception/README.md). Its shell commands assumed the script was at the root; use the commands above.

## For agents

Start at [AGENTS.md](AGENTS.md). Plans live in [plans/](plans/). Shared notes live in [context/](context/).

## Git

Local branch `dev/hand-pipeline`, remote `https://github.com/chengzansong/robotic_hand_tracking.git`, no upstream. Publishing and merging wait for an explicit request. See [context/git-workflow.md](context/git-workflow.md).
