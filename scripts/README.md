# Scripts

Operator scripts will live here. None are checked in yet.

The tracker still runs as a module entry point. From the repository root, after `pip install -r requirements.txt`:

```bash
curl -L -o perception/hand_landmarker.task https://storage.googleapis.com/mediapipe-models/hand_landmarker/hand_landmarker/float16/1/hand_landmarker.task
python perception/main.py
```

`hand_landmarker.task` is gitignored. It has to sit next to `perception/main.py` because the script loads `Path(__file__).with_name("hand_landmarker.task")`.
