# Robotic Hand Tracking

Real-time hand and finger joint tracking using OpenCV and MediaPipe, developed as the control interface for an Ecoflex-based pneumatic soft robotic hand.

## Features

- Real-time webcam hand tracking
- Detection and visualization of 21 hand landmarks
- Measurement of finger joint angles
- Normalized finger curl values from 0 to 1

## Installation

1. Clone the repository.
2. Create and activate a Python virtual environment.
3. Install the dependencies:

   ```bash
   pip install -r requirements.txt
   ```

4. Download the MediaPipe Hand Landmarker model:

   ```bash
   curl -L -o hand_landmarker.task https://storage.googleapis.com/mediapipe-models/hand_landmarker/hand_landmarker/float16/1/hand_landmarker.task
   ```

5. Run the program:

   ```bash
   python main.py
   ```

## Controls

Press `Q` to close the webcam window.

## Future Development

Integrate finger curl measurements with an ESP32-based pneumatic controller to operate individually actuated Ecoflex robotic fingers.
