
import cv2
import mediapipe as mp
import numpy as np
import time
from pathlib import Path

# ---------- Configuration ----------

MODEL_PATH = str(Path(__file__).with_name("hand_landmarker.task"))

BaseOptions = mp.tasks.BaseOptions
HandLandmarker = mp.tasks.vision.HandLandmarker
HandLandmarkerOptions = mp.tasks.vision.HandLandmarkerOptions
RunningMode = mp.tasks.vision.RunningMode

# MediaPipe hand skeleton connections
CONNECTIONS = [
    (0, 1), (1, 2), (2, 3), (3, 4),
    (0, 5), (5, 6), (6, 7), (7, 8), 
    (5, 9), (9, 10), (10, 11), (11, 12),
    (9, 13), (13, 14), (14, 15), (15, 16),
    (13, 17), (17, 18), (18, 19), (19, 20),
    (0, 17)
]

# Three landmarks used to calculate one joint angle
JOINTS = {
    "Thumb": (2, 3, 4),
    "Index": (5, 6, 7),
    "Middle": (9, 10, 11),
    "Ring": (13, 14, 15),
    "Pinky": (17, 18, 19)
}


def joint_angle(a, b, c):
    """
    Compute the angle ABC in degrees.
    a, b, c are 3D landmark coordinates.
    """
    ba = a - b
    bc = c - b

    denominator = np.linalg.norm(ba) * np.linalg.norm(bc)

    if denominator < 1e-8:
        return 180.0

    cosine = np.dot(ba, bc) / denominator
    cosine = np.clip(cosine, -1.0, 1.0)

    return np.degrees(np.arccos(cosine))

def normalize_curl(angle):
    return float(np.clip((180 - angle) / 90, 0, 1))

def draw_hand(frame, landmarks):
    """Draw connections and numbered landmarks."""
    h, w = frame.shape[:2]

    points = [
        (int(lm.x * w), int(lm.y * h))
        for lm in landmarks
    ]

    for start, end in CONNECTIONS:
        cv2.line(
            frame, points[start], points[end],
            (0, 255, 0), 2
        )

    for i, point in enumerate(points):
        cv2.circle(frame, point, 4, (0, 0, 255), -1)
        cv2.putText(
            frame, str(i),
            (point[0] + 5, point[1] - 5),
            cv2.FONT_HERSHEY_SIMPLEX,
            0.4, (255, 255, 255), 1
        )


def main():
    options = HandLandmarkerOptions(
        base_options=BaseOptions(model_asset_path=MODEL_PATH),
        running_mode=RunningMode.VIDEO,
        num_hands=1,
        min_hand_detection_confidence=0.6,
        min_hand_presence_confidence=0.6,
        min_tracking_confidence=0.6
    )

    cap = cv2.VideoCapture(0)

    if not cap.isOpened():
        raise RuntimeError("Cannot open webcam")

    # MediaPipe VIDEO mode requires increasing timestamps
    start_time = time.monotonic()
    last_timestamp = -1

    with HandLandmarker.create_from_options(options) as landmarker:
        while True:
            success, frame = cap.read()

            if not success:
                break

            # Mirror the image for intuitive movement
            frame = cv2.flip(frame, 1)

            rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
            image = mp.Image(
                image_format=mp.ImageFormat.SRGB,
                data=np.ascontiguousarray(rgb)
            )

            timestamp = int(
                (time.monotonic() - start_time) * 1000
            )
            timestamp = max(timestamp, last_timestamp + 1)
            last_timestamp = timestamp

            result = landmarker.detect_for_video(
                image, timestamp
            )

            if result.hand_landmarks:
                landmarks = result.hand_landmarks[0]
                draw_hand(frame, landmarks)

                # Prefer world landmarks for 3D joint angles
                coords = np.array([
                    [lm.x, lm.y, lm.z]
                    for lm in result.hand_world_landmarks[0]
                ])

                for i, (finger, (a, b, c)) in enumerate(JOINTS.items()):
                    angle = joint_angle(
                        coords[a], coords[b], coords[c]
                    )

                    curl = normalize_curl(angle)

                    # Display angle in degrees
                    cv2.putText(
                        frame,
                        f"{finger}: {angle:.1f} deg",
                        (20, 35 + i * 30),
                        cv2.FONT_HERSHEY_SIMPLEX,
                        0.7, (255, 255, 0), 2
                    )

                    # Display normalized curl separately
                    cv2.putText(
                        frame,
                        f"Curl: {curl:.2f}",
                        (260, 35 + i * 30),
                        cv2.FONT_HERSHEY_SIMPLEX,
                        0.7, (0, 255, 255), 2
                    )
            else:
                cv2.putText(
                    frame, "No hand detected",
                    (20, 35),
                    cv2.FONT_HERSHEY_SIMPLEX,
                    0.8, (0, 0, 255), 2
                )

            cv2.imshow("Robotic Hand Tracking", frame)

            # Press Q to quit
            if cv2.waitKey(1) & 0xFF == ord("q"):
                break
            

    cap.release()
    cv2.destroyAllWindows()


if __name__ == "__main__":
    main()
