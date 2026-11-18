import serial
import time

def send_info(arduino, values):
    """
    Send five finger curvature values to Arduino.
    """
    message = ",".join(
        f"{values[finger]:.3f}"
        for finger in ["Thumb", "Index", "Middle", "Ring", "Pinky"]
    )

    arduino.write((message + "\n").encode())

    response = arduino.readline().decode().strip()
    if response:
        print(response)