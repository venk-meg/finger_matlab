import mediapipe as mp
from mediapipe.tasks import python
from mediapipe.tasks.python import vision
from mediapipe import solutions
from mediapipe.framework.formats import landmark_pb2
import numpy as np
import cv2
import csv
import os
from datetime import datetime
import pandas as pd


#drawing 
MARGIN = 10  # pixels
FONT_SIZE = 5
FONT_THICKNESS = 10
HANDEDNESS_TEXT_COLOR = (88, 205, 54) # vibrant green

def draw_landmarks_on_image(rgb_image, detection_result):
  hand_landmarks_list = detection_result.hand_landmarks
  handedness_list = detection_result.handedness
  annotated_image = np.copy(rgb_image)

  # Loop through the detected hands to visualize.
  for idx in range(len(hand_landmarks_list)):
    hand_landmarks = hand_landmarks_list[idx]
    handedness = handedness_list[idx]

    # Draw the hand landmarks.
    hand_landmarks_proto = landmark_pb2.NormalizedLandmarkList()
    hand_landmarks_proto.landmark.extend([
      landmark_pb2.NormalizedLandmark(x=landmark.x, y=landmark.y, z=landmark.z) for landmark in hand_landmarks
    ])
    solutions.drawing_utils.draw_landmarks(
      annotated_image,
      hand_landmarks_proto,
      solutions.hands.HAND_CONNECTIONS,
      solutions.drawing_styles.get_default_hand_landmarks_style(),
      solutions.drawing_styles.get_default_hand_connections_style())

    # Get the top left corner of the detected hand's bounding box.
    height, width, _ = annotated_image.shape
    x_coordinates = [landmark.x for landmark in hand_landmarks]
    y_coordinates = [landmark.y for landmark in hand_landmarks]
    text_x = int(min(x_coordinates) * width)
    text_y = int(min(y_coordinates) * height) - MARGIN

    # Draw handedness (left or right hand) on the image.
    cv2.putText(annotated_image, f"{handedness[0].category_name}",
                (text_x, text_y), cv2.FONT_HERSHEY_DUPLEX,
                FONT_SIZE, HANDEDNESS_TEXT_COLOR, FONT_THICKNESS, cv2.LINE_AA)

  return annotated_image

def showimage(image, required_width):
    """shows image"""
    if isinstance(image, str):
      img = cv2.imread(image) 
    else:
       img = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)

    #resizes image to given dimensions
    img_height, img_width = img.shape[:2]
    scaling_factor = required_width/img_width
    img = cv2.resize(img, (int(scaling_factor*img_width), int(scaling_factor*img_height)))

    #show image
    cv2.imshow('image',img)
    # cv2.waitKey(0)
    # cv2.destroyAllWindows()

    return img

def create_handlandmarker_obj():
    base_options = python.BaseOptions(model_asset_path='hand_landmarker.task')
    options = vision.HandLandmarkerOptions(base_options=base_options, num_hands=1)
    detector = vision.HandLandmarker.create_from_options(options)
    return detector

def detect_handlandmarks_from_frame(frame, detector):
    """Detect hand landmarks from a single video frame."""
    # Convert frame to MediaPipe image format
    rgb_frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
    mp_image = mp.Image(image_format=mp.ImageFormat.SRGB, data=rgb_frame)
    detection_result = detector.detect(mp_image)
    return detection_result

def writeheader_tocsv(filename):
  # Open the file in append mode
    with open(filename, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["frame_count", 
                          "knuckle_x", "knuckle_y", "knuckle_z", 
                          "PIP_x", "PIP_y", "PIP_z", 
                          "DIP_x", "DIP_y", "DIP_z", 
                          "tip_x", "tip_y", "tip_z" ])
   

def writefinger_tocsv(detection_result, filename, frame_count):
    # Extract landmarks for the first detected hand
    hand_landmarks_list = detection_result.hand_landmarks[0]
    index_landmarks_list = [hand_landmarks_list[i] for i in [5,6,7,8]]
    
    # Open the file in append mode
    with open(filename, mode="a", newline="") as file:
        writer = csv.writer(file)
      
        # Write landmark data to CSV   
        row = [frame_count]
        for landmark in index_landmarks_list:
          row.extend([landmark.x, landmark.y, landmark.z])
        writer.writerow(row)

def adjustcsv(filename):
  """
    Adjust the knuckle position in each row to be the origin (0, 0, 0).
    Save the adjusted file.
    """
  # Load the original CSV
  data = pd.read_csv(filename)
    
  # Adjust coordinates relative to the wrist
  for col in ['x', 'y', 'z']:
      knuckle_col = f'knuckle_{col}'
      for joint in ['PIP', 'DIP', 'tip']:
          joint_col = f'{joint}_{col}'
          data[joint_col] = data[joint_col] - data[knuckle_col]
      
      # Set wrist column to 0
      data[knuckle_col] = 0
    
  # Save to new CSV
  data.to_csv(f'{filename[:-4]}_adjusted.csv', index=False)


def loadandextract(videodir):
    """loads video, goes through all the frames and outputs required csv files and annotated video"""

    #YOU CAN CHANGE THESE VALUES
    # videodir = "indexflick_2.mp4" #test video directory
    filename = f"{videodir[:-4]}_index_finger_landmarks_test"

    #LOAD AND INITIALISE
    #loadvideo
    cap = cv2.VideoCapture(videodir)
    if not cap.isOpened():
        print("Error: Could not open video.")
        return
    #initialise csv
    filename = f"{filename}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
    writeheader_tocsv(filename)
    #create handlandmarker object
    detector = create_handlandmarker_obj()
    # Get video properties for saving output
    frame_width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    frame_height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    fps = int(cap.get(cv2.CAP_PROP_FPS))

    # Initialize VideoWriter for annotated frames
    output_video = cv2.VideoWriter(
        f"{filename[:-4]}_annotated.mp4",
        cv2.VideoWriter_fourcc(*'mp4v'),
        fps,
        (frame_width, frame_height)
    )

    #iterate through frames and extract handlandmarkers to csv
    frame_count = 1
    while cap.isOpened():
        ret, frame = cap.read()
        if not ret:
            break

        # Detect hand landmarks
        detection_result = detect_handlandmarks_from_frame(frame, detector)

        if detection_result and detection_result.hand_landmarks:
          writefinger_tocsv(detection_result, filename, frame_count)

          # Annotate the frame with landmarks
          annotated_frame = draw_landmarks_on_image(frame, detection_result)
          output_video.write(annotated_frame)  # Write annotated frame to output video

        frame_count += 1

        # # Display the video with landmarks (optional for debugging)
        # if detection_result:
        #   showimage(annotated_frame,300)
        # # Exit video on 'q' key press
        # if cv2.waitKey(1) & 0xFF == ord('q'):
        #     cv2.destroyAllWindows()
        #     break
        
    cap.release()
    output_video.release()

    adjustcsv(filename)

    print(f"Landmark positions saved to {filename} and {filename[:-4]}_adjusted.csv.")
    print(f"Annotated video saved as {filename[:-4]}_annotated.mp4.")









# execution:
videodirlist= ["earlier/og_videos/20241130_193019.mp4","earlier\og_videos\IMG_0266.MOV"]
for videodir in videodirlist:
   loadandextract(videodir)
