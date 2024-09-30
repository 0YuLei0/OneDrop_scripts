###### ###### ###### ###### ###### ###### ###### ###### ###### ######
import os
import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fastsam import FastSAM, FastSAMPrompt
from tqdm import tqdm
import torch

# Set device for model inference
DEVICE = 'cpu'
# Initialize the FastSAM model
model = FastSAM('/import/home/lyuah/od_image/model/FastSAM-x.pt') ## accurate
#model = FastSAM('/import/home/lyuah/od_image/SAM_model/FastSAM-s.pt')
# Function to calculate CV
def calculate_cv(x):
    return x.std() / x.mean() * 100 if x.mean() != 0 else 0

def _format_results(result, filter=0):
    annotations = []
    n = len(result.masks.data)
    for i in range(n):
        annotation = {}
        mask = result.masks.data[i] == 1.0
        if torch.sum(mask) < filter:
            continue
        annotation['id'] = i
        annotation['segmentation'] = mask.cpu().numpy()
        annotation['bbox'] = result.boxes.data[i]
        annotation['score'] = result.boxes.conf[i]
        annotation['area'] = annotation['segmentation'].sum()
        annotations.append(annotation)
    return annotations

def process_image(image_path):
    results = model(image_path, device=DEVICE, retina_masks=True, imgsz=(1024, 1280), conf=0.75, iou=0.5)
    formatted_results = _format_results(results[0])
    return formatted_results

def filter_contours(masks, image_shape):
    valid_contours = []
    height, width, _ = image_shape
    edge_threshold = 25
    diameter_threshold = 100  # Minimum diameter threshold
    for mask in masks:
        segmentation = mask['segmentation']
        contours, _ = cv2.findContours(segmentation.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        for contour in contours:
            M = cv2.moments(contour)
            if M["m00"] == 0:
                continue
            area = cv2.contourArea(contour)
            perimeter = cv2.arcLength(contour, True)
            circularity = (4 * np.pi * area) / (perimeter ** 2) if perimeter != 0 else 0
            # Calculate diameter from area
            diameter = np.sqrt(4 * area / np.pi)
            if diameter >= diameter_threshold and circularity >= 0.8:
                x, y, w, h = cv2.boundingRect(contour)
                if (x > edge_threshold and y > edge_threshold
                    and (x + w) < (width - edge_threshold)
                    and (y + h) < (height - edge_threshold)):
                    valid_contours.append((contour, circularity))
    return valid_contours

def process_images_in_folder(folder_path):
    results = []
    object_data = pd.DataFrame()
    #filenames = [f for f in os.listdir(folder_path) if f.endswith('.jpg')]
    filenames = [f for f in os.listdir(folder_path) if f.endswith(('.jpg', '.tif'))]
    for filename in tqdm(filenames, desc="Processing images"):
        image_path = os.path.join(folder_path, filename)
        image_bgr = cv2.imread(image_path)
        image_rgb = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2RGB)
        masks = process_image(image_path)
        valid_contours = []
        valid_contours = filter_contours(masks, image_rgb.shape)
        if valid_contours:  # Only process if there are valid contours
            num_circles = len(valid_contours)
            area = [cv2.contourArea(c[0]) for c in valid_contours]
            diameters = 2*(np.array(area)/np.pi)**(0.5)  # Width as diameter
            diameters = diameters/1.29
            mean_diameter = np.mean(diameters)
            mean_circularity = np.mean([c[1] for c in valid_contours]) if valid_contours else 0
            droplets_cv = np.std(diameters) / mean_diameter
            # Draw contours on the original image
            overlay_image = image_rgb.copy()
            cv2.drawContours(overlay_image, [c[0] for c in valid_contours], -1, (0, 255, 0), 1)
            # Save overlay image
            overlay_path = os.path.join(folder_path, f'overlay_{os.path.splitext(filename)[0]}.png')
            cv2.imwrite(overlay_path, overlay_image)
            # Store results
            results.append({
                'Image': filename,
                'Num Circles': num_circles,
                'Mean Diameter': mean_diameter,
                'Mean Circularity': mean_circularity,
                'CV': droplets_cv
            })
            # store raw data
            df = {
                'Diameter': diameters,
                'Circularity': [c[1] for c in valid_contours],
                'Area': [cv2.contourArea(c[0]) / 1.29 / 1.29  for c in valid_contours],
                'Image': filename
            }
            df = pd.DataFrame(df)
            object_data = pd.concat([object_data, df], ignore_index=True)
        summary_df = pd.DataFrame(results)
        summary_df = summary_df.sort_values(by='Image')
        object_data_df = pd.DataFrame(object_data)  # Create DataFrame for object data
        object_data_df = object_data_df.sort_values(by='Image')  # Create DataFrame for object data
    return summary_df, object_data_df

# Example usage
folder_path = "/import/home/lyuah/od_image/image/20240919/"
summary_results, object_data_results = process_images_in_folder(folder_path)
object_data_results['Assay'] = object_data_results['Image'].str.rsplit('-', n=1).str[0]
object_data_results = object_data_results[object_data_results['Diameter'] > 50]

summary = object_data_results.groupby('Assay')['Diameter'].agg(
    Average_diameter='mean',
    CV_percentage=calculate_cv
).reset_index()
summary['Average_diameter'] = summary['Average_diameter'].round(2)
summary['CV_percentage'] = summary['CV_percentage'].round(2)
