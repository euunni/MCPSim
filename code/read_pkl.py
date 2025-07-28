import pickle
import os
import glob
import numpy as np
import re

inputDir = "./../pbs_output/FromPython/pkl"
outputDir = "./../pbs_output/FromPython/txt"

os.makedirs(outputDir, exist_ok=True)

pattern = "400µm600V_3.5_0.5_saturation_0.8_each_5ps_*.pkl"
inputFiles = glob.glob(os.path.join(inputDir, pattern))

total = len(inputFiles)
count = 0

for inputFile in inputFiles:
    count += 1

    fileName = os.path.basename(inputFile)
    # print(f"Processing filename: {fileName}")
    
    match = re.search(r'400µm600V_3\.5_0\.5_saturation_0\.8_each_5ps_(\d+)\.imgt1_(\d+)\.pkl', fileName)
    
    if match:
        jobID = match.group(1)
        index = match.group(2)
        outputFile = os.path.join(outputDir, f"400µm600V_3.5_0.5_saturation_0.8_each_5ps_{jobID}_{index}.txt")
    else:
        try:
            jobID = fileName.split(".imgt1.pkl")[0].split("_")[-1]
            outputFile = os.path.join(outputDir, f"400µm600V_3.5_0.5_saturation_0.8_each_5ps_{jobID}.txt")
        except Exception:
            jobID = "unknown"
            print(f"Failed to extract JobID, using: {jobID}")
            outputFile = os.path.join(outputDir, f"400µm600V_3.5_0.5_saturation_0.8_each_5ps_{jobID}.txt")
    
    try:
        with open(inputFile, "rb") as f:
            data = pickle.load(f)
        
        with open(outputFile, "w") as f:
            f.write(f"Total number of secondary electrons: {len(data)}\n")
            
            electron_count = 0
            for electron in data:
                try:
                    if isinstance(electron, np.ndarray) and electron.shape == (3, 3):
                        x, y, z = electron[0][0], electron[0][1], electron[0][2]
                        vx, vy, vz = electron[1][0], electron[1][1], electron[1][2]
                        time = electron[2][0]
                    else:
                        print(f"Unknown electron data format: {type(electron)}")
                        continue
                    
                    f.write(f"{x:.1f} {y:.6f} {z:.6f} {vx:.6f} {vy:.6f} {vz:.6f} {time:.3f} 0 0\n")
                    electron_count += 1
                except Exception as e:
                    print(f"Error processing electron: {str(e)}")
                    
        print(f"[{count}/{total}] Saved {outputFile}")
        
    except Exception as e:
        print(f"Error in {inputFile} - {str(e)}")

print("Completed!")