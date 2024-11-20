"""
for generating dataset make sure that the Sionna Cnfig is at the highset as possible. 
TODO
1. size of elevation map and covaerage map are not the same 
2. boundaries are not set at their largest 
"""

import os
import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from sionna.rt import Transmitter, PlanarArray
from bostontwin.classes.BostonTwin import BostonTwin

# Generate dataset for training
def generate_dataset(
                    bostwin: BostonTwin, 
                    num_samples: int = 100, 
                    output_dir: str = "dataset", 
                    resolution: float = 5.0,                        # in meters
                    area_radius: float = 100.0,                     # it should be in meters, we should check 
                ) -> None:
    
    os.makedirs(output_dir, exist_ok=True)
    for i in range(num_samples):
        scene_name: str = f"scene"
        center: list[float] = [np.random.uniform(-71.09, -71.07),
                               np.random.uniform(42.33, 42.34)]      # can we find this boundaries from Boston Twin? 

        # elevation map
        try: 
            bostwin.generate_scene_from_radius(scene_name=scene_name,
                                                center_lon=center[0],
                                                center_lat=center[1],
                                                side_m=area_radius,
                                                load=True,
                                                )
            sionna_scene = bostwin.load_scene(scene_name)
        except:
            print(f'can not produce the 3d map for {center}')
            continue

        elevation_map = bostwin.get_elevation_map(resolution=resolution)
        np.save(f"{output_dir}/elevation_map_{i}.npy", elevation_map)
        
        antenna_array = PlanarArray(
            num_rows=1,
            num_cols=1,
            vertical_spacing=0.5,
            horizontal_spacing=0.5,
            pattern="iso",
            polarization="V",
        )
        sionna_scene.tx_array = antenna_array
        sionna_scene.rx_array = antenna_array
        tx0 = Transmitter(
            name="tx0",
            position=[0, 0, 20],  # Center of the scene
            orientation=[np.pi * 5 / 6, 0, 0],
            power_dbm=44,
        )
        sionna_scene.add(tx0)
        coverage_map = sionna_scene.coverage_map(
                                                max_depth=10,
                                                diffraction=True,
                                                cm_cell_size=(resolution, resolution),
                                                combining_vec=None,
                                                precoding_vec=None,
                                                num_samples=int(1e6),
                                                )
        path_gain: np.ndarray = coverage_map.path_gain.numpy().squeeze()
        np.save(f"{output_dir}/coverage_map_{i}.npy", path_gain) 

        # Plot elevation and coverage map
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))
        
        # Elevation map
        axs[0].imshow(elevation_map, cmap='terrain', origin='lower')
        axs[0].set_title('Elevation Map')
        axs[0].set_xlabel('X (pixels)')
        axs[0].set_ylabel('Y (pixels)')
        
        # Coverage map
        axs[1].imshow(path_gain, cmap='viridis', origin='lower')
        axs[1].set_title('Coverage Map')
        axs[1].set_xlabel('X (pixels)')
        axs[1].set_ylabel('Y (pixels)')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/map_{i}.png")
        plt.close(fig)

        print(f"Sample {i + 1}/{num_samples} generated.")

# Main function
if __name__ == "__main__":
    os.environ["CUDA_VISIBLE_DEVICES"] = '1'    # GPU ID 
    os.environ["TF_CPP_MIN_LOG_LEVEL"] = "1"
    
    dataset_path = Path("dataset")
    bostwin = BostonTwin(dataset_path)

    generate_dataset(
        bostwin, 
        num_samples=int(1e1), 
        output_dir="training_data", 
        resolution=1, 
        area_radius=100
    )