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


def generate_dataset(
                    bostwin: BostonTwin, 
                    num_samples: int = 100, 
                    output_dir: str = "out", 
                    resolution: float = 5.0,                        # in meters
                    area_radius: float = 100.0                     # it should be in meters, we should check 
                ) -> None:
    
    os.makedirs(output_dir, exist_ok=True)

    for i in range(num_samples):
        scene_name: str = f"scene_dummy"
        # center: list[float] = [np.random.uniform(-71.09, -71.07),
        #                        np.random.uniform(42.33, 42.34)]      # Random center within bounds
        center = [-71.08793547508023, 42.33724270872659]
        try: 
            # Generate the scene around the chosen center
            bostwin.generate_scene_from_radius(scene_name=scene_name,
                                                center_lon=center[0],
                                                center_lat=center[1],
                                                side_m=area_radius,
                                                load=True,
                                                )
            sionna_scene, _ = bostwin.load_scene(scene_name)
        except Exception as e:
            print(f"Cannot produce the 3D map for {center}. Error: {e}")
            continue

        # Elevation Map
        elevation_map, transform = bostwin.get_elevation_map(resolution=resolution)

        # Convert center to local CRS to align tx position
        tx_position_local = bostwin.lonlat2local([center])
        tx_position_local = tx_position_local[0]
        tx = Transmitter(
            name="tx",
            position=[tx_position_local[0], tx_position_local[1], 30],  # Use local CRS for center
            orientation=[0, 0, 0]
        )
        sionna_scene.add(tx)

        # Coverage Map
        try:
            coverage_map = sionna_scene.coverage_map(
                max_depth=10,
                diffraction=True,
                cm_cell_size=(resolution, resolution),
                combining_vec=None,
                precoding_vec=None,
                num_samples=int(1e6),
            )
            path_gain: np.ndarray = coverage_map.path_gain.numpy().squeeze()
        except Exception as e:
            print(f"Cannot compute coverage map for {center}. Error: {e}")
            continue

        # Save Elevation and Coverage Maps
        np.save(f"{output_dir}/elevation_map_{i}.npy", elevation_map)
        np.save(f"{output_dir}/coverage_map_{i}.npy", path_gain)

        # Plot Elevation and Coverage Maps
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))
        
        # Elevation Map
        axs[0].imshow(elevation_map, cmap='terrain', origin='lower', extent=(
            transform[2], transform[2] + transform[0] * elevation_map.shape[1],
            transform[5] + transform[4] * elevation_map.shape[0], transform[5]
        ))
        axs[0].set_title('Elevation Map')
        axs[0].set_xlabel('X (meters)')
        axs[0].set_ylabel('Y (meters)')
        
        # Coverage Map
        axs[1].imshow(10 * np.log10(path_gain), cmap='viridis', origin='lower')
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
        resolution=5, 
        area_radius=500
    )