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
import h5py


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
        center: list[float] = [np.random.uniform(-71.09, -71.07),
                               np.random.uniform(42.33, 42.34)]      # Random center within bounds
        # center = [-71.08793547508023, 42.33724270872659]
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
        elevation_map = bostwin.get_elevation_map(resolution=resolution)

        # Configure antenna array for all transmitters
        sionna_scene.tx_array = PlanarArray(num_rows=8,
                                num_cols=2,
                                vertical_spacing=0.7,
                                horizontal_spacing=0.5,
                                pattern="tr38901",
                                polarization="VH")

        # Configure antenna array for all receivers
        sionna_scene.rx_array = PlanarArray(num_rows=1,
                                num_cols=1,
                                vertical_spacing=0.5,
                                horizontal_spacing=0.5,
                                pattern="dipole",
                                polarization="cross")

        # Convert center to local CRS to align tx position
        tx_position_local = bostwin.lonlat2local([center])
        # tx_position_local = tx_position_local[0]
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
                cm_cell_size=(resolution*2, resolution*2),
                combining_vec=None,
                precoding_vec=None,
                num_samples=int(1e6),
            )
            path_gain: np.ndarray = coverage_map.path_gain.numpy().squeeze()
        except Exception as e:
            print(f"Cannot compute coverage map for {center}. Error: {e}")
            continue

        # Calculate the target size
        target_size = int(area_radius / resolution)

        

        
        def crop_or_pad(map_data, target_size):
            current_size = map_data.shape[0]
            if current_size < target_size:
                pad_size = (target_size - current_size) // 2
                return np.pad(map_data, pad_width=pad_size, mode='constant', constant_values=0)
            elif current_size > target_size:
                start = (current_size - target_size) // 2
                return map_data[start:start + target_size, start:start + target_size]
            else:
                return map_data
            
        # Crop or pad elevation map and path gain
        elevation_map = crop_or_pad(elevation_map, target_size)
        path_gain = crop_or_pad(path_gain, target_size)

        center_str = f"{center[0]:.4f}_{center[1]:.4f}".replace('.', 'p').replace('-', 'm')
        output_file = os.path.join(output_dir, f"data_{center_str}_{i}.h5")
        with h5py.File(output_file, 'w') as hf:
            hf.create_dataset('center', data=center)
            hf.create_dataset('resolution', data=resolution)
            hf.create_dataset('elevation_map', data=elevation_map)
            hf.create_dataset('path_gain', data=path_gain)
        # Plot Elevation and Coverage Maps
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))
        
        # Elevation Map
        axs[0].imshow(elevation_map, cmap='terrain', origin='lower')
        axs[0].set_title('Elevation Map')
        axs[0].set_xlabel('X (meters)')
        axs[0].set_ylabel('Y (meters)')
        
        # Coverage Map
        axs[1].imshow(10 * np.log10(path_gain), cmap='viridis', origin='lower')
        axs[1].set_title('Coverage Map')
        axs[1].set_xlabel('X (pixels)')
        axs[1].set_ylabel('Y (pixels)')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/map_{center_str}_{i}.png")
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