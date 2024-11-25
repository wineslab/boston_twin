import os
import h5py
import torch
from torch.utils.data import Dataset, DataLoader


class BostonTwinDataset(Dataset):
    """
    PyTorch Dataset for BostonTwin-generated data.
    """
    def __init__(self, data_dir: str, transform=None):
        """
        Initialize the dataset.

        Parameters:
        - data_dir (str): Path to the directory containing the .h5 dataset files.
        - transform (callable, optional): Optional transform to apply to the data.
        """
        self.data_dir = data_dir
        self.file_list = [os.path.join(data_dir, f) for f in os.listdir(data_dir) if f.endswith('.h5')]
        self.transform = transform

    def __len__(self):
        return len(self.file_list)

    def __getitem__(self, idx):
        file_path = self.file_list[idx]
        with h5py.File(file_path, 'r') as f:
            elevation_map = f['elevation_map'][:]
            path_gain = f['path_gain'][:]
            center = f['center'][:]
            resolution = f['resolution'][()]
        
        # Normalize inputs if needed (optional)
        elevation_map = torch.tensor(elevation_map, dtype=torch.float32).unsqueeze(0)  # Add channel dimension
        path_gain = torch.tensor(path_gain, dtype=torch.float32).unsqueeze(0)          # Add channel dimension
        
        # Apply transformations
        if self.transform:
            elevation_map, path_gain = self.transform((elevation_map, path_gain))

        return {
            'elevation_map': elevation_map,
            'path_gain': path_gain,
            'center': torch.tensor(center, dtype=torch.float32),
            'resolution': torch.tensor(resolution, dtype=torch.float32),
        }


def create_dataloader(data_dir, batch_size=32, shuffle=True, num_workers=4):
    """
    Create a PyTorch DataLoader for the BostonTwin dataset.

    Parameters:
    - data_dir (str): Path to the directory containing .h5 dataset files.
    - batch_size (int): Number of samples per batch.
    - shuffle (bool): Whether to shuffle the data at every epoch.
    - num_workers (int): Number of subprocesses for data loading.

    Returns:
    - DataLoader: A PyTorch DataLoader for the dataset.
    """
    dataset = BostonTwinDataset(data_dir)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, num_workers=num_workers)
    return dataloader


# Example usage
if __name__ == "__main__":
    data_dir = "training_data"  # Path to your generated dataset
    batch_size = 8

    dataloader = create_dataloader(data_dir, batch_size=batch_size)

    for batch_idx, batch in enumerate(dataloader):
        elevation_maps = batch['elevation_map']  # Shape: [batch_size, 1, H, W]
        path_gains = batch['path_gain']          # Shape: [batch_size, 1, H, W]
        centers = batch['center']                # Shape: [batch_size, 2]
        resolutions = batch['resolution']        # Shape: [batch_size]

        print(f"Batch {batch_idx + 1}:")
        print(f"  Elevation maps shape: {elevation_maps.shape}")
        print(f"  Path gains shape: {path_gains.shape}")
        print(f"  Centers: {centers}")
        print(f"  Resolutions: {resolutions}")
