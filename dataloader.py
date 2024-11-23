"""
Pytorch dataloader 
"""
import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np
import os

class CoverageDataset(Dataset):
    """
    PyTorch Dataset for elevation maps and coverage maps.
    """

    def __init__(self, data_dir, transform=None):
        """
        Args:
            data_dir (str): Path to the dataset directory.
            transform (callable, optional): Optional transform to be applied
                on a sample (e.g., normalization, data augmentation).
        """
        self.data_dir = data_dir
        self.transform = transform

        # List all elevation and coverage map files
        self.elevation_files = sorted(
            [f for f in os.listdir(data_dir) if f.startswith("elevation_map") and f.endswith(".npy")]
        )
        self.coverage_files = sorted(
            [f for f in os.listdir(data_dir) if f.startswith("coverage_map") and f.endswith(".npy")]
        )
        assert len(self.elevation_files) == len(self.coverage_files), (
            "Mismatch between number of elevation maps and coverage maps"
        )

    def __len__(self):
        return len(self.elevation_files)

    def __getitem__(self, idx):
        """
        Fetch a sample from the dataset.
        Args:
            idx (int): Index of the sample to fetch.

        Returns:
            dict: A dictionary containing 'elevation' and 'coverage' tensors.
        """
        # Load elevation and coverage maps
        elevation_path = os.path.join(self.data_dir, self.elevation_files[idx])
        coverage_path = os.path.join(self.data_dir, self.coverage_files[idx])

        elevation_map = np.load(elevation_path).astype(np.float32)
        coverage_map = np.load(coverage_path).astype(np.float32)

        # Add channel dimension for U-Net compatibility
        elevation_map = np.expand_dims(elevation_map, axis=0)  # Shape: (1, H, W)
        coverage_map = np.expand_dims(coverage_map, axis=0)    # Shape: (1, H, W)

        sample = {
            "elevation": torch.from_numpy(elevation_map),
            "coverage": torch.from_numpy(coverage_map),
        }

        # Apply any transformations if provided
        if self.transform:
            sample = self.transform(sample)

        return sample

# Define a DataLoader
def get_dataloader(data_dir, batch_size=16, shuffle=True, num_workers=4, transform=None):
    """
    Creates a DataLoader for the dataset.
    Args:
        data_dir (str): Path to the dataset directory.
        batch_size (int): Number of samples per batch.
        shuffle (bool): Whether to shuffle the data.
        num_workers (int): Number of subprocesses to use for data loading.
        transform (callable, optional): Optional transform to be applied on each sample.

    Returns:
        DataLoader: PyTorch DataLoader instance.
    """
    dataset = CoverageDataset(data_dir, transform=transform)
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
    )
    return dataloader

# Example Usage
if __name__ == "__main__":
    # Path to the directory where the dataset is stored
    dataset_dir = "training_data"
    
    # Get the DataLoader
    dataloader = get_dataloader(dataset_dir, batch_size=8, shuffle=True, num_workers=2)
    
    # Iterate through the DataLoader
    for batch_idx, batch in enumerate(dataloader):
        elevation = batch["elevation"]  # Shape: (B, 1, H, W)
        coverage = batch["coverage"]    # Shape: (B, 1, H, W)
        print(f"Batch {batch_idx + 1}:")
        print(f"  Elevation shape: {elevation.shape}")
        print(f"  Coverage shape: {coverage.shape}")
        break