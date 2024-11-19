from pathlib import Path

from bostontwin.utils.BostonModelDownloader import BostonModelDownloader

# Download the data
in_model_dir = Path("tmp_data/downloaded_data")
in_model_dir.mkdir(parents=True, exist_ok=True)
out_dataset_dir = Path("dataset", "scenes")
bos_downloader = BostonModelDownloader(
        in_model_dir,
        out_dataset_dir)
# bos_downloader.download_data(save_dir=in_model_dir, extract_objs=True)
bos_downloader.generate_dataset(create_xml=True)
