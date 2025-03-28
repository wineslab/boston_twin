from src.classes.BostonTwin import BostonTwin
from argparse import ArgumentParser
from pathlib import Path


def main(args):
    scene_name = args.scene_name
    # Load the model
    bostwin = BostonTwin("dataset")

    # Generate the new scene
    bostwin.generate_scene_from_radius(
        scene_name=scene_name,
        center_lon=args.center_lon,
        center_lat=args.center_lat,
        radius=args.radius,
        load=True,
    )
    if args.out_dir:
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        bostwin.export_scene_models(scene_name=scene_name, out_dir=out_dir)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--scene_name", type=str, required=True, help="Name of the new scene")
    parser.add_argument("--center_lon", type=float, required=True, help="Longitude of the center of the new scene")
    parser.add_argument("--center_lat", type=float, required=True, help="Latitude of the center of the new scene")
    parser.add_argument("--radius", type=float, required=True, help="Radius of the new scene")
    parser.add_argument("--out_dir", type=str, default=None, required=False, help="Output directory to save the new scene")

    args = parser.parse_args()
    main(args)
