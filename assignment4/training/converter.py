#!/usr/bin/env python3
import argparse
import os,sys
from dataclasses import dataclass
from pathlib import Path
from random import randint
from typing import List

@dataclass(frozen=True)
class CompressionLevel:
    in_color: bool
    level: int

    def __post_init__(self):
        if self.level not in [0, 1, 2]:
            raise RuntimeError(f"CompressionLevel level was {self.level}, must be in [0, 1, 2]")

    @staticmethod
    def make_random(in_color: bool):
        return CompressionLevel(
            in_color=in_color,
            level=randint(0,2)
        )

    def cmd_str(self):
        cmd = "pnmtopng"

        if self.level == 0 and not self.in_color:
            cmd = "convert -intensity average -colorspace gray ppm:- ppm:- | " + cmd

        if not self.in_color:
            cjpeg_flag = "-grayscale"
        else:
            cjpeg_flag = ""
        for q in reversed(Q_VALS[self]):
            cmd = f"cjpeg {cjpeg_flag} -quality {q} | djpeg | " + cmd
        
        return cmd

Q_VALS = {
    CompressionLevel(True, 0): [],
    CompressionLevel(False, 0): [],
    CompressionLevel(True, 1): [90],
    CompressionLevel(True, 2): [90, 80],
    CompressionLevel(False, 1): [75],
    CompressionLevel(False, 2): [80, 70],   
}

@dataclass(frozen=True)
class CompressionJob:
    in_path: Path
    out_path: Path
    compression: CompressionLevel

    def cmd_str(self):
        return f"dcraw -c {self.in_path} | convert -scale 50% ppm:- ppm:- | {self.compression.cmd_str()} > {self.out_path}"

    def csv_str(self):
        return f"{self.out_path.name},{'color' if self.compression.in_color else 'gray'},{self.compression.level}"

    def __str__(self):
        return self.__repr__() + "\n" + self.cmd_str()

if __name__ == '__main__':
    p = argparse.ArgumentParser("DSP Image Generator", description="Given a set of raw images, randomly converts them to PNG with 0,1, or 2 levels of JPEG compression. Outputs a CSV of image-name, in-color, compression-level to stdout")

    p.add_argument("--dry_run", "-d", action='store_true', help="Print compression jobs and exit")
    p.add_argument("--output_dir", type=str, help="Folder to output images to", default=".")
    p.add_argument("--output_csv", type=str, help="CSV file to output compression stats to", default=None)
    p.add_argument("raw_images", nargs='+', help="Raw images in DNG format")

    args = p.parse_args()

    for r in args.raw_images:
        if not os.path.isfile(r):
            raise argparse.ArgumentError(f"{r} isn't a file/doesn't exist")
    if not os.path.isdir(args.output_dir):
        raise argparse.ArgumentError(f"{args.output_dir} isn't a directory/doesn't exist")

    raw_image_paths = [Path(r) for r in args.raw_images]

    compression_jobs: List[CompressionJob] = []
    for i, r_p in enumerate(raw_image_paths):
        compression_jobs.append(CompressionJob(
            in_path=r_p,
            out_path=Path(args.output_dir) / f"test{2*i:02d}.png",
            compression=CompressionLevel.make_random(False)
        ))
        compression_jobs.append(CompressionJob(
            in_path=r_p,
            out_path=Path(args.output_dir) / f"test{2*i+1:02d}.png",
            compression=CompressionLevel.make_random(True)
        ))

    # Output CSV
    if args.output_csv:
        with open(args.output_csv, "w") as f:
            f.write("\n".join([c.csv_str() for c in compression_jobs]))

    if args.dry_run:
        print("\n".join([c.cmd_str() for c in compression_jobs]))
        sys.exit(0)

    for i,job in enumerate(compression_jobs):
        print(f"Converting image {i} of {len(compression_jobs)}", sep="\r")
        os.system(job.cmd_str())
    print("")
