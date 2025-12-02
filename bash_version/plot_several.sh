#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
import argparse

def main():
    parser = argparse.ArgumentParser(description='Plot RMSD or other trajectory-based data.')
    parser.add_argument('sim_name', help='Simulation folder name')
    parser.add_argument('property', help='Property name, e.g., RMSD or Gyration')

    args = parser.parse_args()

    data_path = os.path.join(os.getcwd(), args.sim_name, 'npt', 'cell', 'analysis', args.property, 'gyrate.txt')

    if not os.path.isfile(data_path):
        print(f"Error: File not found at {data_path}")
        return

    data = np.loadtxt(data_path)
    time = data[:, 0]
    values = data[:, 1]

    mean_y = np.mean(values)
    stddev_y = np.std(values)

    stats_path = os.path.join(os.getcwd(), args.sim_name, 'npt', 'cell', 'analysis', args.property, 'stats.txt')
    with open(stats_path, 'w') as f:
        f.write(f"Mean {args.property}: {mean_y:.4f}\n")
        f.write(f"Standard deviation {args.property}: {stddev_y:.4f}\n")

    plt.figure(figsize=(10, 6))
    plt.plot(time, values, label=args.property, linewidth=2, color='black')
    plt.axhline(mean_y, color='blue', linestyle='--', linewidth=2, label='Mean')
    plt.fill_between(time, mean_y - stddev_y, mean_y + stddev_y, color='lightgray', alpha=0.5, label='±1 SD')

    plt.xlabel('Time (ns)', fontsize=14)
    plt.ylabel(f'{args.property} (nm)', fontsize=14)
    plt.title(f'{args.sim_name} - {args.property}', fontsize=16)
    plt.legend()
    plt.grid(True)

    output_file = os.path.join(os.getcwd(), args.sim_name, 'npt', 'cell', 'analysis', args.property, f'{args.property}.png')
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved to {output_file}")

if __name__ == '__main__':
    main()
