"""Recover the in-plane mean, through-thickness fluctuation and shear work from saved GP fields."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


def in_plane(values, engineering_shear):
    a, b, shear = values
    c = shear / 2 if engineering_shear else shear
    return np.array([[a, c], [c, b]])


def contributions(element):
    work = np.zeros(3)
    for resultant in element['resultants']:
        mid_axes = np.array(resultant['mid_axes'])[:2]
        points = [p for p in element['fields'] if p['rs'] == resultant['rs']]
        weights = np.array([p['weight'] for p in points])
        strains, stresses = [], []
        for point in points:
            rotation = mid_axes @ np.array(point['material_axes'])[:2].T
            strains.append(rotation @ in_plane(point['strain'], True) @ rotation.T)
            stresses.append(rotation @ in_plane(point['stress_kpa'], False) @ rotation.T)
            work[2] += point['weight'] * np.dot(point['shear_kpa'], point['shear_strain']) / 2
        strains, stresses = np.array(strains), np.array(stresses)
        mean_strain = np.average(strains, axis=0, weights=weights)
        mean_stress = np.average(stresses, axis=0, weights=weights)
        work[0] += weights.sum() * np.sum(mean_strain * mean_stress) / 2
        work[1] += np.einsum('i,ijk,ijk', weights, strains - mean_strain, stresses - mean_stress) / 2
    return work


def collect(folders, output):
    with output.open('x', encoding='utf-8', newline='') as fp:
        writer = csv.writer(fp, delimiter='\t', lineterminator='\n')
        writer.writerow(['series', 'case', 'in_plane_mean', 'in_plane_variation', 'transverse_shear',
                         'artificial_drilling', 'physical_sum_minus_recorded', 'N_max', 'M_max', 'Q_max'])
        for folder in folders:
            for input_path in sorted(folder.glob('*.input.json')):
                path = input_path.with_name(input_path.name.removesuffix('.input.json') + '.json')
                if not path.exists():
                    continue  # Failed solve has no accepted endpoint; its diagnostic is in the series.
                elements = json.loads(path.read_bytes())['frames'][-1]['active']['elements']
                work = sum((contributions(e) for e in elements), np.zeros(3))
                physical = sum(e['physical_energy_kn_m'] for e in elements)
                artificial = sum(e['drilling_energy_kn_m'] for e in elements)
                sections = [r for e in elements for r in e['resultants']]
                maxima = [max(np.linalg.norm(r[key]) for r in sections) for key in ('n_kn_m','m_kn','q_kn_m')]
                writer.writerow([folder.name, path.stem, *work, artificial, work.sum() - physical, *maxima])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    parser.add_argument('series', nargs='+', type=Path)
    args = parser.parse_args()
    collect(args.series, args.output)
