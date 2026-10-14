"""Read-only aggregate validation of all eight reduced-tree jobs."""
import argparse
import json
from pathlib import Path
import uproot

ERAS = ('2016preAPV', '2016postAPV', '2017', '2018')
CATEGORIES = ('resolved', 'merged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir', type=Path)
    args = parser.parse_args()
    reports = []
    for era in ERAS:
        for category in CATEGORIES:
            stem = f'{era}_{category}'
            if not (args.output_dir / f'{stem}.done').is_file():
                raise RuntimeError(f'Not finalized: {stem}')
            report = json.loads((args.output_dir / f'validation_{stem}.json').read_text())
            assert report['status'] == 'passed' and not report['pilot']
            assert report['era'] == era and report['category'] == category
            assert report['files'] == len(report['outputs']) == 72
            assert len(report['sources']) == 11
            assert len(report['schema']) == 16
            assert sum(row['bridge_entries'] for row in report['outputs'] if row['region'] == 'CR') > 0
            for row in report['outputs']:
                with uproot.open(args.output_dir / row['file']) as handle:
                    tree = handle['tree']
                    assert tree.num_entries == row['entries']
                    assert set(tree.keys()) == set(report['schema'])
            reports.append(report)
    assert len({r['runner_sha256'] for r in reports}) == 1
    assert len({r['producer_sha256'] for r in reports}) == 1
    for era in ERAS:
        assert len({r['config_sha256'] for r in reports if r['era'] == era}) == 1
    assert len(list(args.output_dir.glob('*.root'))) == 576
    assert not list(args.output_dir.glob('*.partial'))
    print(json.dumps(dict(status='passed', jobs=8, files=576,
                          eras=list(ERAS), per_job=[dict(era=r['era'], category=r['category'],
                          files=r['files'], elapsed_seconds=r['elapsed_seconds']) for r in reports]), indent=2))


if __name__ == '__main__':
    main()
