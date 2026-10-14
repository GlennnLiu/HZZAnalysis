"""Check every production input/header and region definition without writing trees."""
from argparse import Namespace
from pathlib import Path
import yaml
import produce_reduced_stream as runner


if __name__ == '__main__':
    for era in runner.producer.ERAS:
        config = yaml.safe_load((Path(__file__).resolve().parent/'cards' / f'config_UL{era}.yml').read_text())
        sf = yaml.safe_load((Path(__file__).resolve().parent/'sf' / f'NetSF_signal_{era}.yml').read_text())
        assert 'ZvsQCD' in sf and 'ZvsLight' in sf
        for category in runner.producer.CATEGORIES:
            inputs = runner.preflight(config, Namespace(era=era, cat=category))
            print(f'{era} {category}: passed, {len(inputs)} combined inputs', flush=True)
