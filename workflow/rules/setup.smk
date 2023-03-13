from os.path import join
from pathlib import Path
from snakebids import generate_inputs, bids
import itertools
import pandas as pd

def get_labels(label):
    cli = config.get(label, None)
    vals = cli if isinstance(cli, int) else cli.split(",") if cli is not None else None
    return vals

def get_participants(participant_file):
    subs = pd.read_csv(participant_file, sep='\t')
    return subs['participant_id'].map(lambda s: s[4:])

participant_label = get_labels("participant_label")
exclude_participant_label = (
    get_labels("exclude_participant_label") if participant_label is None else None
)
if participant_label is None and exclude_participant_label is None:
    try:
        participant_label = list(get_participants(
            Path(config['bids_dir'], 'participants.tsv')
        ))
    except FileNotFoundError:
        pass

###
# Input Globals
###
inputs = generate_inputs(
    bids_dir=config['bids_dir'],
    pybids_inputs=config['pybids_inputs'],
    derivatives=True,
    participant_label=participant_label,
    exclude_participant_label=exclude_participant_label,
    use_bids_inputs=True,
    pybids_database_dir=config.get("pybids_database_dir"),
    pybids_reset_database=config.get("pybids_reset_database"),
)

work = Path(
    config.get(
        'workdir',
        tempfile.mkdtemp(
            prefix="greedy-template.",
            dir=workflow.default_resources.get("tmpdir")
        )
    )
)
output_dir = Path(config['output_dir'])


#get channels using keys in in_images
channels = list(inputs.keys())