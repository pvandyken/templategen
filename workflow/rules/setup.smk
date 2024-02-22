from os.path import join
from pathlib import Path
from snakebids import generate_inputs, bids, set_bids_spec
import itertools as it
import more_itertools as itx
import pandas as pd
import shutil as sh
import tempfile
import subprocess as sp
import templateflow.api as tflow

set_bids_spec("v0_0_0")

def get_labels(label):
    cli = config.get(label, None)
    return cli.split(",") if isinstance(cli, str) else cli

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
    derivatives=config["derivatives"],
    participant_label=participant_label,
    exclude_participant_label=exclude_participant_label,
    pybidsdb_dir=config.get("pybidsdb_dir"),
    pybidsdb_reset=config.get("pybidsdb_reset"),
)

def _get_default_tmpdir():
    _dir = workflow.default_resources.parsed.get("tmpdir")
    if callable(_dir):
        return _dir({}, '', 0, 0, '')
    return _dir

work = Path(
    config.get(
        'workdir',
        tempfile.mkdtemp(
            prefix="greedy-template.",
            dir=_get_default_tmpdir(),
        )
    )
).resolve()
output_dir = Path(config['output_dir'])
qc = os.path.join(config["output_dir"], "qc")

def _uid(entities = None):
    return Path(
        bids(
            **{
                **inputs.subj_wildcards,
                **{
                    entity: f"{{{entity}}}" for entity in entities or []
                }
            }
        )
    ).name

def tempout(name):
    path = work.joinpath(name)
    if config.get('workdir'):
        return path
    return temp(path)

def log(rulename: str, *entities):
    uid = _uid(entities)
    return Path("code", "logs", rulename, uid).with_suffix(".log")

def benchmark(rulename: str, *entities):
    uid = _uid(entities)
    return Path("code", "benchmarks", rulename,  uid).with_suffix(".tsv")

def resource(path):
    return os.path.join(workflow.basedir, "..", "resources", path)


channels = list(inputs.keys())

def get_entity_filters(filters):
    return {
        key: val for key, val in filters.items()
        if key != "scope" and "_" not in key
    }

output_warps = [
    bids(
        output_dir/"template",
        **{
            **get_entity_filters(config['pybids_inputs'][channel]["filters"]),
            **inputs[channel].wildcards,
            "space": "{template}",
            "suffix": (
                config['pybids_inputs'][channel]['filters']['suffix']
                + ''.join(Path(inputs[channel].path).suffixes)
            )
        }
    )
    for channel in channels
]

output_templates = [
    {
        **get_entity_filters(config['pybids_inputs'][channel]["filters"]),
        "space": "{template}",
        "suffix": (
            config['pybids_inputs'][channel]['filters']['suffix']
            + ''.join(Path(inputs[channel].path).suffixes)
        )
    } if channel not in config['pybids_outputs'] else {
        **config['pybids_outputs'],
        "space": "{template}",
    }
    for channel in channels

]
