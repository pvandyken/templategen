from os.path import join
import itertools
import pandas as pd


if config['run_cohort'] == None:
    cohorts = config['cohorts']
else:
    cohorts = [config['run_cohort']]

#subjects is a dict, indexed by the cohort
#load participants.tsv file, and strip off sub- from participant_id column
subjects = dict()
all_subjects = list()
for cohort in config['cohorts']:
    subjects[cohort] = [ s.strip('sub-') for s in pd.read_table(config['participants_tsv'][cohort]).participant_id.to_list() ] 
    all_subjects = all_subjects + subjects[cohort]

#get channels using keys in in_images
channels = list(config['in_images'].keys())