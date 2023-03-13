rule cp_nii_to_templateflow_naming:
    input: 
        template =  lambda wildcards: 'results/cohort-{cohort}/iter_{iteration}/template_{channel}.nii.gz'.format(iteration=config['num_iters'],channel=wildcards.channel,cohort=cohorts[int(wildcards.cohort_ind)-1])
    output:
        template =  'results/tpl-{name}/cohort-{{cohort_ind}}/tpl-{name}_cohort-{{cohort_ind}}_res-{res}_{{channel}}.nii.gz'.format(name=config['template_name'],res=config['resolution_index'])
    group: 'tf'
    shell: 'cp {input} {output}'


rule gen_brainmask_from_brain:
    input:
        t1w =  'results/tpl-{name}/cohort-{{cohort_ind}}/tpl-{name}_cohort-{{cohort_ind}}_res-{res}_T1w.nii.gz'.format(name=config['template_name'],res=config['resolution_index'])
    output:
        mask =  'results/tpl-{name}/cohort-{{cohort_ind}}/tpl-{name}_cohort-{{cohort_ind}}_res-{res}_desc-brain_mask.nii.gz'.format(name=config['template_name'],res=config['resolution_index'])
    group: 'tf'
    container: config['singularity']['prepdwi']
    shell: 'fslmaths {input} -bin {output}'

rule create_template_json:
    input: 
        #needs a template nii to get origin, space, zooms - use all of them as inputs, so the json gets created last
        expand('results/tpl-{name}/cohort-{cohort_ind}/tpl-{name}_cohort-{cohort_ind}_res-{res}_{channel}.nii.gz',name=config['template_name'],res=config['resolution_index'],channel=channels,cohort_ind=[i+1 for i in range(len(cohorts))] )
    params:
        cohorts = cohorts,
        subjects = subjects
    group: 'tf'
    output:
        json = 'results/tpl-{name}/template_description.json'
    script: '../scripts/create_template_description_json.py'
