rule cp_nii_to_templateflow_naming:
    input: 
        expand(
            rules.apply_template_update.output,
            iteration=config['num_iters'],
            channel=channels,
        )
    output:
        [
            bids(
                output_dir/"templateflow/tpl-{name}",
                prefix="tpl-{name}",
                res="{res}",
                **wcards
            ).format(
                name=config['template_name'],
                res=config['resolution_index'],
            )
            for wcards in output_templates
        ]
    group: 'tf'
    run:
        for src, dest in zip(input, output):
            sh.copyfile(src, dest)



rule create_template_json:
    input: 
        #needs a template nii to get origin, space, zooms - use all of them as inputs, so the json gets created last
        expand('results/tpl-{name}/cohort-{cohort_ind}/tpl-{name}_cohort-{cohort_ind}_res-{res}_{channel}.nii.gz',name=config['template_name'],res=config['resolution_index'],channel=channels,cohort_ind=[i+1 for i in range(len(cohorts))] )
    params:
        cohorts = cohorts,
        subjects = subjects
    group: 'tf'
    output:
        json = output_dir/'templateflow/tpl-{name}/template_description.json'
    script: '../scripts/create_template_description_json.py'
