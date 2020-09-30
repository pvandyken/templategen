rule cp_nii_to_templateflow_naming:
    input: 
        template =  'results/iter_{iteration}/template_{{channel}}.nii.gz'.format(iteration=config['max_iters'])
    output:
        template =  'results/tpl-{name}/tpl-{name}_res-{res}_{{channel}}.nii.gz'.format(name=config['template_name'],res=config['resolution_index'])
    shell: 'cp {input} {output}'


rule create_template_json:
    input: 
        #needs a template nii to get origin, space, zooms
        nii = 'results/tpl-{name}/tpl-{name}_res-{res}_{channel}.nii.gz'.format(name=config['template_name'],res=config['resolution_index'],channel=channels[0])
    output:
        json = 'results/tpl-{name}/template_description.json'
    script: '../scripts/create_template_description_json.py'
