
rule gen_init_avg_template:
    input: lambda wildcards: expand(
        inputs[wildcards.channel].path,
        **inputs.subjects
    )
    output: work/'iter_0/init/init_avg_template_{channel}.nii.gz'
    output: work/'iter_0/template_{channel}.nii.gz'
    params:
        dim = config['ants']['dim'],
        use_n4 = '2',
        vox_dims=' '.join([str(d) for d in config['resample_vox_dims']])
    log: 'logs/gen_init_avg_template_{channel}.log'
    container: config['singularity']['ants']
    group: 'init_template'
    shadow: 'minimal'
    shell:
        "AverageImages {params.dim} {output} {params.use_n4} {input} &> {log}"


rule get_existing_template:
    input: lambda wildcards: config['init_template'][wildcards.channel]
    output: 'results/cohort-{cohort}/iter_0/init/existing_template_{channel}.nii.gz'
    output: work/'iter_0/init/existing_template_{channel}.nii.gz'
    log: 'logs/get_existing_template_{channel}_{cohort}.log'
    group: 'init_template'
    shell: 'cp -v {input} {output} &> {log}'


def choose_init_template(wcards):
    if config.get("init_template"):
        return rules.get_existing_template.output[0].format(**wcards)
    return rules.get_init_avg_template.output[0].format(**wcards)

rule resample_init_template:
    input: choose_init_template
    output: work/'iter_0/init/init_resampled_template_{channel}.nii.gz'
    log: 'logs/resample_init_avg_template_{channel}.log'
    container: config['singularity']['ants']
    params:
        dim = config['ants']['dim'],
        vox_dims=' '.join([str(d) for d in config['resample_vox_dims']])
    shell:
        "ResampleImageBySpacing {params.dim} {input} {output} {params.vox_dims}"

def choose_init_template_sampling(wcards):
    if config.get('resample_init_template'):
        return rules.resample_init_template.output[0].format(**wcards)
    return choose_init_template(wcards)


rule set_init_template:
    input: choose_init_template_sampling
    output: work/'iter_0/template_{channel}.nii.gz'
    group: 'init_template'
    shell: "mv {input} {output}"
