
# these rules are for registering the newly-generated cohort templates back 
# to the initial (MNI152) template, and composing the warps from template-building


rule reg_cohort_template_to_std:
    input:
        std_template = [ config['init_template'][channel]  for channel in channels ],
        cohort_template =  expand('results/cohort-{{cohort}}/iter_{iteration}/template_{channel}.nii.gz',iteration=config['max_iters'],channel=channels)
    params:
        input_fixed_moving = lambda wildcards, input: [f'-i {fixed} {moving}' for fixed,moving in zip(input.std_template, input.cohort_template) ],
        input_moving_warped = lambda wildcards, input, output: [f'-rm {moving} {warped}' for moving,warped in zip(input.cohort_template,output.warped) ],
    output:
        warp = 'results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_1Warp.nii.gz',
        invwarp = 'results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_1InverseWarp.nii.gz',
        affine_xfm_ras = 'results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_affine_ras.txt',
        warped = expand('results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_WarpedToTemplate_{channel}.nii.gz',channel=channels,allow_missing=True)
        
    log: 'logs/reg_cohort_template_to_std/cohort-{cohort}_{std_template}.log'
    threads: 4
    container: config['singularity']['itksnap']
    resources:
        # this is assuming 1mm
        mem_mb = 16000,
        time = 30
    shell: 
        #affine first
        'greedy -d 3 -threads {threads} -a -m NCC 2x2x2 {params.input_fixed_moving} -o {output.affine_xfm_ras} -ia-image-centers -n 100x50x10 &> {log} && '
        #then deformable:
        'greedy -d 3 -threads {threads} -m NCC 2x2x2 {params.input_fixed_moving} -it {output.affine_xfm_ras} -o {output.warp} -oinv {output.invwarp} -n 100x50x10 &>> {log} && '

        #and finally warp the moving image
        'greedy -d 3 -threads {threads} -rf {input.std_template[0]} {params.input_moving_warped} -r {output.warp} {output.affine_xfm_ras} &>> {log}'



