

rule gen_init_avg_template:
    input: lambda wildcards: expand(config['in_images'][wildcards.channel],subject=subjects[wildcards.cohort])
    output: 'results/cohort-{cohort}/iter_0/init/init_avg_template_{channel}.nii.gz'
    params:
        dim = config['ants']['dim'],
        use_n4 = '2'
    log: 'logs/gen_init_avg_template_{channel}_{cohort}.log'
    container: config['singularity']['ants']
    group: 'init_template'
    shell:
        'AverageImages {params.dim} {output} {params.use_n4} {input} &> {log}'

rule get_existing_template:
    input: lambda wildcards: config['init_template'][wildcards.channel]
    output: 'results/cohort-{cohort}/iter_0/init/existing_template_{channel}.nii.gz'
    log: 'logs/get_existing_template_{channel}_{cohort}.log'
    group: 'init_template'
    shell: 'cp -v {input} {output} &> {log}'


rule set_init_template:
    input:
        'results/cohort-{cohort}/iter_0/init/init_avg_template_{channel}.nii.gz' if config['init_template'] == None else 'results/cohort-{cohort}/iter_0/init/existing_template_{channel}.nii.gz'
    params: 
        cmd = lambda wildcards, input, output:
                'ResampleImageBySpacing {dim} {input} {output} {vox_dims}'.format(
                        dim = config['ants']['dim'], input = input, output = output,
                        vox_dims=' '.join([str(d) for d in config['resample_vox_dims']]))
                     if config['resample_init_template'] else f"cp -v {input} {output}"
    output: 'results/cohort-{cohort}/iter_0/template_{channel}.nii.gz'
    log: 'logs/set_init_template_{channel}_{cohort}.log'
    group: 'init_template'
    container: config['singularity']['ants']
    shell: '{params.cmd} &> {log}'

rule reg_to_template:
    input: 
        template = lambda wildcards: ['results/cohort-{cohort}/iter_{iteration}/template_{channel}.nii.gz'.format(
                                iteration=iteration,channel=channel,cohort=wildcards.cohort) for iteration,channel in itertools.product([int(wildcards.iteration)-1],channels)],
        target = lambda wildcards: [config['in_images'][channel] for channel in channels]
    params:
        input_fixed_moving = lambda wildcards, input: [
            f'-i {fixed} {moving}'
            for fixed,moving in zip(input.template, input.target)
        ],
        input_moving_warped = lambda wildcards, input, output: [
            f'-rm {moving} {warped}'
            for moving,warped in zip(input.target,output.warped)
        ],
    output:
        warp = 'results/cohort-{cohort}/iter_{iteration}/sub-{subject}_1Warp.nii.gz',
        invwarp = 'results/cohort-{cohort}/iter_{iteration}/sub-{subject}_1InverseWarp.nii.gz',
        affine = 'results/cohort-{cohort}/iter_{iteration}/sub-{subject}_0GenericAffine.mat',
        affine_xfm_ras = 'results/cohort-{cohort}/iter_{iteration}/sub-{subject}_affine_ras.txt',
        warped = expand('results/cohort-{cohort}/iter_{iteration}/sub-{subject}_WarpedToTemplate_{channel}.nii.gz',channel=channels,allow_missing=True)
    log: 'logs/reg_to_template/cohort-{cohort}/iter_{iteration}_sub-{subject}.log'
    threads: 4
    group: 'reg'
    container: config['singularity']['itksnap']
    resources:
        # this is assuming 1mm
        mem_mb = 16000,
        time = 30
    shell: 
        #affine first
        'greedy -d 3 -threads {threads} -a -m NCC 2x2x2 '
        '{params.input_fixed_moving} -o {output.affine_xfm_ras} '
        '-ia-image-centers -n 100x50x10 &> {log} && '
        #then deformable:
        'greedy -d 3 -threads {threads} -m NCC 2x2x2 '
        '{params.input_fixed_moving} -it {output.affine_xfm_ras} '
        '-o {output.warp} -oinv {output.invwarp} -n 100x50x10 &>> {log} && '
        #then convert affine to itk format that ants uses
        'c3d_affine_tool {output.affine_xfm_ras} -oitk {output.affine} &>> {log} && '
        #and finally warp the moving image
        'greedy -d 3 -threads {threads} -rf {input.template[0]} '
        '{params.input_moving_warped} -r {output.warp} {output.affine_xfm_ras} '
        '&>> {log}'

rule avg_warped:
    input: 
        targets = lambda wildcards: expand('results/cohort-{cohort}/iter_{iteration}/sub-{subject}_WarpedToTemplate_{channel}.nii.gz',subject=subjects[wildcards.cohort],iteration=wildcards.iteration,channel=wildcards.channel,cohort=wildcards.cohort,allow_missing=True)
    params:
        dim = config['ants']['dim'],
        use_n4 = '0'  # changed to no normalization
    output: 'results/cohort-{cohort}/iter_{iteration}/shape_update/avg_warped_{channel}.nii.gz'
    group: 'shape_update'
    log: 'logs/avg_warped/cohort-{cohort}/iter_{iteration}_{channel}.log'
    container: config['singularity']['ants']
    shell:
        'AverageImages {params.dim} {output} {params.use_n4} {input} &> {log}'
       
rule avg_inverse_warps:
    input:
        warps = lambda wildcards: expand('results/cohort-{cohort}/iter_{iteration}/sub-{subject}_1Warp.nii.gz',subject=subjects[wildcards.cohort],iteration=wildcards.iteration,cohort=wildcards.cohort,allow_missing=True),
    params:
        dim = config['ants']['dim'],
        use_n4 = '0'
    output: 
        invwarp = 'results/cohort-{cohort}/iter_{iteration}/shape_update/avg_inverse_warps.nii.gz'
    group: 'shape_update'
    log: 'logs/avg_inverse_warps/cohort-{cohort}/iter_{iteration}.log'
    container: config['singularity']['ants']
    shell:
        'AverageImages {params.dim} {output} {params.use_n4} {input} &> {log}'
         
rule scale_by_gradient_step:
    input: 'results/cohort-{cohort}/iter_{iteration}/shape_update/avg_inverse_warps.nii.gz'
    params:
        dim = config['ants']['dim'],
        gradient_step = '-{gradient_step}'.format(gradient_step = config['ants']['shape_update']['gradient_step'])
    output: 'results/cohort-{cohort}/iter_{iteration}/shape_update/avg_inverse_warps_scaled.nii.gz'
    group: 'shape_update'
    log: 'logs/scale_by_gradient_step/cohort-{cohort}/iter_{iteration}.log'
    container: config['singularity']['ants']
    shell:
        'MultiplyImages {params.dim} {input} {params.gradient_step} {output} &> {log}' 

rule avg_affine_transforms:
    input: 
        affine = lambda wildcards: expand('results/cohort-{cohort}/iter_{iteration}/sub-{subject}_0GenericAffine.mat',subject=subjects[wildcards.cohort],iteration=wildcards.iteration,cohort=wildcards.cohort,allow_missing=True),
    params:
        dim = config['ants']['dim']
    output:
        affine = 'results/cohort-{cohort}/iter_{iteration}/shape_update/avg_affine.mat'
    group: 'shape_update'
    log: 'logs/avg_affine_transforms/cohort-{cohort}/iter_{iteration}.log'
    container: config['singularity']['ants']
    shell:
        'AverageAffineTransformNoRigid {params.dim} {output} {input} &> {log}'

rule transform_inverse_warp:
    input:
        affine = 'results/cohort-{cohort}/iter_{iteration}/shape_update/avg_affine.mat',
        invwarp = 'results/cohort-{cohort}/iter_{iteration}/shape_update/avg_inverse_warps_scaled.nii.gz',
        ref = lambda wildcards: 'results/cohort-{cohort}/iter_{iteration}/shape_update/avg_warped_{channel}.nii.gz'.format(iteration=wildcards.iteration,channel=channels[0],cohort=wildcards.cohort) #just use 1st channel as ref
    params:
        dim = '-d {dim}'.format(dim = config['ants']['dim'])
    output: 
        invwarp = 'results/cohort-{cohort}/iter_{iteration}/shape_update/avg_inverse_warps_scaled_transformed.nii.gz'
    group: 'shape_update'
    log: 'logs/transform_inverse_warp/cohort-{cohort}/iter_{iteration}.log'
    container: config['singularity']['ants']
    shell:
        'antsApplyTransforms {params.dim} -e vector -i {input.invwarp} -o {output} -t [{input.affine},1] -r {input.ref} --verbose 1 &> {log}'

rule apply_template_update:
    input:
        template =  'results/cohort-{cohort}/iter_{iteration}/shape_update/avg_warped_{channel}.nii.gz',
        affine = 'results/cohort-{cohort}/iter_{iteration}/shape_update/avg_affine.mat',
        invwarp = 'results/cohort-{cohort}/iter_{iteration}/shape_update/avg_inverse_warps_scaled_transformed.nii.gz'
    params:
        dim = '-d {dim}'.format(dim = config['ants']['dim'])
    output:
        template =  'results/cohort-{cohort}/iter_{iteration}/template_{channel}.nii.gz'
    log: 'logs/apply_template_update/cohort-{cohort}/iter_{iteration}_{channel}_{cohort}.log'
    group: 'shape_update'
    container: config['singularity']['ants']
    shell:
        'antsApplyTransforms {params.dim} --float 1 --verbose 1 -i {input.template} -o {output.template} -t [{input.affine},1] '
        ' -t {input.invwarp} -t {input.invwarp} -t {input.invwarp} -t {input.invwarp} -r {input.template} &> {log}' #apply warp 4 times

rule get_final_xfm:
    input:
        template=expand(
            rules.apply_template_update.output['template'],
            iteration=config['max_iters']
        ),
        affine=expand(
            rules.reg_to_template.output['affine_xfm_ras'],
            iteration=config['max_iters']
        ),
        warp=expand(
            rules.reg_to_template.output['warp']
            iteration=config['max_iters']
        ),
    output:
        bids(
            output_dir/"template",
            datatype="xfm",
            mode="image",
            from_="individual",
            to="{template}",
            suffix="xfm.nii.gz"
            **inputs.subj_wildcards,
        )

    container: config['singularity']['itksnap']
    shell: 'greedy -d 3 -rf {input.template} '
          ' -r {input.warp} {input.affine} '
          ' -rc {output}'


rule get_final_inv_xfm:
    input:
        template=expand(
            rules.apply_template_update.output['template'],
            iteration=config['max_iters']
        ),
        affine=expand(
            rules.reg_to_template.output['affine_xfm_ras'],
            iteration=config['max_iters']
        ),
        warp=expand(
            rules.reg_to_template.output['invwarp']
            iteration=config['max_iters']
        ),
    output:
        bids(
            output_dir/"template",
            datatype="xfm",
            mode="image",
            from_="{template}",
            to="individual",
            suffix="xfm.nii.gz"
            **inputs.subj_wildcards,
        )

    container: config['singularity']['itksnap']
    shell: 'greedy -d 3 -rf {input.template} '
          ' -r {input.affine},-1 {input.warp}'
          ' -rc {output}'