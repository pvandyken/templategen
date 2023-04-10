rule reg_to_template:
    input: 
        template = lambda wildcards: expand(
            work/'iter_{iteration}/template_{channel}.nii.gz',
            iteration=int(wildcards.iteration) - 1,
            channel=channels,
        ),
        target = lambda wildcards: [
            inputs[channel].path.format(**wildcards) for channel in channels
        ]
    output:
        warp = work/'iter_{iteration}/sub-{subject}_1Warp.nii.gz',
        invwarp = work/'iter_{iteration}/sub-{subject}_1InverseWarp.nii.gz',
        affine = work/'iter_{iteration}/sub-{subject}_0GenericAffine.mat',
        affine_xfm_ras = work/'iter_{iteration}/sub-{subject}_affine_ras.txt',
        warped = expand(
            work/'iter_{iteration}/sub-{subject}_WarpedToTemplate_{channel}.nii.gz',
            channel=channels,
            allow_missing=True
        )
    log: 'logs/reg_to_template/iter_{iteration}_sub-{subject}.log'
    benchmark:  'benchmark/reg_to_template/iter_{iteration}_sub-{subject}.tsv'
    threads: 4
    group: 'reg'
    container: config['singularity']['itksnap']
    resources:
        # this is assuming 1mm
        mem_mb = 16000,
        runtime = 30
    params:
        input_fixed_moving = lambda wildcards, input: [
            f'-i {fixed} {moving}'
            for fixed,moving in zip(input.template, input.target)
        ],
        input_moving_warped = lambda wildcards, input, output: [
            f'-rm {moving} {warped}'
            for moving,warped in zip(input.target,output.warped)
        ],
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
        targets = expand(
           rules.reg_to_template.output['warped'],
           subject=inputs.subjects,
           allow_missing=True,
        )
    params:
        dim = config['ants']['dim'],
        use_n4 = '0'  # changed to no normalization
    output: work/'iter_{iteration}/shape_update/avg_warped_{channel}.nii.gz'
    group: 'shape_update'
    log: 'logs/avg_warped/iter_{iteration}_{channel}.log'
    container: config['singularity']['ants']
    shell:
        'AverageImages {params.dim} {output} {params.use_n4} {input} &> {log}'
       
rule avg_inverse_warps:
    input:
        warps = expand(
            rules.reg_to_template.output['warp'],
            subject=inputs.subjects,
            allow_missing=True
        ),
    params:
        dim = config['ants']['dim'],
        use_n4 = '0'
    output: 
        invwarp = work/'iter_{iteration}/shape_update/avg_inverse_warps.nii.gz'
    group: 'shape_update'
    log: 'logs/avg_inverse_warps/iter_{iteration}.log'
    container: config['singularity']['ants']
    shell:
        'AverageImages {params.dim} {output} {params.use_n4} {input} &> {log}'
         
rule scale_by_gradient_step:
    input: rules.avg_inverse_warps.output
    params:
        dim = config['ants']['dim'],
        gradient_step = f"-{config['ants']['shape_update']['gradient_step']}"
    output: work/'iter_{iteration}/shape_update/avg_inverse_warps_scaled.nii.gz'
    group: 'shape_update'
    log: 'logs/scale_by_gradient_step/iter_{iteration}.log'
    container: config['singularity']['ants']
    shell:
        'MultiplyImages {params.dim} {input} {params.gradient_step} {output} &> {log}' 

rule avg_affine_transforms:
    input: 
        affine = expand(
            rules.reg_to_template.output['affine'],
            subject=inputs.subjects,
            allow_missing=True
        ),
    params:
        dim = config['ants']['dim']
    output:
        affine = work/'iter_{iteration}/shape_update/avg_affine.mat'
    group: 'shape_update'
    log: 'logs/avg_affine_transforms/iter_{iteration}.log'
    container: config['singularity']['ants']
    shell:
        'AverageAffineTransformNoRigid {params.dim} {output} {input} &> {log}'

rule transform_inverse_warp:
    input:
        affine = rules.avg_affine_transforms.output,
        invwarp = rules.scale_by_gradient_step.output,
        ref = expand(
            rules.avg_warped.output,
            channel=channels[0],  # just use 1st channel as ref
            allow_missing=True
        )
    params:
        dim = '-d {dim}'.format(dim = config['ants']['dim'])
    output: 
        invwarp = work/'iter_{iteration}/shape_update/avg_inverse_warps_scaled_transformed.nii.gz'
    group: 'shape_update'
    log: 'logs/transform_inverse_warp/iter_{iteration}.log'
    container: config['singularity']['ants']
    shell:
        'antsApplyTransforms {params.dim} -e vector -i {input.invwarp} -o {output} -t [{input.affine},1] -r {input.ref} --verbose 1 &> {log}'

rule apply_template_update:
    input:
        template=rules.avg_warped.output,
        affine=rules.avg_affine_transforms.output,
        invwarp=rules.transform_inverse_warp.output,
    params:
        dim = '-d {dim}'.format(dim = config['ants']['dim'])
    output:
        template = work/'iter_{iteration}/template_{channel}.nii.gz'
    log: 'logs/apply_template_update/iter_{iteration}_{channel}.log'
    group: 'shape_update'
    container: config['singularity']['ants']
    shell:
        'antsApplyTransforms {params.dim} --float 1 --verbose 1 -i {input.template} -o {output.template} -t [{input.affine},1] '
        ' -t {input.invwarp} -t {input.invwarp} -t {input.invwarp} -t {input.invwarp} -r {input.template} &> {log}' #apply warp 4 times

rule get_final_xfm:
    input:
        template=expand(
            rules.apply_template_update.output['template'],
            iteration=config['num_iters'],
            channel=channels[0],
        ),
        affine=expand(
            rules.reg_to_template.output['affine_xfm_ras'],
            iteration=config['num_iters'],
            allow_missing=True,
        ),
        warp=expand(
            rules.reg_to_template.output['warp'],
            iteration=config['num_iters'],
            allow_missing=True,
        ),
    output:
        bids(
            output_dir/"template",
            datatype="xfm",
            mode="image",
            from_="individual",
            to="{template}",
            suffix="xfm.nii.gz",
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
            iteration=config['num_iters'],
            channel=channels[0],
        ),
        affine=expand(
            rules.reg_to_template.output['affine_xfm_ras'],
            iteration=config['num_iters'],
            allow_missing=True,
        ),
        warp=expand(
            rules.reg_to_template.output['invwarp'],
            iteration=config['num_iters'],
            allow_missing=True,
        ),
    output:
        bids(
            output_dir/"template",
            datatype="xfm",
            mode="image",
            from_="{template}",
            to="individual",
            suffix="xfm.nii.gz",
            **inputs.subj_wildcards,
        )

    container: config['singularity']['itksnap']
    shell: 'greedy -d 3 -rf {input.template} '
          ' -r {input.affine},-1 {input.warp}'
          ' -rc {output}'


rule get_final_warped:
    input:
        expand(
            rules.reg_to_template.output['warped'],
            iteration=config['num_iters'],
            channel=channels,
            allow_missing=True,
        )
    output: output_warps
    run:
        for src, dest in zip(input, output):
            sh.copyfile(src, dest)