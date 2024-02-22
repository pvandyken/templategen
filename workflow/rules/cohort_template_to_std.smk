
# these rules are for registering the newly-generated cohort templates back 
# to the initial (MNI152) template, and composing the warps from template-building

rule get_avg_channel:
    input:
        inputs["FA"].expand(
            rules.get_final_warped.output,
            template=config["template_name"]
        )
    output:
        bids(
            output_dir/"tpl-{template}",
            prefix="tpl-{template}",
            suffix="{channel,FA}.nii.gz",
        )
    container: config['singularity']['itksnap']
    shell:
        "c3d {input} -mean -o {output}"
        

def _get_reference_template(wcards):
    return tflow.get(
        template=wcards["ref_template"],
        resolution="02",
        desc=None,
        suffix="T1w",
        extension=".nii.gz",
    )

rule reg_cohort_template_to_std:
    input:
        std_template = _get_reference_template,
        cohort_template = expand(
            rules.get_avg_channel.output,
            channel=channels,
            allow_missing=True,
        )
    output:
        warp = tempout('reg_to_{ref_template}/to-{ref_template}_from-{template}_1Warp.nii.gz'),
        invwarp = tempout('reg_to_{ref_template}/to-{template}_from-{ref_template}_1InverseWarp.nii.gz'),
        affine_xfm_ras = tempout('reg_to_{ref_template}/to-{ref_template}_from-{template}_affine_ras.txt'),
    log: 'logs/reg_cohort_template_to_std/{ref_template}.{template}.log'
    benchmark:  'benchmark/reg_cohort_template_to_std/{ref_template}.{template}.tsv'
    threads: 8
    container: config['singularity']['itksnap']
    resources:
        # this is assuming 1mm
        mem_mb = 32000,
        runtime = 60
    group: 'reg_to_std'
    params:
        input_fixed_moving = lambda wildcards, input: [
            f'-i {fixed} {moving}'
            for fixed,moving in it.zip_longest(
                itx.always_iterable(input.std_template),
                input.cohort_template
            )
        ],
    shell: 
        #affine first
        """
        greedy -d 3 -threads {threads} -a -m NCC 2x2x2 -n 100x50x10 -ia-image-centers \\
            {params.input_fixed_moving} -o {output.affine_xfm_ras} \\
            &> {log}
        """
        #then deformable:
        """
        greedy -d 3 -threads {threads} -m NCC 2x2x2 -n 100x50x10 \\
            {params.input_fixed_moving} -it {output.affine_xfm_ras} \\
            -o {output.warp} -oinv {output.invwarp} &>> {log}
        """

rule get_template_xfm:
    input:
        template=_get_reference_template,
        affine=rules.reg_cohort_template_to_std.output['affine_xfm_ras'],
        warp=rules.reg_cohort_template_to_std.output['warp'],
    output:
        bids(
            output_dir/"tpl-{template}",
            prefix="tpl-{template}",
            mode="image",
            from_="{template}",
            to="{ref_template}",
            suffix="xfm.nii.gz",
        )

    container: config['singularity']['itksnap']
    shell: 
        """
        greedy -d 3 -rf {input.template} \\
            -r {input.warp} {input.affine} -rc {output}
          """


rule get_template_inv_xfm:
    input:
        template=_get_reference_template,
        affine=rules.reg_cohort_template_to_std.output['affine_xfm_ras'],
        warp=(rules.reg_cohort_template_to_std.output['invwarp']),
    output:
        bids(
            output_dir/"tpl-{template}",
            prefix="tpl-{template}",
            mode="image",
            from_="{ref_template}",
            to="{template}",
            suffix="xfm.nii.gz",
        )

    container: config['singularity']['itksnap']
    shell: 
        """
        greedy -d 3 -rf {input.template} \\
            -r {input.affine},-1 {input.warp} -rc {output}
        """


def get_inputs_composite_subj_to_std (wildcards):
    """ Function for setting all the inputs
        Needed since cohort isn't in the output filename, and is determined
        by looking at the input lists
    """
    for c in cohorts:
        if wildcards.subject in subjects[c]:
            cohort = c
    std_template = wildcards.std_template
    subject = wildcards.subject
    iteration=config['num_iters']

    return {
        'cohort2std_warp': f'results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_1Warp.nii.gz',
        'cohort2std_affine_xfm_ras': f'results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_affine_ras.txt',
        'subj2cohort_warp': f'results/cohort-{cohort}/iter_{iteration}/sub-{subject}_1Warp.nii.gz',
        'subj2cohort_affine_xfm_ras': f'results/cohort-{cohort}/iter_{iteration}/sub-{subject}_affine_ras.txt',
        'ref_std': config['init_template'][channels[0]] }
      

rule create_composite_subj_to_std:
    """ This concatenates the subject to cohort to mni warps/affines to get a single warp from subject to mni """
    input: unpack(get_inputs_composite_subj_to_std)
    output:
        subj2std_warp = 'results/composite/sub-{subject}_to-{std_template}_via-cohort_CompositeWarp.nii.gz'
    group: 'composite'
    shell: 'greedy -d 3 -rf {input.ref_std} '
          ' -r {input.cohort2std_warp} {input.cohort2std_affine_xfm_ras} '
          '  {input.subj2cohort_warp} {input.subj2cohort_affine_xfm_ras} '
          ' -rc {output.subj2std_warp}'



def get_inputs_composite_subj_to_std_inverse (wildcards):
    """ Function for setting all the inputs
        Needed since cohort isn't in the output filename, and is determined
        by looking at the input lists
    """
    for c in cohorts:
        if wildcards.subject in subjects[c]:
            cohort = c
    std_template = wildcards.std_template
    subject = wildcards.subject
    iteration=config['num_iters']

    return {
        'cohort2std_invwarp': f'results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_1InverseWarp.nii.gz',
        'cohort2std_affine_xfm_ras': f'results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_affine_ras.txt',
        'subj2cohort_invwarp': f'results/cohort-{cohort}/iter_{iteration}/sub-{subject}_1InverseWarp.nii.gz',
        'subj2cohort_affine_xfm_ras': f'results/cohort-{cohort}/iter_{iteration}/sub-{subject}_affine_ras.txt',
        'ref_subj':   config['in_images'][channels[0]] }
      



rule create_composite_subj_to_std_inverse:
    """ This concatenates the subject to cohort to mni warps/affines to get a single warp from subject to std_template"""
    input: unpack(get_inputs_composite_subj_to_std_inverse)
    output:
        subj2std_invwarp = 'results/composite/sub-{subject}_to-{std_template}_via-cohort_CompositeInverseWarp.nii.gz'
    group: 'composite'
    shell: 'greedy -d 3 -rf {input.ref_subj} -r '
          ' {input.subj2cohort_affine_xfm_ras},-1 '
          ' {input.subj2cohort_invwarp}'
          ' {input.cohort2std_affine_xfm_ras},-1 '
          ' {input.cohort2std_invwarp} '
          ' -rc {output.subj2std_invwarp}'



