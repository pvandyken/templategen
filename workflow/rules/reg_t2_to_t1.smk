rule reg_aladin_t2w_to_t1:
    input: 
        ref = config['in_preproc_T1w'],
        flo = config['in_raw_T2w']
    output: 
        warped = bids(root='results',suffix='T2w.nii.gz',space='T1w',subject='{subject}'),
        xfm_ras = bids(root='results',suffix='xfm.txt',from_='T2w',to='T1w',type_='ras',subject='{subject}'),
    container: config['singularity']['prepdwi']
    shell:
        'reg_aladin -rigOnly -flo {input.flo} -ref {input.ref} -res {output.warped} -aff {output.xfm_ras}'

rule convert_xfm_ras2itk:
    input:
        xfm_ras = bids(root='results',suffix='xfm.txt',from_='T2w',to='T1w',type_='ras',subject='{subject}'),
    output:
        xfm_itk = bids(root='results',suffix='xfm.txt',from_='T2w',to='T1w',type_='itk',subject='{subject}'),
    container: config['singularity']['prepdwi']
    shell:
        'c3d_affine_tool {input.xfm_ras}  -oitk {output.xfm_itk}'

#run n4 on t2w image in t1w space, using t1w mask
rule n4biasfield_t2w:
    input: 
        img = bids(root='results',suffix='T2w.nii.gz',space='T1w',subject='{subject}'),
        mask = config['in_brainmask_T1w']
    output:
        img = bids(root='results',suffix='T2w.nii.gz',space='T1w',desc='n4',subject='{subject}'),
    threads: 8
    container: config['singularity']['prepdwi']
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'N4BiasFieldCorrection -d 3 -i {input.img} -x {input.mask} -o {output}'


rule mask_t1w:
    input: 
        mask = config['in_brainmask_T1w'],
        img = config['in_preproc_T1w']
    output:
        img = bids(root='results',suffix='T1w.nii.gz',desc='brain',subject='{subject}')
    container: config['singularity']['prepdwi']
    shell:
        'fslmaths {input.img} -mas {input.mask} {output.img}'

rule mask_t2w:
    input: 
        mask = config['in_brainmask_T1w'],
        img = bids(root='results',suffix='T2w.nii.gz',space='T1w',desc='n4',subject='{subject}'),
    output:
        img = bids(root='results',suffix='T2w.nii.gz',space='T1w',desc='brain',subject='{subject}')
    container: config['singularity']['prepdwi']
    shell:
        'fslmaths {input.img} -mas {input.mask} {output.img}'

    
