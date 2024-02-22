rule qc_get_template_edges:
    input:
        brain=expand(
            rules.apply_template_update.output.template,
            allow_missing=True,
            iteration=config["num_iters"],
            channel=channels[0],
        ),
        # mask=bids(
        #     root=root,
        #     datatype="anat",
        #     desc="brain",
        #     suffix="mask.nii.gz",
        #     **subj_wildcards,
        # ),
    output:
        temp(tempout("templateEdges.nii.gz")),
    shadow:
        "minimal"
    container:
        config["singularity"]["itksnap"]
    shell:
        "c3d {input} -canny 0mm 0 0.2 {output}"


rule qc_reg:
    input:
        overlay=rules.get_final_warped.output[0],
        lines=rules.qc_get_template_edges.output[0],
    output:
        png=bids(
            root=output_dir,
            datatype="qc",
            from_="individual",
            to="{template}",
            suffix="reg.png",
            **inputs.subj_wildcards,
        ),
    container:
        config["singularity"]["python"]
    script:
        "../scripts/qc_reg.py"


rule compile_qc_reg_manifest:
    input:
        reg_qc=inputs["FA"].expand(
            rules.qc_reg.output.png,
            template=config["template_name"],
        ),
    output:
        os.path.join(qc, "data", "reg.json"),
    run:
        with open(output[0], "w") as f:
            json.dump(
                {
                    "title": "Registration",
                    "images": sorted(
                        map(
                            lambda p: str(Path(p).relative_to(config["output_dir"])),
                            input["reg_qc"]
                        ),
                    ),
                },
                f,
            )

rule qc:
    input:
        reg_qc=rules.compile_qc_reg_manifest.output[0],
    output:
        os.path.join(qc, "data.json"),
    run:
        with open(output[0], "w") as f:
            json.dump(
                {
                    "reg": json.loads(Path(input["reg_qc"]).read_text()),
                },
                f,
            )


_qc_app = os.path.join(workflow.basedir, "..", "resources", "qc-app.tar.gz")


def _get_tar_contents(file):
    try:
        return [
            p
            for p in sp.check_output(["tar", "-tf", _qc_app]).decode().splitlines()
            if p[-1] != "/"
        ]
    except sp.CalledProcessError as err:
        raise Exception("Unable to find qc-app.tar.gz...") from err


rule unpack_qc_app:
    input:
        os.path.join(workflow.basedir, "..", "resources", "qc-app.tar.gz"),
    output:
        _get_tar_contents(_qc_app),
    shell:
        "tar -xvzf {input}"
