from pathlib import Path


out = "results/apoe4/"
config["hap_files"] = [
    x for x in Path(config["hap_files"]).glob("**/*")
    if x.is_file() and x.suffix == ".hap" and x.name[:-4] in config["haps"]
]
config["samples"] = [x.name[:-4] for x in config["hap_files"]]
config["hap_files"] = dict(zip(config["samples"], config["hap_files"]))
config["genotypes"] = Path(config["genotypes"])

rule all:
    input:
        expand(
            "h{heritability}/b{beta}/manhattan.pdf",
            beta=config["betas"], heritability=config["heritabilities"],
        )

rule transform:
    input:
        gts = config["hap_vars"],
        hap = lambda wildcards: str(config["hap_files"][wildcards.samp]),
    output:
        gts = out+"transform/{samp}.pgen"
    resources:
        runtime="0:01:00"
    log:
        out+"logs/transform/{samp}.log"
    benchmark:
        out+"bench/transform/{samp}.txt"
    conda:
        "../envs/haptools.yml"
    shell:
        "haptools transform -v INFO -o {output.gts} {input.gts} {input.hap} &> {log}"

rule merge:
    input:
        gts = config["hap_vars"],
        hps = rules.transform.output.gts,
    params:
        prefix = lambda w, output: Path(output.gts).with_suffix(""),
    output:
        log = temp(out+"merge/{samp}.log"),
        gts = out+"merge/{samp}.pgen",
    resources:
        runtime="0:01:00"
    log:
        out+"logs/merge/{samp}.log"
    benchmark:
        out+"bench/merge/{samp}.txt"
    conda:
        "../envs/haptools.yml"
    shell:
        "plink2 --pfile {input.gts} --pmerge {input.hps} --out {params.prefix} &> {log}"

rule sim_pts:
    input:
        gts = rules.transform.output.gts,
        hap = lambda wildcards: str(config["hap_files"][wildcards.samp]),
    params:
        beta = lambda wildcards: wildcards.beta,
        h2 = lambda wildcards: wildcards.heritability,
    output:
        pts = out+"h{heritability}/b{beta}/{samp}.pheno",
    resources:
        runtime="0:01:00"
    log:
        out+"logs/sim_pts/h{heritability}/b{beta}/{samp}.log"
    benchmark:
        out+"bench/sim_pts/h{heritability}/b{beta}/{samp}.txt"
    conda:
        "../envs/default.yml"
    shell:
        "haptools simphenotype -o {output.pts} -h {params.h2} -v DEBUG {input.gts} "
        "<( sed 's/EUR\\t0.99$/EUR\\t{params.beta}/' {input.hap} ) &>{log}"

rule gwas:
    input:
        gts = str(config["genotypes"]),
        pts = rules.sim_pts.output.pts,
    params:
        pgen_prefix = lambda w, input: Path(input.gts).with_suffix(""),
        out_prefix = lambda w, output: Path(output.log).with_suffix(""),
        name = lambda wildcards: wildcards.samp,
    output:
        log = temp(out+"h{heritability}/b{beta}/{samp}.log"),
        linear = out+"h{heritability}/b{beta}/{samp}.{samp}.glm.linear",
    resources:
        runtime="0:05:00"
    log:
        out+"logs/gwas/h{heritability}/b{beta}/{samp}.log"
    benchmark:
        out+"bench/gwas/h{heritability}/b{beta}/{samp}.txt"
    threads: 1
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --glm dominant --variance-standardize {params.name} "
        "--pheno {input.pts} --pfile {params.pgen_prefix} --out {params.out_prefix} "
        "--threads {threads} &>{log} || true"

rule manhattan:
    input:
        linear = expand(
            out+"h{heritability}/b{beta}/{samp}.{samp}.glm.linear",
            samp=config["samples"], allow_missing=True,
        ),
    params:
        linear = lambda wildcards, input: [f"-l {i}" for i in input.linear],
        ids = [f"-i {i}" for i in config["samples"][0].split("-")],
    output:
        png = out+"h{heritability}/b{beta}/manhattan.pdf",
    resources:
        runtime="0:02:00"
    log:
        out+"logs/manhattan/h{heritability}/b{beta}/manhattan.log"
    benchmark:
        out+"bench/manhattan/h{heritability}/b{beta}/manhattan.txt"
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/manhattan.py -o {output.png} {params.linear} {params.ids} "
        "&>{log}"
