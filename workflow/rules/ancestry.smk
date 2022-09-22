from pathlib import Path


out = "ancestry"
# create a region param that encodes a 1 Mbp region around the given site
config["region"] = tuple(map(int, config["snp"].split(":")))
config["region"] = tuple(map(str, (
    config["region"][0],
    int(config["region"][1] - 1e6),
    int(config["region"][1] + 1e6),
)))


rule sim_gts:
    input:
        ref = config["reference"],
        samps = config["sample_info"],
        model = lambda wildcards: config["models"][wildcards.samp],
        mapdir = config["mapdir"],
    params:
        chroms = config["region"][0],
        region = config["region"][0]+":"+"-".join(config["region"][1:]),
        out_prefix = out+"sim_gts/gts",
        nsamps = 10*10000,
    output:
        gts = out+"sim_gts/{samp}.vcf",
        bkpt = out+"sim_gts/{samp}.bp",
    resources:
        runtime="2:00:00"
    log:
        out+"logs/sim_gts/{samp}.log"
    benchmark:
        out+"bench/sim_gts/{samp}.txt"
    conda:
        "../envs/default.yml"
    shell:
        "haptools simgenotype --invcf {input.ref} --sample_info {input.samps} "
        "--model {input.model} --mapdir {input.mapdir} --chroms {params.chroms} "
        "--out {params.out_prefix} --region {params.region} --popsize {params.nsamps}"

rule transform:
    input:
        gts = rules.sim_gts.output.gts,
        bkpt = rules.sim_gts.output.bkpt,
        hap = config["hap"],
    params:
        ancs = lambda wildcards: ["", "--ancestry"][wildcards.type == "ancestry"]
    output:
        gts = out+"transform/{type}/{samp}.vcf.gz",
    resources:
        runtime="1:00:00"
    log:
        out+"logs/transform/{type}/{samp}.log"
    benchmark:
        out+"bench/transform/{type}/{samp}.txt"
    conda:
        "../envs/default.yml"
    shell:
        "haptools transform {params.ancs} -o {output.gts} "
        "{input.gts} {input.hap}"

rule merge:
    input:
        gts = rules.transform.input.gts,
        hps = rules.transform.output.gts,
    output:
        gts = out+"merge/{type}/{samp}.vcf.gz",
    resources:
        runtime="0:30:00"
    log:
        out+"logs/merge/{type}/{samp}.log"
    benchmark:
        out+"bench/merge/{type}/{samp}.txt"
    conda:
        "../envs/default.yml"
    shell:
        "bcftools concat -Oz -o {output.gts} -a {input.gts} {input.hps}"

rule sim_pts:
    input:
        gts = rules.transform.output.gts,
        hap = config["hap"],
    params:
        beta = lambda wildcards: wildcards.beta,
    output:
        pts = out+"b{beta}/{type}/{samp}.pheno",
    resources:
        runtime="0:01:00"
    log:
        out+"logs/sim_pts/b{beta}/{type}/{samp}.log"
    benchmark:
        out+"bench/sim_pts/b{beta}/{type}/{samp}.txt"
    conda:
        "../envs/default.yml"
    shell:
        "haptools simphenotype -o {output.pts} -v DEBUG {input.gts} "
        "<( sed 's/YRI\\t0.99$/YRI\\t{params.beta}/' {input.hap} ) &>{log}"

rule gwas:
    input:
        gts = rules.merge.output.gts,
        pts = rules.sim_pts.output.pts,
    params:
        out_prefix = out+"/b{beta}/{type}/{samp}",
        name = lambda wildcards: wildcards.samp,
    output:
        log = temp(out+"/b{beta}/{type}/{samp}.log"),
        linear = out+"/b{beta}/{type}/{samp}.{samp}.glm.linear",
    resources:
        runtime="0:05:00"
    log:
        out+"/logs/gwas/b{beta}/{type}/{samp}.log"
    benchmark:
        out+"/bench/gwas/b{beta}/{type}/{samp}.txt"
    threads: 1
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --glm dominant --variance-standardize {params.name} "
        "--pheno {input.pts} --vcf {input.gts} --out {params.out_prefix} "
        "--threads {threads} &>{log} || true"

rule manhattan:
    input:
        linear = expand(
            out+"/b{beta}/{samp}.{samp}.glm.linear",
            samp=config["models"].keys(), allow_missing=True,
        ),
    params:
        linear = lambda wildcards, input: [f"-l {i}" for i in input.linear],
        ids = [f"-i {i}" for i in list(config["models"].keys())[0].split("-")],
    output:
        png = out+"/b{beta}/manhattan.pdf",
    resources:
        runtime="0:02:00"
    log:
        out+"/logs/manhattan/b{beta}/manhattan.log"
    benchmark:
        out+"/bench/manhattan/b{beta}/manhattan.txt"
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/manhattan.py -o {output.png} {params.linear} {params.ids} "
        "&>{log}"
