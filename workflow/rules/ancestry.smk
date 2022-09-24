from pathlib import Path


out = "results/ancestry/"
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
        out_prefix = lambda w, output: Path(output.gts).with_suffix(""),
        nsamps = 10*10000,
    output:
        gts = out+"sim_gts/{samp}.vcf.gz",
        bkpt = out+"sim_gts/{samp}-nochr.bp",
    resources:
        runtime="2:00:00"
    log:
        out+"logs/sim_gts/{samp}.log"
    benchmark:
        out+"bench/sim_gts/{samp}.txt"
    conda:
        "../envs/haptools.yml"
    shell:
        "haptools simgenotype --invcf {input.ref} --sample_info {input.samps} "
        "--model {input.model} --mapdir {input.mapdir} --chroms {params.chroms} "
        "--out {params.out_prefix} --region {params.region} --popsize {params.nsamps}"
        "&> {log} && mv {params.out_prefix}/{wildcards.samp}.bp {output.bkpt} "
        "&>> {log} && bgzip {params.out_prefix}.vcf &>> {log}"

rule index:
    input:
        gts = rules.sim_gts.output.gts,
    output:
        idx = out+"sim_gts/{samp}.vcf.gz.tbi",
    resources:
        runtime="0:01:00"
    log:
        out+"logs/index/{samp}.log"
    benchmark:
        out+"bench/index/{samp}.txt"
    conda:
        "../envs/default.yml"
    shell:
        "tabix -p vcf {input.gts} &>> {log}"

rule add_chr_to_bp:
    input:
        bp = rules.sim_gts.output.bkpt,
    output:
        bkpt = out+"sim_gts/{samp}.bp",
    resources:
        runtime="2:00:00"
    log:
        out+"logs/add_chr_to_bp/{samp}.log"
    benchmark:
        out+"bench/add_chr_to_bp/{samp}.txt"
    conda:
        "../envs/default.yml"
    shell:
        "awk -F $'\\t' -v 'OFS=\\t' "
        "'$1 !~ /^Sample.*/ {{print $1, \"chr\"$2, $3, $4; next}}1' "
        "{input.bp} > {output.bkpt} &> {log}"

rule transform:
    input:
        gts = rules.sim_gts.output.gts,
        bkpt = rules.add_chr_to_bp.output.bkpt,
        hap = config["hap"],
    params:
        ancs = lambda wildcards: ["", "--ancestry"][wildcards.type == "ancestry"],
        region = "chr"+config["region"][0]+":"+"-".join(config["region"][1:]),
    output:
        gts = out+"transform/{type}/{samp}.vcf.gz",
        idx = out+"transform/{type}/{samp}.vcf.gz.tbi",
    resources:
        runtime="1:00:00"
    log:
        out+"logs/transform/{type}/{samp}.log"
    benchmark:
        out+"bench/transform/{type}/{samp}.txt"
    conda:
        "../envs/haptools.yml"
    shell:
        "haptools transform -v INFO {params.ancs} -o {output.gts} "
        "--region {params.region} {input.gts} {input.hap} &> {log} && "
        "tabix -p vcf {output.gts} &>> {log}"

rule merge:
    input:
        gts = rules.transform.input.gts,
        hps = rules.transform.output.gts,
        gts_idx = rules.index.output.idx,
        hps_idx = rules.transform.output.idx,
    params:
        region = "chr"+config["region"][0]+":"+"-".join(config["region"][1:]),
    output:
        gts = out+"merge/{type}/{samp}.vcf.gz",
    resources:
        runtime="0:10:00"
    log:
        out+"logs/merge/{type}/{samp}.log"
    benchmark:
        out+"bench/merge/{type}/{samp}.txt"
    conda:
        "../envs/default.yml"
    shell:
        "bcftools concat -Oz -o {output.gts} -r {params.region} -a "
        "{input.gts} {input.hps} &> {log}"

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
        out_prefix = lambda w, output: Path(output.log).with_suffix(""),
        name = lambda wildcards: wildcards.samp,
    output:
        log = temp(out+"b{beta}/{type}/{samp}.log"),
        linear = out+"b{beta}/{type}/{samp}.{samp}.glm.linear",
    resources:
        runtime="0:05:00"
    log:
        out+"logs/gwas/b{beta}/{type}/{samp}.log"
    benchmark:
        out+"bench/gwas/b{beta}/{type}/{samp}.txt"
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
            out+"b{beta}/{samp}.{samp}.glm.linear",
            samp=config["models"].keys(), allow_missing=True,
        ),
    params:
        linear = lambda wildcards, input: [f"-l {i}" for i in input.linear],
        ids = [f"-i {i}" for i in list(config["models"].keys())[0].split("-")],
    output:
        png = out+"b{beta}/manhattan.pdf",
    resources:
        runtime="0:02:00"
    log:
        out+"logs/manhattan/b{beta}/manhattan.log"
    benchmark:
        out+"bench/manhattan/b{beta}/manhattan.txt"
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/manhattan.py -o {output.png} {params.linear} {params.ids} "
        "&>{log}"
