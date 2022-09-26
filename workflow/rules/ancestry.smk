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
        gts = temp(out+"sim_gts/{samp}.vcf"),
        bkpt = out+"sim_gts/{samp}.bp",
    resources:
        runtime="2:30:00"
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
        "&> {log}"

rule vcf2pgen:
    input:
        gts = rules.sim_gts.output.gts,
    params:
        prefix = lambda w, output: Path(output.pgen).with_suffix(""),
    output:
        log = temp(out+"sim_gts/{samp}.log"),
        pgen = out+"sim_gts/{samp}.pgen",
        pvar = out+"sim_gts/{samp}.pvar",
        psam = out+"sim_gts/{samp}.psam",
    resources:
        runtime="0:02:00"
    log:
        out+"logs/vcf2pgen/{samp}.log"
    benchmark:
        out+"bench/vcf2pgen/{samp}.txt"
    threads: 12
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --vcf {input.gts} --make-pgen --snps-only --max-alleles 2 "
        "--threads {threads} --out {params.prefix} &> {log}"

rule transform:
    input:
        pgen = rules.vcf2pgen.output.pgen,
        pvar = rules.vcf2pgen.output.pvar,
        psam = rules.vcf2pgen.output.psam,
        bkpt = rules.sim_gts.output.bkpt,
        hap = config["hap"],
    params:
        ancs = lambda wildcards: ["", "--ancestry"][wildcards.type == "ancestry"],
        region = config["region"][0]+":"+"-".join(config["region"][1:]),
    output:
        pgen = temp(out+"transform/{type}/{samp}.pgen"),
        pvar = temp(out+"transform/{type}/{samp}.pvar"),
        psam = temp(out+"transform/{type}/{samp}.psam"),
    resources:
        runtime="0:00:30"
    log:
        out+"logs/transform/{type}/{samp}.log"
    benchmark:
        out+"bench/transform/{type}/{samp}.txt"
    conda:
        "../envs/haptools.yml"
    shell:
        "haptools transform -v INFO {params.ancs} -o {output.pgen} "
        "--region {params.region} {input.pgen} {input.hap} &> {log}"

rule sim_pts:
    input:
        pgen = rules.transform.output.pgen,
        pvar = rules.transform.output.pvar,
        psam = rules.transform.output.psam,
        hap = config["hap"],
    params:
        beta = lambda wildcards: wildcards.beta,
    output:
        pts = out+"b{beta}/{type}/{samp}.pheno",
    resources:
        runtime="0:00:30"
    log:
        out+"logs/sim_pts/b{beta}/{type}/{samp}.log"
    benchmark:
        out+"bench/sim_pts/b{beta}/{type}/{samp}.txt"
    conda:
        "../envs/haptools.yml"
    shell:
        "haptools simphenotype -o {output.pts} -v DEBUG {input.pgen} "
        "<( sed 's/YRI\\t0.99$/YRI\\t{params.beta}/' {input.hap} ) &>{log}"

rule merge:
    input:
        gts = rules.transform.input.pgen,
        gts_pvar = rules.transform.input.pvar,
        gts_psam = rules.transform.input.psam,
        hps = rules.transform.output.pgen,
        hps_pvar = rules.transform.output.pvar,
        hps_psam = rules.transform.output.psam,
    params:
        gts_prefix = lambda w, input: Path(input.gts).with_suffix(""),
        hps_prefix = lambda w, input: Path(input.hps).with_suffix(""),
        prefix = lambda w, output: Path(output.pgen).with_suffix(""),
    output:
        log = temp(out+"merge/{type}/{samp}.log"),
        pgen = out+"merge/{type}/{samp}.pgen",
        pvar = out+"merge/{type}/{samp}.pvar",
        psam = out+"merge/{type}/{samp}.psam",
    resources:
        runtime="0:10:00"
    log:
        out+"logs/merge/{type}/{samp}.log"
    benchmark:
        out+"bench/merge/{type}/{samp}.txt"
    threads: 1
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --pfile {params.gts_prefix} --pmerge {params.hps_prefix} "
        "--threads {threads} --out {params.prefix} &> {log}"

rule gwas:
    input:
        pgen = rules.merge.output.pgen,
        pvar = rules.merge.output.pvar,
        psam = rules.merge.output.psam,
        pts = rules.sim_pts.output.pts,
    params:
        in_prefix = lambda w, input: Path(input.pgen).with_suffix(""),
        out_prefix = lambda w, output: Path(output.log).with_suffix(""),
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
        "plink2 --linear allow-no-covars --variance-standardize "
        "--pheno {input.pts} --pfile {params.in_prefix} --out {params.out_prefix} "
        "--threads {threads} &>{log}"

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

# rule power:
#     input:
#         pts = 
